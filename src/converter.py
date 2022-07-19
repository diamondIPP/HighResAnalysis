#!/usr/bin/env python
# --------------------------------------------------------
#       adds clustering and charge to trees created with pXar
# created on August 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from numpy import all, ones, count_nonzero

from src.proteus import Proteus
from utility.utils import *
from src.run import Run, Analysis
from src.calibration import Calibration


class Converter:
    """
    Converts EUDAQ2 raw files in several steps into hdf5 files.
    STEP  0: raw -> root                  (EUDAQ2)
    STEP  1: noisescan                    (proteus)
    STEP  2: alignment                    (proteus)
    STEP  3: track reconstruction         (proteus)
    STEP  4: root -> hdf5                 (python)
    """

    def __init__(self, data_dir: Path, run_number):

        self.Run = Run(run_number, 0, data_dir, single_mode=True)

        # DIRECTORIES
        self.DataDir = data_dir
        self.SaveDir = data_dir.joinpath('data')
        self.SoftDir = Path(Analysis.Config.get('SOFTWARE', 'dir')).expanduser()

        self.NTelPlanes = Analysis.Config.getint('TELESCOPE', 'planes')
        self.NDUTPlanes = Analysis.Config.getint('DUT', 'planes')

        # PRE-CONVERTER
        self.Raw = self.init_raw()
        self.RawFilePath = self.load_raw_file_path()
        self.Proteus = self.init_proteus()

        # FILES
        self.OutFilePath = self.SaveDir.joinpath(f'run{self.Run:04d}.hdf5')
        self.TrackFile = None
        self.F = None

    def __repr__(self):
        return f'{self.__class__.__name__} of {self.Run!r}'

    def run(self):
        for i, (s, f) in enumerate(self.steps):
            if not f.exists():
                print_banner(f'Start converter step {i}: {s.__doc__}')
                s()
            else:
                info(f'found out file of step {i}, continue with next step ({f})')

    @property
    def first_steps(self):
        return self.Raw.Steps

    @property
    def steps(self):
        return self.first_steps + self.Proteus.Steps + [(self.root_2_hdf5, self.OutFilePath)]

    # ----------------------------------------
    # region INIT
    @classmethod
    def from_run(cls, run: Run):
        cal = cls(run.TCDir, run.Number)
        cal.Run = run
        return cal

    @classmethod
    def from_ana(cls, run_number, dut=0, ana: Analysis = None, single_mode=False):
        ana = choose(ana, Analysis)
        return cls.from_run(Run.from_ana(run_number, dut, ana, single_mode))

    def init_proteus(self):
        soft_dir = self.SoftDir.joinpath(Analysis.Config.get('SOFTWARE', 'proteus'))
        data_dir = self.DataDir.joinpath('proteus')
        conf_dir = Dir.joinpath('proteus')
        return Proteus(soft_dir, data_dir, conf_dir, self.RawFilePath, *[Analysis.Config.getint('align', opt) for opt in ['max events', 'skip events']])

    def init_raw(self):
        from src.raw import Raw
        return Raw(self)

    def load_calibration(self, dut_nr=None):
        return Calibration(self.Run if dut_nr is None else Run(self.Run.Number, dut_nr, self.Run.TCDir, single_mode=True))

    def load_raw_file_path(self):
        return self.Raw.OutFilePath
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region HDF5
    def root_2_hdf5(self):
        """ convert tracked root file to hdf5 file. """
        remove_file(self.OutFilePath)  # remove hdf5 file if it exists
        start_time = info('Start root->hdf5 conversion ...')

        self.F = h5py.File(self.OutFilePath, 'w')
        self.TrackFile = TFile(str(self.Proteus.OutFilePath))

        self.add_tracks()
        self.add_planes()

        add_to_info(start_time, '\nFinished conversion in')

    def add_tracks(self):
        info('add track information ...')
        g = self.F.create_group('Tracks')
        b = array([['evt_frame', 'Events', 'u4'], ['evt_ntracks', 'N', 'u1'],
                   ['trk_size', 'Size', 'f2'],
                   ['trk_chi2', 'Chi2', 'f2'], ['trk_dof', 'Dof', 'u1'],
                   ['trk_du', 'SlopeX', 'f2'], ['trk_dv', 'SlopeY', 'f2']]).T
        tree = self.TrackFile.Get('C0').Get('tracks_clusters_matched')  # same for all plane dirs
        tree.SetEstimate(-1)
        self.add_time_stamp(tree)
        self.add_data(tree, g, b)

    def add_time_stamp(self, tree):
        t = get_tree_vec(tree, var='evt_timestamp', dtype='u8')
        self.F['Tracks'].create_dataset('Time', data=((t - t[0]) / 1e9).astype('f4'))  # time stamp is just a number counting up

    def add_planes(self):
        n = self.NTelPlanes + self.NDUTPlanes
        info(f'add {self.NTelPlanes + self.NDUTPlanes} planes ... ')
        PBAR.start(n * 2)
        for pl in range(n):
            self.add_plane(pl)

    def add_plane(self, i):
        g = self.F.create_group(f'Plane{i}')
        key = self.TrackFile.GetListOfKeys()[i].GetName()

        # mask
        m_tree = self.TrackFile.Get(key).Get('masked_pixels')
        g.create_dataset('Mask', data=get_tree_vec(m_tree, ['col', 'row'], dtype='u2'))

        # cluster & track interpolations
        tree = self.TrackFile.Get(key).Get('tracks_clusters_matched')
        tree.SetEstimate(-1)
        self.add_plane_tracks(g, tree)
        self.add_clusters(g, tree)
        if i >= self.NTelPlanes:
            self.add_trigger_info(g)

    @update_pbar
    def add_plane_tracks(self, group, tree):
        g = group.create_group('Tracks')
        b = array([['trk_u', 'U', 'f2'], ['trk_v', 'V', 'f2'], ['trk_col', 'X', 'f2'], ['trk_row', 'Y', 'f2'], ['trk_std_u', 'eU', 'f2'], ['trk_std_v', 'eV', 'f2']]).T
        self.add_data(tree, g, b)

    @update_pbar
    def add_clusters(self, group, tree):
        g = group.create_group('Clusters')
        cs = self.add_data(tree, g, array([['clu_size', 'Size', 'u2'], ['evt_nclusters', 'N', 'u2']]).T)[0]
        b = array([['clu_u', 'U', 'f2'], ['clu_v', 'V', 'f2'], ['clu_col', 'X', 'f2'], ['clu_row', 'Y', 'f2']]).T
        self.add_data(tree, g, b, cut=cs > 0)  # filter out the nan events

    def add_trigger_info(self, group):
        f = TFile(str(self.Raw.OutFilePath))
        t = f.Get(group.name.strip('/')).Get('Hits')
        t.SetEstimate(-1)
        g = group.create_group('Trigger')
        self.add_data(t, g, array([['TriggerPhase', 'Phase', 'u1'], ['TriggerCount', 'Count', 'u1']]).T)

    @staticmethod
    def add_data(tree, g, b, cut=...):
        data = get_tree_vec(tree, b[0], dtype=b[2])
        for i, n in enumerate(b[1]):
            g.create_dataset(n, data=data[i][cut])
        return data
    # endregion HDF5
    # ----------------------------------------

    # ----------------------------------------
    # region ALIGN
    def align_tree(self):
        import uproot
        hf = h5py.File(self.OutFilePath)
        with uproot.recreate('bla.root') as f:
            n = [array(hf[f'Plane{pl}']['Clusters']['Size']) > 0 for pl in range(self.NTelPlanes)]
            cut = all(n, axis=0) & (array(hf['Tracks']['Size']) == self.NTelPlanes)
            for pl in range(self.NTelPlanes):
                d = f.mkdir(f'Plane{pl}')
                ci = cut[n[pl]]
                i = count_nonzero(ci)
                d['Clusters'] = {'NClusters': ones(i, 'i'),
                                 'Col': array(hf[f'Plane{pl}']['Clusters']['X'])[ci].astype('d').reshape((-1, 1)),
                                 'Row': array(hf[f'Plane{pl}']['Clusters']['Y'])[ci].astype('d').reshape((-1, 1)),
                                 'VarCol': ones((i, 1)),
                                 'VarRow': ones((i, 1)),
                                 'Value': ones((i, 1)),
                                 'CovColRow': zeros((i, 1)),
                                 'Timing': zeros((i, 1)),
                                 'Track': ones((i, 1), 'i')}
            f['Plane0']['Clusters'].show()
    # endregion ALIGN
    # ----------------------------------------

    # ----------------------------------------
    # region MISC
    @staticmethod
    def check_root_version():
        v = gROOT.GetVersion()
        return True if v.startswith('6') else critical(f'ROOT 6 required for the conversion! Current version: {v}')

    def remove_files(self, all_=False):
        for s, f in self.steps:
            if f.suffix == '.root' or f.suffix == '.hdf5' or all_:
                remove_file(f)
    # endregion MISC
    # ----------------------------------------


if __name__ == '__main__':
    from argparse import ArgumentParser

    parser = ArgumentParser()
    parser.add_argument('run', nargs='?', default=41)
    parser.add_argument('dut', nargs='?', default=1)
    pargs = parser.parse_args()

    z = Converter.from_ana(pargs.run, pargs.dut)
    r = z.Raw
    p = z.Proteus
    c = z.load_calibration()
