#!/usr/bin/env python
# --------------------------------------------------------
#       adds clustering and charge to trees created with pXar
# created on August 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from dataclasses import field

import uproot
from numpy import all, ones, count_nonzero  # noqa
from uproot import ReadOnlyDirectory
from uproot.models import TTree

from src.proteus import Proteus
from src.run import Run, Analysis
from utility.utils import *
from plotting.utils import download_file, remove_file


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
        self.NDUTPlanes = self.Run.NDUTs

        # PRE-CONVERTER
        self.Raw = self.init_raw()
        self.Proteus = self.init_proteus()

        # FILES
        self.OutFilePath = self.SaveDir.joinpath(f'run{self.Run:04d}.hdf5')
        self.RawFile: ReadOnlyDirectory = field(init=False)  # output from proteus
        self.F = None

    def __repr__(self):
        return f'{self.__class__.__name__} of {self.Run!r}'

    def run(self, force=False, steps=None, rm=True):
        if force or not self.OutFilePath.exists():
            t0 = info(f'Starting {self!r}\n')
            for i, (s, f) in enumerate(choose(steps, self.steps)):
                if not f.exists():
                    print_banner(f'Start converter step {i}: {s.__doc__}')
                    s()
                else:
                    info(f'found outfile of step {i}, continue with next step ({f})')
            if rm:
                self.remove_aux_files()
            add_to_info(t0, f'\nFinished {self!r} in ', color=GREEN)
            return True
        return False

    def raw2root(self):
        self.run(force=True, steps=self.first_steps, rm=False)

    def realign(self):
        self.raw2root()
        self.Proteus.align(force=True)
        self.remove_aux_files()

    @property
    def first_steps(self):
        return self.Raw.Steps

    @property
    def steps(self):
        return self.first_steps + self.Proteus.Steps + [(self.root2hdf5, self.OutFilePath)]

    @property
    def aux_files(self):
        return [self.Proteus.RawFilePath, self.Proteus.OutFilePath, self.Proteus.HistFilePath]

    @property
    def raw_files(self):
        return [self.Raw.RawFilePath]

    def remove_aux_files(self):
        remove_file(*self.aux_files)

    @staticmethod
    def download_raw_file(f: Path, out=True):
        server, loc = [Analysis.Config.get('data', n) for n in ['server', 'server dir']]
        info(f'downloading DUT raw file from {server}:{loc}')
        f.parent.mkdir(exist_ok=True)
        download_file(server, Path(loc).expanduser().joinpath(*f.parts[-4:]), f, out)

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
        conf_dir = Dir.joinpath('proteus', self.DataDir.stem)
        me, se = [Analysis.Config.getint('align', opt) for opt in ['max events', 'skip events']]
        return Proteus(soft_dir, data_dir, conf_dir, raw_file=self.proteus_raw_file_path(), max_events=me, skip_events=se, dut_pos=self.Run.Positions, duts=self.Run.DUTs)

    def init_raw(self):
        from src.raw import Raw
        return Raw(self)

    def init_event_alignment(self):
        return None

    @property
    def calibration(self):
        from src.calibration import Calibration
        return Calibration

    def load_calibration(self, dut_nr=None):
        return self.calibration(self.Run if dut_nr is None else Run(self.Run.Number, dut_nr, self.Run.TCDir, single_mode=True))

    def proteus_raw_file_path(self):
        return self.Raw.OutFilePath

    def load_file(self, opt='w'):
        return h5py.File(self.OutFilePath, opt)
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region HDF5
    def root2hdf5(self):
        """ convert tracked root file to hdf5 file. """
        remove_file(self.OutFilePath)  # remove hdf5 file if it exists
        start_time = info('Start root -> hdf5 conversion ...')

        self.F = self.load_file()
        self.RawFile = uproot.open(self.Proteus.OutFilePath)

        self.add_tracks()
        self.add_planes()

        add_to_info(start_time, '\nFinished hdf5 conversion in')
        self.F.close()

    def add_tracks(self):
        t0 = info('add track information ...', endl=False)
        g = self.F.create_group('Tracks')
        b = array([['evt_frame', 'Events', 'u4'], ['evt_ntracks', 'N', 'u1'],
                   ['trk_size',  'Size',   'f2'],
                   ['trk_chi2',  'Chi2',   'f2'], ['trk_dof', 'Dof', 'u1'],
                   ['trk_du',    'SlopeX', 'f2'], ['trk_dv', 'SlopeY', 'f2']]).T
        tree = self.RawFile['C0/tracks_clusters_matched']
        self.add_time_stamp(tree)
        self.add_data(tree, g, b)
        add_to_info(t0, color=GREEN)

    @staticmethod
    def get_time_stamp(tree: TTree):
        t = array(tree['evt_timestamp'])
        return ((t - t[0]) / 1e9).astype('f4')  # time stamp is just a number counting up

    def add_time_stamp(self, tree: TTree):
        self.F.create_group('Event').create_dataset('Time', data=self.get_time_stamp(tree))

    def add_planes(self):
        n = len(self.RawFile.keys(recursive=False))  # noqa (raises stupid warning)
        info(f'add {n} planes ... ')
        PBAR.start(n * 2)
        for pl in range(n):
            self.add_plane(pl)

    def add_plane(self, i):
        g = self.F.create_group(f'Plane{i}')
        key = list(self.RawFile.keys(recursive=False))[i]

        # mask
        m_tree = self.RawFile[f'{key}/masked_pixels']
        g.create_dataset('Mask', data=array(list(m_tree.arrays(library='np').values()), 'u2'))

        # cluster & track interpolations
        tree = self.RawFile[f'{key}/tracks_clusters_matched']
        self.add_plane_tracks(g, tree)
        self.add_clusters(g, tree)
        if i >= self.NTelPlanes + self.Proteus.NRefPlanes:
            self.add_trigger_info(g)

    @update_pbar
    def add_plane_tracks(self, group, tree: TTree):
        g = group.create_group('Tracks')
        b = array([['trk_u', 'U', 'f2'], ['trk_v', 'V', 'f2'], ['trk_col', 'X', 'f2'], ['trk_row', 'Y', 'f2'], ['trk_std_u', 'eU', 'f2'], ['trk_std_v', 'eV', 'f2']]).T
        self.add_data(tree, g, b)

    @update_pbar
    def add_clusters(self, group, tree: TTree):
        g = group.create_group('Clusters')
        cs = self.add_data(tree, g, array([['clu_size', 'Size', 'u2'], ['evt_nclusters', 'N', 'u2']]).T)['clu_size']
        b = array([['clu_u', 'U', 'f2'], ['clu_v', 'V', 'f2'], ['clu_col', 'X', 'f2'], ['clu_row', 'Y', 'f2'], ['clu_value', 'Charge', 'i']]).T
        self.add_data(tree, g, b, cut=cs > 0)  # filter out the nan events

    def trigger_info_file(self):
        return self.Raw.OutFilePath

    def add_trigger_info(self, group):
        tree = uproot.open(self.trigger_info_file())[f'{group.name}/Hits']
        g = group.create_group('Trigger')
        b = array([[key, key.replace('Trigger', ''), 'u1'] for key in tree.keys(filter_name='Trigger*')]).T
        self.add_data(tree, g, b)

    @staticmethod
    def add_data(tree: TTree, g, b, cut=...):
        data = tree.arrays(b[0], library='np')
        for i, (old_name, new_name, dtype) in enumerate(b.T):
            g.create_dataset(new_name, data=data[old_name][cut].astype(dtype))
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
    from cern.converter import CERNConverter
    from src.analysis import Analysis
    import awkward as aw  # noqa

    parser = ArgumentParser()
    parser.add_argument('run', nargs='?', default=232)
    parser.add_argument('dut', nargs='?', default=0)
    parser.add_argument('-tc', nargs='?', default=None)
    pargs = parser.parse_args()

    ana_ = Analysis(pargs.tc)
    tc = ana_.BeamTest

    z = (Converter if tc.Location == 'DESY' else CERNConverter).from_ana(pargs.run, pargs.dut, ana_, single_mode=True)
    r = z.Raw
    p = z.Proteus
    c = z.load_calibration()
    rn = z.Run
    draw = c.Draw
    if hasattr(z, 'Adc2Vcal'):
        adc = z.Adc2Vcal
