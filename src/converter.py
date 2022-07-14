#!/usr/bin/env python
# --------------------------------------------------------
#       adds clustering and charge to trees created with pXar
# created on August 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from numpy import sum, append, delete, average, all, ones, count_nonzero

from src.proteus import Proteus
from utility.utils import *
from src.desy_run import DESYRun
from src.calibration import Calibration


class Converter:

    def __init__(self, data_dir: Path, run_number, config):
        """ Converts EUDAQ2 raw files in several steps into hdf5 files.
            STEP 0: raw -> root             (EUDAQ2)
            STEP 1: noisescan               (proteus)
            STEP 2: alignment               (proteus)
            STEP 3: track reconstruction    (proteus)
            STEP 4: root -> hdf5            (python) """

        self.RunNumber = run_number
        self.Config = config

        # DIRECTORIES
        self.DataDir = data_dir
        self.SaveDir = data_dir.joinpath('data')
        self.SoftDir = Path(self.Config.get('SOFTWARE', 'dir')).expanduser()

        self.NTelPlanes = self.Config.getint('TELESCOPE', 'planes')
        self.NDUTPlanes = self.Config.getint('DUT', 'planes')

        # PRE-CONVERTER
        self.Raw = self.init_raw()
        self.Proteus = self.init_proteus()

        # FILES
        self.OutFileName = self.SaveDir.joinpath(f'run{self.RunNumber:04d}.hdf5')
        self.TrackFile = None
        self.F = None

        self.Steps = self.Raw.Steps + self.Proteus.Steps + [(self.root_2_hdf5, self.OutFileName)]
        self.PBar = PBar()

    def run(self):
        for i, (s, f) in enumerate(self.Steps):
            if not f.exists():
                print_banner(f'Start converter step {i}: {s.__doc__}')
                s()

    # ----------------------------------------
    # region INIT
    def init_proteus(self):
        soft_dir = self.SoftDir.joinpath(self.Config.get('SOFTWARE', 'proteus'))
        data_dir = self.DataDir.joinpath('proteus')
        conf_dir = Dir.joinpath('proteus')
        return Proteus(soft_dir, data_dir, conf_dir, self.Raw.OutFilePath, *[self.Config.getint('align', opt) for opt in ['max events', 'skip events']])

    def init_raw(self):
        from src.raw import Raw
        return Raw(self)

    def get_calibration(self, dut_number=0):
        return Calibration(DESYRun(self.RunNumber, dut_number, self.DataDir, self.Config))
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region HDF5
    def root_2_hdf5(self):
        """ convert tracked root file to hdf5 file. """
        remove_file(self.OutFileName)  # remove hdf5 file if it exists
        start_time = info('Start root->hdf5 conversion ...')

        self.F = h5py.File(self.OutFileName, 'w')
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
        self.PBar.start(n * 2)
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
        hf = h5py.File(self.OutFileName)
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
        for s, f in self.Steps:
            if f.suffix == '.root' or f.suffix == '.hdf5' or all_:
                remove_file(f)

    @staticmethod
    def clusterise(hits):
        if not hits.size:
            return []
        # sort hits into clusters
        clusters = [Cluster(hits[0])]
        hits = delete(hits, 0, axis=0)
        while hits.size:
            in_existing_cluster = False
            n_deleted = 0
            for i, hit in enumerate(hits):
                for cluster in clusters:
                    if cluster.hit_is_adjacent(hit):
                        cluster.add_hit(hit)
                        hits = delete(hits, i - n_deleted, axis=0)
                        n_deleted += 1
                        in_existing_cluster = True
                        break
            if not in_existing_cluster:  # make a new cluster if none of the hits is adjacent to any of the existing Clusters
                clusters.append(Cluster(hits[0]))
                hits = delete(hits, 0, axis=0)
        return clusters
    # endregion MISC
    # ----------------------------------------


class Cluster:

    def __init__(self, seed_hit):
        self.Hits = array([seed_hit])

    def add_hit(self, hit):
        self.Hits = append(self.Hits, [hit], axis=0)

    def hit_is_adjacent(self, hit):
        """ returns: if any of the existing cluster hits has a distance of 1 to the given hit"""
        return any(sqrt(sum(abs(hit[:2] - self.Hits[:, :2]), axis=1)) <= 1.1)

    def get_charge(self):
        return sum(self.Hits[:, 2])

    def get_size(self):
        return self.Hits.shape[0]

    def get_x(self):
        weights = self.Hits[:, 2]
        return average(self.Hits[:, 0], weights=weights) if weights.nonzero()[0].size else mean(self.Hits[:, 0])

    def get_y(self):
        weights = self.Hits[:, 2]
        return average(self.Hits[:, 1], weights=weights) if weights.nonzero()[0].size else mean(self.Hits[:, 1])

    def __repr__(self):
        return 'Cluster with {} hits at [{:.1f}, {:.1f}] with charge {:.1f} vcals'.format(self.get_size(), self.get_x(), self.get_y(), self.get_charge())


if __name__ == '__main__':
    from argparse import ArgumentParser
    from analysis import Analysis

    parser = ArgumentParser()
    parser.add_argument('run', nargs='?', default=11)
    pargs = parser.parse_args()
    a = Analysis()

    z = Converter(a.BeamTest.Path, pargs.run, a.Config)
    r = z.Raw
    p = z.Proteus
    c = z.get_calibration()
