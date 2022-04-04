#!/usr/bin/env python
# --------------------------------------------------------
#       convert EUDAQ-2 raw files to hdf5
# created on August 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from argparse import ArgumentParser
from glob import glob
from os import chdir
from subprocess import check_call

from numpy import cumsum, split, sum, nan, min, invert

from plotting.draw import *
from src.calibration import Calibration
from src.converter import Converter
from src.desy_run import DESYRun
from src.dut import Plane
from src.utils import print_banner, file_exists, basename, h5py, get_tree_vec, update_pbar
from ROOT import TFile
import toml


class DESYConverter(Converter):
    """ Converts EUDAQ2 raw files in several steps into hdf5 files.
        STEP 0: raw -> root     (EUDAQ2)
        STEP 1: noisescan       (proteus)
        STEP 2: alignment       (proteus)
        STEP 3: tracking        (proteus)
        STEP 4: track-matching  (proteus)
        STEP 5: root -> hdf5    (python) """

    def __init__(self, data_dir, run_number, config):
        Converter.__init__(self, data_dir, run_number, config)

        # DIRECTORIES
        self.EUDAQDir = self.SoftDir.joinpath(self.Config.get('SOFTWARE', 'eudaq2'))
        self.ProteusSoftDir = self.SoftDir.joinpath(self.Config.get('SOFTWARE', 'proteus'))
        self.ProteusDataDir = self.DataDir.joinpath('proteus')

        # FILENAMES
        self.ROOTFileName = join(self.SaveDir, 'run{:06d}.root'.format(self.RunNumber))
        self.AlignSteps = self.get_align_steps()
        self.TrackFileName = join(self.ProteusDataDir, 'track-{:04d}-data.root'.format(self.RunNumber))
        self.MatchFileName = join(self.ProteusDataDir, 'match-{:04d}-trees.root'.format(self.RunNumber))
        self.FinalFileName = join(self.SaveDir, 'run{:04d}.hdf5'.format(self.RunNumber))
        self.FileNames = [self.ROOTFileName] + self.get_align_files() + [self.TrackFileName, self.MatchFileName, self.FinalFileName]

        self.NSteps = len(self.FileNames)
        self.NTelPlanes = self.Config.getint('TELESCOPE', 'planes')
        self.NDUTPlanes = self.Config.getint('DUT', 'planes')

        # FILES
        self.F = None
        self.TrackFile = TFile()
        self.MatchFile = TFile()

        self.Draw = Draw()

    # TODO: add check for alignment or align every run?
    # TODO: remove root files after conversion

    def run(self):
        for f in self.get_steps():
            print_banner('Starting{}'.format(f.__doc__))
            f()

    # ----------------------------------------
    # region INIT
    def load_raw_file_name(self):
        names = glob(join(self.DataDir, 'raw', 'run{:06d}*.raw'.format(self.RunNumber)))
        return names[0] if names else None

    def get_steps(self):
        if file_exists(self.FinalFileName):
            return []
        steps = [self.convert_raw_to_root, self.noise_scan, self.align, self.track, self.match, self.convert_root_to_hdf5]
        final_steps = []
        for step, f in zip(steps, self.FileNames):
            if not file_exists(f):
                final_steps.append(step)
        return final_steps

    @staticmethod
    def check_root_version():
        version = gROOT.GetVersion()
        if not version.startswith('6'):
            warning(f'ROOT6 required proteus! Current version: {version}')
        return version.startswith('6')
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_z_positions(self, raw=False):
        d = toml.load(self.ProteusDataDir.joinpath('geometry.toml' if raw else self.make_toml_name()))
        return array([s['offset'][-1] for s in d['sensors']])

    def get_align_steps(self):
        with open(join(self.ProteusDataDir, 'analysis.toml')) as f:
            return [line.split(']')[0].replace('[align.', '') for line in f.readlines() if line.startswith('[align.')]

    def get_align_files(self):
        return [join(self.ProteusDataDir, name) for name in [self.make_toml_name('all', 'mask', 'mask'), self.make_toml_name()]]

    def get_alignment(self, step=-1):
        return toml.load(self.ProteusDataDir.joinpath(self.make_toml_name(self.AlignSteps[step])))

    def get_calibration(self, dut_number):
        return Calibration(DESYRun(self.RunNumber, dut_number, self.DataDir, self.Config))
    # endregion GET
    # ----------------------------------------

    def convert_raw_to_root(self):
        """ step 0: convert raw file to root file with eudaq"""
        if self.RawFilePath is None:
            critical('raw file does not exist for run: {}'.format(self.RunNumber))
        ensure_dir(self.SaveDir)
        chdir(self.SaveDir)
        cmd = '{} -i {} -o {}'.format(join(self.EUDAQDir, 'bin', 'euCliConverter'), self.RawFilePath, self.ROOTFileName)
        info('Convert {} to {}.root using eudaq-2\n'.format(basename(self.RawFilePath), self.ROOTFileName))
        check_call(cmd.split())
        for filename in glob('AutoDict_vector*'):
            remove(filename)

    # ----------------------------------------
    # region HDF5 Converter
    def convert_root_to_hdf5(self):
        """ step 5: convert root file to hdf5 file. """
        if file_exists(self.FinalFileName):
            remove_file(self.FinalFileName)
        start_time = info('Start root->hdf5 conversion ...')
        # files
        self.F = h5py.File(self.FinalFileName, 'w')
        self.TrackFile = TFile(self.TrackFileName)
        self.MatchFile = TFile(self.MatchFileName)

        self.add_track_vars()
        self.add_time_stamp()
        self.add_tel_planes()
        self.add_dut_planes()

        add_to_info(start_time, '\nFinished conversion in')

    def add_track_vars(self):
        info('add track information ...')
        g = self.F.create_group('Tracks')

        track_types = {'Chi2': 'f2', 'Dof': 'u1', 'X': 'f2', 'Y': 'f2', 'SlopeX': 'f2', 'SlopeY': 'f2'}
        match_types = {'evt_frame': 'u4', 'evt_ntracks': 'u1', 'trk_size': 'f2'}
        self.PBar.start(len(track_types) + len(match_types), counter=True, t='s')

        # from tracking tree
        t_tree = self.TrackFile.Get('Tracks')
        t_tree.SetEstimate(t_tree.GetEntries())
        g.create_dataset('NTracks', data=get_tree_vec(t_tree, var='NTracks', dtype='u1'))
        t_tree.SetEstimate(sum(g['NTracks']))
        for n, t in track_types.items():
            g.create_dataset(n, data=get_tree_vec(t_tree, n, dtype=t))
            self.PBar.update()

        # from match tree
        m_tree = self.MatchFile.Get('C0').Get('tracks_clusters_matched')  # same for all plane dirs
        m_tree.SetEstimate(m_tree.GetEntries())
        for n, t in match_types.items():
            g.create_dataset(n.replace('trk_', '').title().replace('_', ''), data=get_tree_vec(m_tree, n, dtype=t))
            self.PBar.update()

    def add_time_stamp(self):
        tree = self.TrackFile.Get('Event')
        tree.SetEstimate(tree.GetEntries())
        t_vec = get_tree_vec(self.TrackFile.Get('Event'), var='TriggerTime', dtype='u8')
        g = self.F.create_group('Event')
        g.create_dataset('Time', data=((t_vec - t_vec[0]) / 1e9).astype('f4'))  # convert time vec to seconds

    def add_plane(self, i):
        group = self.F.create_group(f'Plane{i}')
        i_dut = i - self.NTelPlanes
        pl_tag = f'C{i_dut}' if i_dut >= 0 else f'M{i}'

        # mask
        m_tree = self.MatchFile.Get(pl_tag).Get('masked_pixels')
        m_tree.SetEstimate(m_tree.GetEntries())
        group.create_dataset('Mask', data=get_tree_vec(m_tree, ['col', 'row'], dtype='u2'))

        # cluster
        tree = self.MatchFile.Get(pl_tag).Get('tracks_clusters_matched')
        tree.SetEstimate(-1)
        self.add_tracks(group, tree)
        self.add_clusters(group, tree)
        if i_dut >= 0:
            self.add_cluster_charge(tree, group, i_dut)
            self.remove_doubles(group, tree)

            # trigger info
            self.add_trigger_info(group)

    def add_tel_planes(self):
        info(f'reading branches for {self.NTelPlanes} telescope planes ... ')
        self.PBar.start(self.NTelPlanes * 2)
        for i in range(self.NTelPlanes):
            self.add_plane(i)

    def add_dut_planes(self):
        info(f'reading branches for {self.NDUTPlanes} DUT planes ... ')
        for i in arange(self.NDUTPlanes) + self.NTelPlanes:
            self.add_plane(i)

    @staticmethod
    def remove_doubles(group, tree, only_1cluster=False):
        """ find the tracks which are closer to the cluster in case there are more than 1 track per cluster. """
        s = array(group['Clusters']['Size'])
        cut = s > 0
        nt = get_tree_vec(tree, var='evt_ntracks', dtype='u1')[cut]      # number of tracks
        nc = get_tree_vec(tree, var='evt_nclusters', dtype='u1')[cut]    # number of clusters
        tu, tv = array(group['Tracks']['U'])[cut], array(group['Tracks']['V'])[cut]
        u, v = array(group['Clusters']['U']), array(group['Clusters']['V'])
        r = sqrt((tu - u) ** 2 + (tv - v) ** 2)
        good = nt == 1
        c2 = (nt > 1) & (nc == 2)  # selector for one cluster for two tracks, if nc == 2 and nt > 1 there is always a double
        r2 = r[c2][::2] < r[c2][1::2]
        r2c = concatenate(array([r2, invert(r2)]).T)  # merge with zipper method
        good[c2] = r2c
        # more than 1 cluster
        if not only_1cluster:
            for nt_i in arange(2, 5):
                for nc_i in arange(4, 7, 2):
                    c = (nt == nt_i) & (nc == nc_i)
                    ri = array(split(r[c], arange(nt_i, r[c].size, nt_i)))
                    good[c] = concatenate(ri) == min(ri, axis=1).repeat(nt_i)
        # write changes to the file
        s0 = s[cut]
        s0[invert(good)] = 0
        s[cut] = s0  # set all sizes of the removed clusters to 0
        group['Clusters']['Size'][...] = s
        for name in ['U', 'V', 'X', 'Y', 'Charge']:
            data = array(group['Clusters'][name])[good]
            del group['Clusters'][name]
            group['Clusters'].create_dataset(name, data=data)

    @update_pbar
    def add_tracks(self, group, tree):
        gr = group.create_group('Tracks')
        data = get_tree_vec(tree, [f'trk_{n}' for n in ['u', 'v', 'col', 'row', 'du', 'dv']], dtype='f2')
        for i, name in enumerate(['U', 'V', 'X', 'Y', 'dU', 'dV']):
            gr.create_dataset(name, data=data[i])

    @update_pbar
    def add_clusters(self, group, tree):
        gr = group.create_group('Clusters')

        # tracking
        t = self.TrackFile.Get(group.name.strip('/')).Get('Clusters')
        t.SetEstimate(-1)
        gr.create_dataset('N', data=get_tree_vec(t, 'NClusters', dtype='u1'))

        # matching
        branches = [f'clu_{n}' for n in ['size', 'col', 'row', 'u', 'v']]
        data = get_tree_vec(tree, branches, dtype=['u2'] + ['f2'] * len(branches))
        cluster_size = data[0]
        gr.create_dataset('Size', data=cluster_size)
        for i, name in enumerate(['X', 'Y', 'U', 'V'], 1):
            gr.create_dataset(name, data=data[i][cluster_size > 0])  # filter out the nan events

    def add_trigger_info(self, group):
        f = TFile(self.ROOTFileName)
        tree = f.Get(group.name.strip('/')).Get('Hits')
        tree.SetEstimate(tree.GetEntries())
        data = get_tree_vec(tree, ['TriggerPhase', 'TriggerCount'], dtype='u1')
        for i, name in enumerate(['TriggerPhase', 'TriggerCount']):
            group.create_dataset(name, data=data[i])

    def add_cluster_charge(self, match_tree, group, dut_nr):
        plane = Plane(dut_nr + self.NTelPlanes, self.Config('DUT'))
        hit_list = self.get_hits(match_tree, group, dut_nr)
        info(f'calculating cluster charges for DUT Plane {dut_nr} ... ')
        clusters = []
        self.PBar.start(group['Clusters']['X'].size)
        for hits in split(hit_list, cumsum(group['Clusters']['Size'])[:-1].astype('i')):
            clusters.append(self.clusterise(hits))
            self.PBar.update()
        t = info(f'recalculating cluster positions according to new charge weights for DUT {dut_nr} ... ', endl=False)
        clusters = concatenate(clusters)
        group = group['Clusters']
        group.create_dataset('Charge', data=array([nan if cluster is None else cluster.get_charge() for cluster in clusters], dtype='f2'))
        x = array([cluster.get_x() for cluster in clusters])
        y = array([cluster.get_y() for cluster in clusters])
        group['U'][...] = array(group['U']) + (x - array(group['X'])) * plane.PX  # correct for the new charge weights
        group['V'][...] = array(group['V']) + (y - array(group['Y'])) * plane.PY
        group['X'][...] = x
        group['Y'][...] = y
        add_to_info(t)

    def get_hits(self, match_tree, group, dut_nr):
        t = info('getting hit charges for DUT Plane {} ...'.format(dut_nr), endl=False)
        calibration = self.get_calibration(dut_nr)
        x, y, adc = get_tree_vec(match_tree, ['hit_col', 'hit_row', 'hit_value'], dtype=['u2', 'u2', 'u1'])
        hits = array([x, y, calibration(x, y, adc)]).T
        group.create_dataset('CalChiSquare', data=calibration.get_chi2s())
        add_to_info(t)
        return hits
    # endregion HDF5 Converter
    # ----------------------------------------

    # ----------------------------------------
    # region PROTEUS
    def run_proteus(self, prog, out, geo=None, sub_section=None, f=None):
        if not self.check_root_version():
            return
        chdir(self.ProteusDataDir)
        sub_section = '' if sub_section is None else ' -u {}'.format(sub_section)
        geo = '' if geo is None else ' -g {}'.format(geo)
        cmd = '{} {} {}{}{}'.format(join(self.ProteusSoftDir, 'bin', prog), join(self.SaveDir, self.ROOTFileName) if f is None else f, out, geo, sub_section)
        info(cmd)
        check_call(cmd.split())

    def noise_scan(self):
        """ step 1: find noisy pixels with proteus. """
        ensure_dir('mask')
        self.run_proteus('pt-noisescan', join('mask', 'all'))

    def align(self):
        """ step 2: align the telescope in several steps with proteus. """
        ensure_dir('align')
        for i in range(len(self.AlignSteps)):
            self.align_step(i)

    def align_step(self, step, force=False):
        step_name = self.AlignSteps[step]
        if not file_exists(self.make_toml_name(step_name)) or force:
            self.run_proteus('pt-align', join('align', step_name), self.make_toml_name(self.AlignSteps[step - 1] if step else None), step_name)
        else:
            warning('geo file already exists!')

    def track(self):
        """ step 3: based on the alignment generate the tracks with proteus. """
        self.run_proteus('pt-track', 'track-{:04d}'.format(self.RunNumber), self.make_toml_name())

    def match(self):
        """ step 4: match the tracks to the hits in the DUT with proteus. """
        self.run_proteus('pt-match', 'match-{:04d}'.format(self.RunNumber), self.make_toml_name(), f='track-{:04d}-data.root'.format(self.RunNumber))
    # endregion PROTEUS
    # ----------------------------------------

    def check_alignment(self, step=0, p1=0, p2=1):
        f = TFile(join(self.ProteusDataDir, 'align', '{}-hists.root'.format(self.AlignSteps[step])))
        tdir = f.Get('validation').Get('correlations').Get('M{}-M{}'.format(p1, p2))
        c = self.Draw.canvas('c_cor{}'.format(step), 'Correlations Plane {} & {}'.format(p1, p2), divide=(2, 2))
        for i, h in enumerate([tdir.Get(name) for name in ['correlation_x', 'correlation_y', 'difference_x', 'difference_y']], 1):
            self.Draw(h, canvas=c.cd(i), draw_opt='colz' if i in [1, 2] else '')
        self.Draw.Objects.append([f])

    def make_toml_name(self, name=None, d='align', typ='geo'):
        return join(d, '{}-{}.toml'.format(self.AlignSteps[-1] if name is None else name, typ))

    def draw_track(self, ev=10):
        from ROOT import TFile
        from draw import Draw, format_histo
        d = Draw()
        f = TFile(join(self.SaveDir, self.ROOTFileName))
        trees = [f.Get(key.GetName()).Get('Hits') for key in f.GetListOfKeys() if 'Plane' in key.GetName()]
        for tree in trees:
            tree.GetEntry(ev)
        x, y, n = [array([array(getattr(tree, att)) for tree in trees]) for att in ['PixX', 'PixY', 'NHits']]
        zed = concatenate([[i_z] * i_n for i_z, i_n in zip(self.get_z_positions(), n)])
        sx = concatenate([([.0184] if i < 6 else [.15]) * i_n for i, i_n in enumerate(n)])
        sy = concatenate([([.0184] if i < 6 else [.1]) * i_n for i, i_n in enumerate(n)])
        gx = d.make_tgrapherrors('gx', 'x', x=zed, y=concatenate(x) * sx)
        gy = d.make_tgrapherrors('gy', 'y', x=zed, y=concatenate(y) * sy)
        for g in [gx, gy]:
            format_histo(g, x_tit='z [mm]')
            d.histo(g, draw_opt='ap')
        self.Draw.add(trees, f, d)

    def remove_root_files(self):
        for f in self.FileNames:
            if f.endswith('.root') or f.endswith('.hdf5'):
                remove_file(f)


if __name__ == '__main__':
    from analysis import Analysis

    parser = ArgumentParser()
    parser.add_argument('run', nargs='?', default=11)
    # noinspection PyTypeChecker
    parser.add_argument('dut', nargs='?', default=1, type=int)
    pargs = parser.parse_args()
    a = Analysis()
    run = DESYRun(pargs.run, pargs.dut, a.BeamTest.Path, a.Config, single_mode=True)
    z = DESYConverter(run.TCDir, run.Number, a.Config)
    # z.run()
