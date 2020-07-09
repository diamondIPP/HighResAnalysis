#!/usr/bin/env python
# --------------------------------------------------------
#       convert EUDAQ-2 raw files to hdf5
# created on August 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from __future__ import print_function

from draw import *
from os import chdir, remove
from os.path import basename, join
from argparse import ArgumentParser
from converter import Converter
from subprocess import check_call
from glob import glob
from numpy import concatenate, cumsum, split, sum, in1d, where, argmax
import toml
from calibration import Calibration
from desy_run import DESYRun


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
        self.EUDAQDir = join(self.SoftDir, self.Config.get('SOFTWARE', 'eudaq2'))
        self.ProteusSoftDir = join(self.SoftDir, self.Config.get('SOFTWARE', 'proteus'))
        self.ProteusDataDir = join(self.DataDir, 'proteus')

        # FILENAMES
        self.ROOTFileName = join(self.SaveDir, 'run{:06d}.root'.format(self.RunNumber))
        self.NoiseName = join('mask', 'all')
        self.AlignSteps = self.get_align_steps()
        self.ProteusROOTFiles = [join(self.ProteusDataDir, '{}-{:04d}-{}.root'.format(n, self.RunNumber, m)) for n, m in [('track', 'data'), ('match', 'trees')]]
        self.FinalFileName = join(self.SaveDir, 'run{:04d}.hdf5'.format(self.RunNumber))
        self.FileNames = [self.ROOTFileName] + self.get_align_files() + self.ProteusROOTFiles + [self.FinalFileName]

        self.NSteps = len(self.FileNames)
        self.NTelPlanes = self.Config.getint('TELESCOPE', 'planes')
        self.NDUTPlanes = self.Config.getint('DUT', 'planes')

        self.Draw = Draw()

    # TODO: add check for alignment or align every run?

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
            warning('you require ROOT6 to run proteus! Current version: {}'.format(version))
        return version.startswith('6')

    def load_mask(self):
        d = toml.load(self.FileNames[1])['sensors']
        return {dic['id']: array(dic['masked_pixels']) for dic in d}
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_z_positions(self):
        with open(join(self.ProteusDataDir, 'geometry.toml')) as f:
            return array([float(line.split(',')[-1].strip(']\n\r')) for line in f.readlines() if line.startswith('offset')])

    def get_align_steps(self):
        with open(join(self.ProteusDataDir, 'analysis.toml')) as f:
            return [line.split(']')[0].replace('[align.', '') for line in f.readlines() if line.startswith('[align.')]

    def get_align_files(self):
        return [join(self.ProteusDataDir, name) for name in [self.make_toml_name('all', 'mask', 'mask'), self.make_toml_name()]]

    def get_calibration(self, dut_number):
        return Calibration(DESYRun(self.RunNumber, dut_number, self.DataDir, self.Config))
    # endregion GET
    # ----------------------------------------

    def convert_raw_to_root(self):
        """ step 0: convert raw file to root file with eudaq"""
        if self.RawFilePath is None:
            critical('raw file does not for run: {}'.format(self.RunNumber))
        ensure_dir(self.SaveDir)
        chdir(self.SaveDir)
        cmd = '{} -i {} -o {}'.format(join(self.EUDAQDir, 'bin', 'euCliConverter'), self.RawFilePath, self.ROOTFileName)
        info('Convert {} to {}.root using eudaq-2\n'.format(basename(self.RawFilePath), self.ROOTFileName))
        check_call(cmd.split())
        for filename in glob('AutoDict_vector*'):
            remove(filename)

    def load_root_file(self):
        return TFile(self.FileNames[-3])

    def convert_root_to_hdf5(self):
        """ step 5: convert root file to hdf5 file. """
        if file_exists(self.FinalFileName):
            remove_file(self.FinalFileName)
        start_time = info('Start root->hdf5 conversion ...')
        # files
        f = h5py.File(self.FinalFileName, 'w')
        root_file = self.load_root_file()

        self.add_tracks(root_file, f)
        t_vec = get_root_vec(root_file.Get('Event'), var='TriggerTime', dtype='u8')
        f.create_dataset('Time', data=((t_vec - t_vec[0]) / 1e9).astype('f2'))  # convert time vec to seconds
        self.add_planes(root_file, f)

        for i in range(self.NDUTPlanes):
            group = f['Plane{}'.format(i + self.NTelPlanes)]
            self.add_hit_charge(group, i)
            self.add_trigger_info(group)
            self.add_cluster_charge(group)

        add_to_info(start_time, '\nFinished conversion in')

    @staticmethod
    def add_tracks(root_file, hdf5_file):
        t0 = info('adding track information ... ', overlay=True)
        tree = root_file.Get('Tracks')
        names = [n.GetName() for n in tree.GetListOfBranches()]  # NTracks, Chi2, Dof, X, Y, SlopeX, SlopeY, Cov
        hdf5_file.create_dataset('NTracks', data=get_root_vec(tree, var=names[0], dtype='u1'))
        n = tree.Draw(':'.join(names[1:-1]), '', 'goff')
        types = ['f2', 'u1', 'f2', 'f2', 'f2', 'f2']
        for i, (name, typ) in enumerate(zip(names[1:-1], types)):
            hdf5_file.create_dataset(name, data=get_root_vec(tree, n, i, dtype=typ))
        add_to_info(t0)

    def add_planes(self, root_file, hdf5_file):
        names = [key.GetName() for key in root_file.GetListOfKeys() if key.GetName().startswith('Plane')]
        info('reading branches for {} planes ... '.format(len(names)))
        hit_types = {'X': 'u2', 'Y': 'u2', 'ADC': 'u1', 'HitInCluster': 'u1'}
        plane_types = {'X': 'f2', 'Y': 'f2', 'VarX': 'f2', 'VarY': 'f2', 'NTracks': 'u1'}
        intercept_types = {'X': 'f2', 'Y': 'f2', 'SlopeX': 'f2', 'SlopeY': 'f2', 'NTracks': 'u1'}
        self.PBar.start(sum(len(lst) for lst in [hit_types, plane_types, intercept_types]) * len(names))
        mask = self.load_mask()
        for i, name in enumerate(names):
            g = hdf5_file.create_group(name)
            # mask
            if i in mask:
                g.create_dataset('Mask', data=mask[i].astype('u2'))
            # hits: NHits, PixX, PixY, Timing, Value, HitInCluster
            self.convert_plane_data(root_file, g, name, branch_name='Hits', types=hit_types)
            # clusters:  NClusters, Col, Row, VarCol, VarRow, CovColRow, Timing, Value, Track
            self.convert_plane_data(root_file, g, name, 'Clusters', plane_types, exclude='Value')  # adc doesn't make sense for clusters
            # intercepts
            self.convert_plane_data(root_file, g, name, 'Intercepts', intercept_types)  # adc doesn't make sense for clusters

    def add_trigger_info(self, group):
        f = TFile(self.ROOTFileName)
        tree = f.Get(group.name.strip('/')).Get('Hits')
        trigger_phase = get_root_vec(tree, var='TriggerPhase', dtype='u1')
        trigger_count = get_root_vec(tree, var='TriggerCount', dtype='u1')
        group.create_dataset('TriggerPhase', data=trigger_phase)
        group.create_dataset('TriggerCount', data=trigger_count)

    def add_hit_charge(self, group, dut_nr):
        t = info('adding hit charges for {} ... '.format(group.name.strip('/')), overlay=True)
        x, y, adc = [array(group['Hits'][name]) for name in ['X', 'Y', 'ADC']]
        calibration = self.get_calibration(dut_nr)
        group['Hits'].create_dataset('Charge', data=calibration(x, y, adc), dtype='f2')
        group.create_dataset('CalChiSquare', data=calibration.get_chi2s())
        add_to_info(t)

    def add_cluster_charge(self, group):
        t = info('adding cluster charges for {} ... '.format(group.name.strip('/')), overlay=True)
        clusters = concatenate([self.clusterise(hit_list) for hit_list in self.get_hits(group)])
        group['Clusters'].create_dataset('Charge', data=array([cluster.get_charge() for cluster in clusters], dtype='f2'))
        group['Clusters']['X'][...] = array([cluster.get_x() for cluster in clusters], dtype='f2')
        group['Clusters']['Y'][...] = array([cluster.get_y() for cluster in clusters], dtype='f2')
        group['Clusters'].create_dataset('Size', data=array([cluster.get_size() for cluster in clusters], dtype='u1'))
        add_to_info(t)

    @staticmethod
    def get_hits(group):
        data = array([array(group['Hits'][s]) for s in ['X', 'Y', 'Charge']]).T
        hs = cumsum(group['Hits']['NHits'])  # indices to sort hit array for events
        if 'Mask' not in group.keys():
            return split(data, hs[:-1])
        mask = group['Mask'] if 'Mask' in group.keys() else array([])
        keys = data[:, 0].astype('i') * 10000 + data[:, 1].astype('i')
        mask_keys = mask[:, 0].astype('i') * 10000 + mask[:, 1]
        i_mask = where(in1d(keys, mask_keys))[0]
        data = data[in1d(keys, mask_keys, invert=True)]
        i_sub = zeros(hs.size, dtype='i')
        for i in i_mask:
            i_sub[argmax(hs == i + 1):] += 1
        return split(data, (hs - i_sub)[:-1].astype('i'))

    def convert_plane_data(self, root_file, group, plane_name, branch_name, types, exclude=None):
        g0 = group.create_group(branch_name)
        t0 = root_file.Get(plane_name).Get(branch_name)  # NHits, PixX, PixY, Timing, Value, HitInCluster
        n_name = 'N{}'.format(branch_name)
        g0.create_dataset(n_name, data=get_root_vec(t0, var=n_name, dtype='u2'))
        exclude = make_list(exclude, dtype=list) + ['Timing', 'Cov', 'CovColRow']  # Timing & CovColRow is empty, Cov too much data
        names = [b.GetName() for b in t0.GetListOfBranches() if b.GetName() not in exclude][1:]
        n = t0.Draw(':'.join(names), '', 'goff')
        for j, (branch_name, typ) in enumerate(types.items()):
            data = get_root_vec(t0, n, j, typ)
            self.PBar.update()
            if branch_name == 'ADC' and mean(data) == 1:  # don't save empty data... not all planes have pulse height information
                continue
            g0.create_dataset(branch_name, data=data)

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
        self.run_proteus('pt-noisescan', self.NoiseName)

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

    def check_alignment(self, step=0, p1=0, p2=1):
        f = TFile(join(self.ProteusDataDir, 'align', '{}-hists.root'.format(self.AlignSteps[step])))
        tdir = f.Get('validation').Get('correlations').Get('M{}-M{}'.format(p1, p2))
        c = self.Draw.make_canvas('c_cor{}'.format(step), 'Correlations Plane {} & {}'.format(p1, p2), divide=(2, 2))
        self.Draw.format_statbox(all_stat=True)
        for i, h in enumerate([tdir.Get(name) for name in ['correlation_x', 'correlation_y', 'difference_x', 'difference_y']], 1):
            self.Draw.draw_histo(h, canvas=c.cd(i), draw_opt='colz' if i in [1, 2] else '')
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
            d.draw_histo(g, draw_opt='ap')
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
    r = DESYRun(pargs.run, pargs.dut, a.TCDir, a.Config, single_mode=True)
    z = DESYConverter(r.TCDir, r.Number, a.Config)
    # z.run()
