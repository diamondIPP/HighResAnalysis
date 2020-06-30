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
from numpy import concatenate
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
        self.EUDAQDir = join(self.SoftDir, self.Config.get('SOFTWARE', 'eudaq2'))
        self.ProteusSoftDir = join(self.SoftDir, self.Config.get('SOFTWARE', 'proteus'))
        self.ProteusDataDir = join(self.DataDir, 'proteus')
        self.SaveDir = join(self.DataDir, 'data')

        # FILENAMES
        self.ROOTFileName = join(self.SaveDir, 'run{:06d}.root'.format(self.RunNumber))
        self.NoiseName = join('mask', 'all')
        self.AlignSteps = self.get_align_steps()
        self.ProteusROOTFiles = [join(self.ProteusDataDir, '{}.root'.format(n)) for n in ['track-data', 'match-trees']]
        self.FinalFileName = join(self.SaveDir, 'run{:04d}.hdf5'.format(self.RunNumber))
        self.FileNames = [self.ROOTFileName] + [self.get_toml_path(p) for p in [self.NoiseName, join('align', self.AlignSteps[-1])]] + self.ProteusROOTFiles + [self.FinalFileName]

        self.NSteps = len(self.FileNames)
        self.AtStep = self.find_step()

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

    def find_step(self):
        for i in reversed(range(self.NSteps)):
            if file_exists(self.FileNames[i]):
                return i + 1
        return 0

    def get_steps(self):
        steps = [self.convert_raw_to_root, self.noise_scan, self.align, self.track, self.match, self.convert_root_to_hdf5]
        return steps[self.AtStep:]

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

    def get_toml_path(self, path):
        name = basename(path)
        d = None if dirname(path) == name else dirname(path)
        return join(self.ProteusDataDir, self.make_toml_name(name, d, typ='geo' if d is None else d))
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

    def get_tree(self, plane=0):
        root_file = TFile(self.FileNames[-3])
        return root_file, root_file.Get('Plane{}'.format(plane))

    def convert_root_to_hdf5(self):
        """ step 5: convert root file to hdf5 file. """
        start_time = info('Start root->hdf5 conversion ...')
        root_file = TFile(self.FileNames[-3])
        if file_exists(self.FinalFileName):
            remove_file(self.FinalFileName)
        f = h5py.File(self.FinalFileName, 'w')

        # tracks
        t = root_file.Get('Tracks')
        names = [n.GetName() for n in t.GetListOfBranches()]  # NTracks, Chi2, Dof, X, Y, SlopeX, SlopeY, Cov
        f.create_dataset('NTracks', data=get_root_vec(t, var=names[0], dtype='u1'))
        n = t.Draw(':'.join(names[1:-1]), '', 'goff')
        types = ['f2', 'u1', 'f2', 'f2', 'f2', 'f2']
        for i, (name, typ) in enumerate(zip(names[1:-1], types)):
            f.create_dataset(name, data=get_root_vec(t, n, i, dtype=typ))

        # planes
        names = [n.GetName() for n in root_file.GetListOfKeys() if n.GetName().startswith('Plane')]
        info('reading branches for {} planes'.format(len(names)))
        hit_types = {'X': 'u2', 'Y': 'u2', 'ADC': 'u1', 'HitInCluster': 'u1'}
        plane_types = {'X': 'f4', 'Y': 'f4', 'VarX': 'f2', 'VarY': 'f2', 'NTracks': 'u1'}
        intercept_types = {'X': 'f4', 'Y': 'f4', 'SlopeX': 'f2', 'SlopeY': 'f2', 'NTracks': 'u1'}
        self.PBar.start(sum(len(lst) for lst in [hit_types, plane_types, intercept_types]) * len(names))
        mask = self.load_mask()
        for i, name in enumerate(names):
            g = f.create_group(name)

            # hits: NHits, PixX, PixY, Timing, Value, HitInCluster
            self.convert_plane_data(root_file, g, name, branch_name='Hits', types=hit_types)

            # clusters:  NClusters, Col, Row, VarCol, VarRow, CovColRow, Timing, Value, Track
            self.convert_plane_data(root_file, g, name, 'Clusters', plane_types, exclude='Value')  # adc doesn't make sense for clusters

            # intercepts
            self.convert_plane_data(root_file, g, name, 'Intercepts', intercept_types)  # adc doesn't make sense for clusters

            # mask
            if i in mask:
                g.create_dataset('Mask', data=mask[i].astype('u2'))

        add_to_info(start_time, 'Finished conversion in')

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
        self.run_proteus('pt-track', 'track', self.make_toml_name())

    def match(self):
        """ step 4: match the tracks to the hits in the DUT with proteus. """
        self.run_proteus('pt-match', 'match', self.make_toml_name(), f='track-data.root')

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


if __name__ == '__main__':
    from analysis import Analysis
    from desy_run import DESYRun

    parser = ArgumentParser()
    parser.add_argument('run', nargs='?', default=11)
    # noinspection PyTypeChecker
    parser.add_argument('dut', nargs='?', default=1, type=int)
    pargs = parser.parse_args()
    a = Analysis()
    r = DESYRun(pargs.run, pargs.dut, a.TCDir, a.Config, single_mode=True)
    z = DESYConverter(r.TCDir, r.Number, a.Config)
    # z.run()
