#!/usr/bin/env python
# --------------------------------------------------------
#       convert eudaq-2 raw files to hdf5
# created on August 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from __future__ import print_function

from draw import *
from ROOT import TBrowser
from os import chdir, remove
from os.path import basename, join
from argparse import ArgumentParser
from converter import Converter
from subprocess import check_call
from glob import glob
from numpy import concatenate


class DESYConverter(Converter):

    def __init__(self, filepath, config):
        Converter.__init__(self, filepath, config)

        self.EUDAQDir = join(self.SoftDir, self.Config.get('SOFTWARE', 'eudaq2'))
        self.ProteusSoftDir = join(self.SoftDir, self.Config.get('SOFTWARE', 'proteus'))
        self.ProteusDataDir = join(self.TCDir, 'proteus')
        self.SaveDir = join(self.TCDir, 'hdf5')
        self.ROOTFileName = '{}.root'.format(basename(self.RawFilePath).split('_')[0])

        self.AlignSteps = self.get_align_steps()
        self.Draw = Draw()

        self.B = []

    def get_z_positions(self):
        with open(join(self.ProteusDataDir, 'geometry.toml')) as f:
            return array([float(line.split(',')[-1].strip(']\n\r')) for line in f.readlines() if line.startswith('offset')])

    def get_align_steps(self):
        with open(join(self.ProteusDataDir, 'analysis.toml')) as f:
            return [line.split(']')[0].replace('[align.', '') for line in f.readlines() if line.startswith('[align.')]

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
        self.B.append([trees, f, d])

    def convert_raw_to_root(self):
        ensure_dir(self.SaveDir)
        chdir(self.SaveDir)
        cmd = '{} -i {} -o {}'.format(join(self.EUDAQDir, 'bin', 'euCliConverter'), self.RawFilePath, self.ROOTFileName)
        info('Convert {} to {}.root using eudaq-2\n'.format(basename(self.RawFilePath), self.ROOTFileName))
        check_call(cmd.split())
        for filename in glob('AutoDict_vector*'):
            remove(filename)

    def run_proteus(self, prog, out, geo=None, sub_section=None, f=None):
        chdir(self.ProteusDataDir)
        sub_section = '' if sub_section is None else ' -u {}'.format(sub_section)
        geo = '' if geo is None else ' -g {}'.format(geo)
        cmd = '{} {} {}{}{}'.format(join(self.ProteusSoftDir, 'bin', prog), join(self.SaveDir, self.ROOTFileName) if f is None else f, out, geo, sub_section)
        info(cmd)
        check_call(cmd.split())

    def noise_scan(self):
        ensure_dir('mask')
        self.run_proteus('pt-noisescan', join('mask', 'all'))

    def align(self):
        ensure_dir('align')
        for i in xrange(len(self.AlignSteps)):
            self.align_step(i)

    def align_step(self, step, force=False):
        step_name = self.AlignSteps[step]
        if not file_exists(self.make_toml_name(step_name)) or force:
            self.run_proteus('pt-align', join('align', step_name), self.make_toml_name(self.AlignSteps[step - 1] if step else None), step_name)
        else:
            warning('geo file already exists!')

    def track(self):
        self.run_proteus('pt-track', 'track', self.make_toml_name())

    def match(self):
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


if __name__ == '__main__':
    from analysis import Analysis
    from desy_run import DESYRun

    p = ArgumentParser()
    p.add_argument('run', nargs='?', default=11)
    # noinspection PyTypeChecker
    p.add_argument('dut', nargs='?', default=1, type=int)
    args = p.parse_args()
    a = Analysis()
    r = DESYRun(args.run, args.dut, a.TCDir, a.Config, single_mode=True)
    z = DESYConverter(r.RawFileName, a.Config)
    # z.run()
