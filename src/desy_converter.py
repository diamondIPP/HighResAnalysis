#!/usr/bin/env python
# --------------------------------------------------------
#       convert eudaq-2 raw files to hdf5
# created on August 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from __future__ import print_function

from utils import *
from os import chdir, remove
from os.path import basename, join, dirname
from argparse import ArgumentParser
from collections import OrderedDict
from progressbar import Bar, ETA, FileTransferSpeed, Percentage, ProgressBar
from pickle import load, dump
from draw import Draw, ufloat
from converter import Converter
from subprocess import check_call
from glob import glob


class DESYConverter(Converter):

    def __init__(self, filepath, config):
        Converter.__init__(self, filepath, config)

        self.EUDAQDir = join(self.SoftDir, self.Config.get('SOFTWARE', 'eudaq2'))
        self.ProteusDir = join(self.SoftDir, self.Config.get('SOFTWARE', 'proteus'))
        self.SaveDir = join(self.TCDir, 'hdf5')
        self.ROOTFileName = '{}.root'.format(basename(self.RawFilePath).split('_')[0])

    def convert_raw_to_root(self):
        ensure_dir(self.SaveDir)
        chdir(self.SaveDir)
        cmd = '{} -i {} -o {}'.format(join(self.EUDAQDir, 'bin', 'euCliConverter'), self.RawFilePath, self.ROOTFileName)
        info('Convert {} to {}.root using eudaq-2\n'.format(basename(self.RawFilePath), self.ROOTFileName))
        check_call(cmd.split())
        for filename in glob('AutoDict_vector*'):
            remove(filename)

    def run_proteus(self, prog):
        chdir(join(self.TCDir, 'proteus'))
        cmd = '{} {} run{:04d}'.format(join(self.ProteusDir, 'bin', prog), join(self.SaveDir, self.ROOTFileName), int(remove_letters(self.ROOTFileName)))
        info(cmd)
        check_call(cmd.split())

    def noise_scan(self):
        self.run_proteus('pt-noisescan')


if __name__ == '__main__':
    from analysis import Analysis
    from desy_run import DESYRun

    p = ArgumentParser()
    p.add_argument('run', nargs='?', default='')
    # noinspection PyTypeChecker
    p.add_argument('dut', nargs='?', default=1, type=int)
    args = p.parse_args()
    a = Analysis()
    r = DESYRun(args.run, args.dut, a.TCDir, a.Config, single_mode=True)
    z = DESYConverter(r.RawFileName, a.Config)
    # z.run()
