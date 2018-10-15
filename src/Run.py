#!/usr/bin/env python
# --------------------------------------------------------
#       run module that handles the files and information of a single run
# created on October 5th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from utils import *
from os.path import join, basename
from Analysis import Analysis
from subprocess import call
from Converter import Converter
from os import environ, remove, devnull
from collections import OrderedDict


class Run:
    """ Run class containing all the information for a single run from the tree and the json file. """

    def __init__(self, run_number, plane, tc_dir, config, single_mode=False):

        self.RunNumber = run_number
        self.Plane = plane
        self.Config = config
        self.TCDir = tc_dir
        self.SingleMode = single_mode

        self.RunInfo = self.load_run_info()
        self.RunLogs = self.load_run_logs()
        self.DutName = self.load_dut_name()

        # Files and Trees
        self.RawFileName = self.load_raw_file_name()
        self.FileName = join(self.TCDir, self.DutName, 'Clustered_{}.root'.format(str(self.RunNumber).zfill(3)))
        self.File = self.load_file()
        self.Tree = self.File.Get('Hits')

    def load_run_info(self):
        if self.SingleMode:
            info('Using single raw files')
            return None
        data = load_json(join(self.TCDir, self.Config.get('MAIN', 'run plan file')))
        data = data.values()[self.Plane][str(self.RunNumber)]
        data['Batches'] = data['Batches'].split(',')
        return data

    def load_run_logs(self):
        data = load_json(join(self.TCDir, self.Config.get('MAIN', 'run log file')))
        data = {key: value for key, value in data.iteritems() if value['Run Number']['CMS'] not in ['-', '/']}
        if self.RunInfo is None:
            return next(value for value in data.itervalues() if int(value['Run Number']['CMS']) == self.RunNumber)
        return OrderedDict(sorted([(int(value['Run Number']['CMS']), value) for value in data.itervalues() if value['Batch'] in self.RunInfo['Batches']]))

    def load_dut_name(self):
        return load_json(join(self.TCDir, self.Config.get('MAIN', 'run plan file')))['Plane{}'.format(self.Plane)]['Name']

    def load_raw_file_name(self):
        if self.SingleMode:
            return join(self.TCDir, 'cms-raw', 'ljutel_{}.root'.format(str(self.RunNumber).zfill(3)))
        return join(self.TCDir, self.DutName, 'run_{}.root'.format(str(self.RunNumber).zfill(2)))

    def load_file(self):
        if not file_exists(self.FileName):
            self.convert_file()
        return TFile(self.FileName)

    def convert_file(self):
        if not file_exists(self.RawFileName) and not self.SingleMode:
            self.merge_root_files()
        warning('final root file "{}" does not exist. Starting Converter!'.format(self.FileName))
        converter = Converter(self.RawFileName, self.Plane, join(self.TCDir, self.DutName))
        converter.run()
        if not self.SingleMode:
            info('removing raw file "{}"'.format(self.RawFileName))
            remove(self.RawFileName)

    def merge_root_files(self):
        warning('raw file "{}" does not exist. Starting single file merger!'.format(self.RawFileName))
        single_files = [join(self.TCDir, 'cms-raw', 'ljutel_{}.root'.format(n)) for n in self.RunLogs.keys()]
        new_file = join(self.TCDir, self.DutName, 'run_{}.root'.format(str(self.RunNumber).zfill(2)))
        with open(devnull, 'w') as f:
            call([join(environ.get('ROOTSYS'), 'bin', 'hadd'), '-f', new_file] + single_files, stdout=f)
        info('successfully merged the single files to "{}"'.format(basename(new_file)))

    def print_run_info(self):
        for key, value in sorted(load_json(join(self.TCDir, self.Config.get('MAIN', 'run plan file'))).values()[self.Plane].iteritems(), key=lambda k: int(k[0])):
            print key, value


if __name__ == '__main__':
    a = Analysis()
    z = Run(301, 0, a.TCDir, a.Config, single_mode=True)
