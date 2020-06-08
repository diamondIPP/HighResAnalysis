#!/usr/bin/env python
# --------------------------------------------------------
#       run module that handles the files and information of a single run
# created on October 5th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from utils import *
from os.path import join, basename
from analysis import Analysis
from subprocess import call
from converter import Converter
from os import environ, remove, devnull
from dut import DUT


class Run:
    """ Run class containing all the information for a single run from the tree and the json file. """

    def __init__(self, run_number, dut_nr, tc_dir, config, single_mode=False):

        # Main Info
        self.RunNumber = int(run_number)
        self.Config = config
        self.TCDir = tc_dir
        self.SingleMode = single_mode

        # Files
        self.RawFileName = self.load_raw_file_name()
        self.FileName = self.load_file_name()
        self.File = self.load_file()

        # Info
        self.RunInfo = self.load_run_info()
        self.RunLogs = self.load_run_logs()
        self.DUT = DUT(dut_nr, self.RunLogs, config)

        # Times
        # TODO: Fix later with real timestamps from the data
        self.StartTime = self.RunLogs['start']
        self.EndTime = self.RunLogs['end']

    def load_run_info(self):
        pass

    def load_run_logs(self):
        pass

    def load_dut_name(self):
        pass

    def load_raw_file_name(self):
        pass

    def load_file_name(self):
        pass

    def load_file(self):
        if not file_exists(self.FileName):
            self.convert_file()
        return TFile(self.FileName)

    def convert_file(self):
        if not file_exists(self.RawFileName) and not self.SingleMode:
            self.merge_root_files()
        warning('final root file "{}" does not exist. Starting Converter!'.format(self.FileName))
        converter = Converter(self.RawFileName, self.DUT.Number, join(self.TCDir, self.DUT.Name), first_run=int(self.RunLogs.keys()[0]))
        converter.run()
        if not self.SingleMode:
            info('removing raw file "{}"'.format(self.RawFileName))
            remove(self.RawFileName)

    def merge_root_files(self):
        warning('raw file "{}" does not exist. Starting single file merger!'.format(self.RawFileName))
        single_files = [join(self.TCDir, 'cms-raw', 'ljutel_{}.root'.format(n)) for n in self.RunLogs.keys()]
        new_file = join(self.TCDir, self.DUT.Name, 'run_{}.root'.format(str(self.RunNumber).zfill(2)))
        with open(devnull, 'w') as f:
            call([join(environ.get('ROOTSYS'), 'bin', 'hadd'), '-f', new_file] + single_files, stdout=f)
        info('successfully merged the single files to "{}"'.format(basename(new_file)))

    def print_run_info(self):
        for key, value in sorted(load_json(join(self.TCDir, self.Config.get('MAIN', 'run plan file'))).values()[self.DUT.Number].iteritems(), key=lambda k: int(k[0])):
            print key, value


if __name__ == '__main__':
    a = Analysis()
    z = Run(301, 0, a.TCDir, a.Config, single_mode=True)
