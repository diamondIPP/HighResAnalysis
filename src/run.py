#!/usr/bin/env python
# --------------------------------------------------------
#       run module that handles the files and information of a single run
# created on October 5th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from src.utils import *
from os.path import join
from src.analysis import Analysis
from src.dut import DUT


class Run:
    """ Run class containing all the information for a single run from the tree and the json file. """

    def __init__(self, run_number, dut_number, tc_dir, config, single_mode=False):

        # Main Info
        self.Number = int(run_number)
        self.Config = config
        self.TCDir = tc_dir
        self.SingleMode = single_mode

        # Files
        self.FileName = self.load_file_name()

        # Info
        self.Info = self.load_run_info()
        self.Logs = self.load_run_logs()
        self.DUT = DUT(dut_number, self.Logs, self.Config)

        # Times
        # TODO: Fix later with real timestamps from the data
        self.StartTime = self.Logs['start']
        self.EndTime = self.Logs['end']

    def get_dut_nrs(self):
        return [int(remove_letters(key)) for key in self.Logs.keys() if key.startswith('dut')]

    def load_run_info(self):
        pass

    def load_run_logs(self):
        return {}

    def load_dut_name(self):
        pass

    def load_file_name(self):
        pass

    def print_run_info(self):
        for key, value in sorted(load_json(join(self.TCDir, self.Config.get('MAIN', 'run plan file'))).values()[self.DUT.Number].iteritems(), key=lambda k: int(k[0])):
            print(key, value)


if __name__ == '__main__':
    a = Analysis()
    z = Run(301, 1, a.TCDir, a.Config, single_mode=True)
