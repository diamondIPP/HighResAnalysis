#!/usr/bin/env python
# --------------------------------------------------------
#       run module that handles the files and information of a single run
# created on October 5th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from utility.utils import *
from plotting.utils import load_json, Config
from os.path import join
from src.analysis import Analysis
from src.dut import DUT


class Run:
    """ Run class containing all the information for a single run from the tree and the json file. """

    def __init__(self, run_number, dut_number, tc_dir: Path, config: Config, single_mode=False):

        # Main Info
        self.Number = int(run_number)
        self.Config = config
        self.TCDir = tc_dir
        self.SingleMode = single_mode

        # Files
        self.FileName = self.TCDir.joinpath('data', f'run{self.Number:04d}.hdf5')

        # Info
        self.Logs = self.load_logs()
        self.DUT = DUT(dut_number, self.Logs, self.Config)

        # Times  TODO: Fix later with real timestamps from the data
        self.StartTime = self.Logs['start']
        self.EndTime = self.Logs['end']

    def __str__(self):
        return str(self.Number)

    def __repr__(self):
        return f'Run {self}'

    def load_logs(self) -> dict:
        data = load_json(self.TCDir.joinpath(self.Config.get('data', 'runlog file')))
        if self.SingleMode:
            return data[str(self.Number)]
        return {}

    def print_info(self):
        print_table(rows=[[key, str(value)] for key, value in self.Logs.items()])


if __name__ == '__main__':
    a = Analysis()
    z = Run(41, 0, a.BeamTest.Path, a.Config, single_mode=True)
