#!/usr/bin/env python
# --------------------------------------------------------
#       run module that handles the files and information of a single run
# created on October 5th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from utility.utils import print_table
from plotting.utils import load_json
from src.analysis import Analysis, Path, choose
from src.dut import DUT


class Run:
    """ Run class containing all the information for a single run from the tree and the json file. """

    def __init__(self, run_number, dut_number, tc_dir: Path, single_mode=False):

        # Main Info
        self.Number = int(run_number)
        self.TCDir = tc_dir
        self.SingleMode = single_mode

        # Files
        self.FileName = self.TCDir.joinpath('data', f'run{self.Number:04d}.hdf5')

        # Info
        self.Logs = self.load_logs()
        self.DUT = DUT(dut_number, self.Logs)

        # Times  TODO: Fix later with real timestamps from the data
        self.StartTime = self.Logs['start']
        self.EndTime = self.Logs['end']

    def __str__(self):
        return str(self.Number)

    def __repr__(self):
        return f'Run {self}'

    def __format__(self, format_spec):
        return f'{self.Number:{format_spec}}'

    def load_logs(self) -> dict:
        data = load_json(self.TCDir.joinpath(Analysis.Config.get('data', 'runlog file')))
        if self.SingleMode:
            return data[str(self.Number)]
        return {}

    def print_info(self):
        print_table(rows=[[key, str(value)] for key, value in self.Logs.items()])

    @classmethod
    def from_ana(cls, run_number, dut=0, ana: Analysis = None, single_mode=False):
        ana = choose(ana, Analysis)
        return cls(run_number, dut, ana.BeamTest.Path, single_mode)


if __name__ == '__main__':
    z = Run.from_ana(41, 1, single_mode=True)
