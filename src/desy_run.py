#!/usr/bin/env python
# --------------------------------------------------------
#       run module that handles the files and information of a single run taken at DESY
# created on February 4th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from src.run import Run
from plotting.utils import load_json


class DESYRun(Run):
    """ Run class containing all the information for a single run from the tree and the json file. """

    def __init__(self, run_number, dut_nr, tc_dir, config, single_mode=True):

        self.TreeName = 'tree'
        Run.__init__(self, run_number, dut_nr, tc_dir, config, single_mode)

    def load_run_logs(self):
        data = load_json(self.TCDir.joinpath(self.Config.get('data', 'run log file')))
        if self.SingleMode:
            return data[str(self.Number)]

    def load_file_name(self):
        return self.TCDir.joinpath('data', f'run{self.Number:04d}.hdf5')


if __name__ == '__main__':

    from analysis import Analysis
    a = Analysis()
    z = DESYRun(4, 1, a.BeamTest.Path, a.Config, single_mode=True)
