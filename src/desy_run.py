#!/usr/bin/env python
# --------------------------------------------------------
#       run module that handles the files and information of a single run taken at DESY
# created on February 4th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from utils import *
from run import Run
from os.path import join
from glob import glob


class DESYRun(Run):
    """ Run class containing all the information for a single run from the tree and the json file. """

    def __init__(self, run_number, dut_nr, tc_dir, config, single_mode=True):

        self.TreeName = 'tree'
        Run.__init__(self, run_number, dut_nr, tc_dir, config, single_mode)

    def load_run_logs(self):
        data = load_json(join(self.TCDir, self.Config.get('MAIN', 'run log file')), ordered=True)
        if self.SingleMode:
            return data[str(self.RunNumber)]

    def load_raw_file_name(self):
        names = glob(join(self.TCDir, 'raw', 'run{:06d}*.raw'.format(self.RunNumber)))
        return names[0] if names else None

    def load_file(self):
        return


if __name__ == '__main__':

    from analysis import Analysis
    a = Analysis()
    z = DESYRun(4, 1, a.TCDir, a.Config, single_mode=True)
