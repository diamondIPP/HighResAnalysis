#!/usr/bin/env python
# --------------------------------------------------------
#       cuts for the telescope planes
# created on March 30th 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from numpy import all
from mod.dut_cuts import DUTCut


class TelCut(DUTCut):
    def __init__(self, ana):
        super().__init__(ana, meta_sub_dir='tel_cuts')

    def make(self, redo=False):
        pass

    def make_all_cluster(self):
        return all([self.make_cluster(pl) for pl in self.Ana.planes], axis=0)
