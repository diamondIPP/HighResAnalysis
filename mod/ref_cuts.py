#!/usr/bin/env python
# --------------------------------------------------------
#       cuts for the reference detector
# created on March 30th 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from mod.dut_cuts import DUTCut


class RefCut(DUTCut):
    def __init__(self, ana):
        self.Ana = ana
        super().__init__(ana, meta_sub_dir='ref_cuts')
        self.N = self.Ana.N
