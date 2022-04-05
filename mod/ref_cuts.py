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

    def make(self, redo=False):
        self.register('cluster', self.make_cluster(), 90, 'tracks with a cluster')

    def make_trk_residual(self, redo=False):
        return self.pl2trk(self.make_residual(_redo=redo))
