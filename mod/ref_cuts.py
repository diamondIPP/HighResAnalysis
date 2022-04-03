#!/usr/bin/env python
# --------------------------------------------------------
#       cuts for the reference detector
# created on March 30th 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from mod.dut_cuts import DUTCut, save_hdf5, choose
from numpy import ndarray, sqrt, zeros


class RefCut(DUTCut):
    def __init__(self, ana):
        self.Ana = ana
        super().__init__(ana, meta_sub_dir='ref_cuts')
        self.N = self.get('cluster').N

    def make(self, redo=False):
        self.register('cluster', self.make_cluster(), 90, 'tracks with a cluster')
        self.register('corr', self.make_corr(), 91, 'tracks with a cluster in DUT and REF')

    def __call__(self, cut=None, pl=0):
        cut = super().__call__(cut, pl)
        return cut if type(cut) is ndarray and cut.size == self.N or cut is ... else self.trk2ref(cut)

    # ----------------------------------------
    # region GENERATE
    @save_hdf5('Corr', arr=True, dtype='?')
    def make_corr(self, _redo=False):
        return self['cluster'] & self.Ana.DUTCut['cluster']

    @save_hdf5('Res', arr=True, dtype='?', suf_args='all')
    def make_residual(self, v=None, _redo=False):
        x, y, (mx, my) = self.Ana.get_du(), self.Ana.get_dv(), self.Ana.get_res_means()
        return self.ref2trk(sqrt((x - mx.n) ** 2 + (y - my.n) ** 2) < choose(v, self.get_config('residuals', dtype=float)))

    def make_corrs(self, cut=True):
        c = self.make_corr() & cut
        return self.trk2ref(c), c[self.Ana.DUTCut['cluster']]
    # endregion GENERATE
    # ----------------------------------------

    def ref2trk(self, cut):
        c = zeros(self.get('cluster').Size, '?')
        c[self['cluster']] = cut
        return c

    def trk2ref(self, cut):
        return ... if type(cut) is bool else cut[self['cluster']]
