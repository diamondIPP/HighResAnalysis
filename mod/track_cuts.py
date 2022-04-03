#!/usr/bin/env python
# --------------------------------------------------------
#       cuts for the tracks
# created on March 30th 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from mod.dut_cuts import DUTCut
from numpy import ndarray


class TrackCut(DUTCut):
    def __init__(self, ana):
        self.DUTAna = ana.Ana
        super().__init__(ana, meta_sub_dir='track_cuts')

    def __call__(self, cut=None, data=None):
        if data is None:
            return super().__call__(cut)
        if data.size == self.Ana.NEvents:
            data = self.ev2trk(data)
        if type(cut) is not ndarray and (cut is None or cut == -1):
            return data[self.DUTAna.Cut['cluster']][self.DUTAna.Cut(cut)]
        cut = super().__call__(cut)
        return data if cut is ... else data[cut] if cut.size == self.Ana.N else data[self.DUTAna.Cut['cluster']][cut]

    def make(self, redo=False):
        self.register('tp', self.ev2trk(self.make_trigger_phase(_redo=redo)), 10, 'trigger phase')
        self.register('res', self.DUTAna.REF.Cut.make_residual(_redo=redo), 20, 'tracks with a small residual in the REF')
        self.register('fid', self.make_fiducial(_redo=redo), 30, 'tracks in fiducial area')
        self.register('tstart', self.make_start_time(_redo=redo), 40, 'exclude first events')
        self.register('chi2', self.make_chi2(_redo=redo), 50, 'small chi2')

    def get_xy(self):
        # todo: check and use get_coods()!
        return self.Ana.get_xy(cut=False, local=True)

    def make_trk(self, trks):
        return self.make_ev(trks, self.Ana.N)
