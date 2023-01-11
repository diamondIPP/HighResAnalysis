#!/usr/bin/env python
# --------------------------------------------------------
#       cuts for the tracks
# created on March 30th 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from mod.dut_cuts import DUTCut, zeros, save_cut, choose, all, invert, array
from src.dut_analysis import no_trans


class TrackCut(DUTCut):
    def __init__(self, ana):
        super().__init__(ana, meta_sub_dir='track_cuts')

    def __call__(self, cut=None, data=None, pl=None):
        cut = super().__call__(cut)
        if data is None:
            return cut
        if data.size == self.Ana.NEvents:
            return self.ev2trk(data)[cut]
        if cut is not ... and cut.size == self.Ana.N:
            return self.trk2pl(data, pl)[cut]
        return data[cut]

    def make(self, redo=False):
        self.register('tp', self.make_trigger_phase(_redo=redo), 10, 'trigger phase')
        if self.Ana.REF is not None:
            self.register('res', self.Ana.REF.Cut.make_trk_residual(redo), 20, 'tracks with a small residual in the REF')
        if not hasattr(self.Ana, 'IsRef'):
            self.register('fid', self.make_fiducial(redo=redo), 30, 'tracks in fiducial area')
        self.register('mask', self.make_mask(_redo=redo), 31, 'masked pixels for tracks')
        self.register('tstart', self.make_start_time(_redo=redo), 40, 'exclude first events')
        self.register('chi2', self.make_chi2(_redo=redo), 50, f'chi2 < {100 * self.get_config("chi2 quantile", dtype=float)}%q')
        self.register('slope', self.make_slope(_redo=redo), 88, '{}%q < slope < {}%q'.format(*100 * abs(array([0, 1]) - self.get_config('slope quantile', dtype=float))))

    def make_trk(self, trks):
        return self.make_ev(trks, self.Ana.NTracks)

    @save_cut('TrackMask', cfg='track mask')
    def make_mask(self, t=None, _redo=False):
        return self.make_cluster_mask(*self.get_config('track mask', default=zeros((0, 2))).T, t=choose(t, self.get_config('track mask range', default=1.1)))

    @no_trans
    def make_cluster_mask(self, mx, my, t=.5):
        """cluster mask in track space doesn't require the DUT alignment"""  # no duplicate
        x, y = self.Ana.get_xy(local=True, cut=False)  # noqa
        return all([invert((x >= mx[i] - t) & (x <= mx[i] + t) & (y >= my[i] - t) & (y <= my[i] + t)) for i in range(mx.size)], axis=0)

