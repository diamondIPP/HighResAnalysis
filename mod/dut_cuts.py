#!/usr/bin/env python
# --------------------------------------------------------
#       cuts for analysis of a single DUT
# created on March 26th 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from numpy import array, invert, all, zeros, quantile, max, inf, sqrt

from plotting.draw import make_box_args, Draw, prep_kw, TCutG, Config
from src.cut import Cuts
from utility.utils import Dir, save_hdf5, parallel, make_list, choose


class DUTCut(Cuts):
    def __init__(self, ana, meta_sub_dir='dut_cuts'):
        self.Ana = ana
        self.MetaSubDir = meta_sub_dir
        self.N = self.Ana.N
        super().__init__()

    def __call__(self, cut=None, data=None, pl=None):
        cut = super().__call__(cut)
        if data is None:
            return cut
        if data.size == self.Ana.NTracks:
            return self.trk2pl(data, pl)[cut]
        if data.size == self.Ana.NEvents:
            return self.ev2pl(data, pl)[cut] if cut is ... or cut.size == self.Ana.N else self.ev2trk(data)[cut]
        return data if cut is ... else data[cut] if cut.size == data.size else data[self.trk2pl(cut, pl)]

    def register(self, name, values=None, level=None, description=None):
        """ guarantee that the cut array has the correct size """
        return super().register(name, self(cut=False, data=values), level, description)

    def init_config(self):
        return Config(Dir.joinpath('cuts', f'cut{self.Ana.BeamTest.Tag}.ini'), section=self.Ana.DUT.Name)

    def make(self, redo=False):
        self.register('fid', self.make_fiducial(_redo=redo), 10, 'fid cut')
        self.register('mask', self.make_cluster_mask(), 20, 'mask pixels')
        self.register('charge', self.Ana.get_phs(cut=False) != 0, 30, 'events with non-zero charge')
        self.register('cluster', self.make_cluster(_redo=redo), 90, 'tracks with a cluster')

    def make_additional(self, redo=False):
        self.register('res', self.make_ref_residual(redo=redo), 69, 'small residuals to REF plane')
        self.register('tp', self.make_trigger_phase(_redo=redo), 61, 'trigger phase')
        self.register('tstart', self.make_start_time(_redo=redo), 40, 'exclude first events')
        self.register('chi2', self.make_chi2(_redo=redo), 50, 'small chi2')

    # ----------------------------------------
    # region GENERATE
    @save_hdf5('Fid', arr=True, dtype='?', suf_args='all')
    @parallel('point_in_polygon', 'fiducial cut')
    def make_fiducial(self, surface=False, _redo=False):
        x, y = self.Ana.get_xy(local=True, cut=False)
        return array([x, y]).T, self.get_fid(surface=surface)

    def make_cluster_mask(self):
        """ exclude all clusters within half a pixel distance of the masked pixel"""
        x, y = self.Ana.get_xy(local=True, cut=False)
        mx, my = array(self.get_config('mask')).T
        return all([invert((x >= mx[i] - .5) & (x <= mx[i] + .5) & (y >= my[i] - .5) & (y <= my[i] + .5)) for i in range(mx.size)], axis=0)

    @save_hdf5('TP', arr=True, dtype='?')
    def make_trigger_phase(self, _redo=False):
        tp = self.Ana.get_trigger_phase(cut=False)
        low, high = self.get_config('trigger phase')
        return (tp >= low) & (tp <= high)

    @save_hdf5('Clu', arr=True, dtype='?', suf_args='all')
    def make_cluster(self, pl=None, _redo=False):
        return self.Ana.get_data('Clusters', 'Size', cut=False, pl=pl) > 0

    @save_hdf5('Time', arr=True, dtype='?')
    def make_start_time(self, _redo=False):
        t = self.Ana.get_time(cut=False)
        return t >= t[0] + self.get_config('start time', dtype=int) * 60

    @save_hdf5('Chi2', arr=True, dtype='?', suf_args='all')
    def make_chi2(self, q=None, _redo=False):
        x, q = self.Ana.get_chi2(cut=False), choose(q, self.get_config('chi2 quantile', dtype=float))
        return x < quantile(x, q)

    @save_hdf5('Res', arr=True, dtype='?', suf_args='all')
    def make_residual(self, v=None, _redo=False):
        x, y, (mx, my) = self.Ana.Residuals.du(cut=0), self.Ana.Residuals.dv(cut=0), self.Ana.Residuals.get_means(cut=0)
        return sqrt((x - mx.n) ** 2 + (y - my.n) ** 2) < choose(v, self.get_config('residuals', dtype=float))

    def make_ref_residual(self, pl=None, redo=False):
        return self.trk2pl(self.Ana.REF.Cut.make_trk_residual(redo), pl)

    def make_ph(self, xmax, xmin=None):
        x = self.Ana.get_phs(cut=False)
        return (x >= choose(xmin, -inf)) & (x < xmax)

    def make_correlation(self, pl0, pl1=None):
        return self.make_cluster(pl0) & self.make_cluster(pl1)
    # endregion GENERATE
    # ----------------------------------------

    # ----------------------------------------
    # region FIDUCIAL
    def get_fid_config(self, surface=False, name=None):
        p = self.get_config(choose(name, 'surface fiducial' if surface else 'fiducial'))
        p = make_box_args(*p[[0, 2, 1, 3]]) if p.size == 4 else p  # unpack short box notation
        p[p == max(p, axis=1).reshape((-1, 1))] += 1  # extend one pixel to the top and right
        return p - .5  # pixel centre is at the integer

    def get_fid(self, surface=False, name=None, **dkw):
        x, y = self.get_fid_config(surface, name)
        return Draw.polygon(x, y, **prep_kw(dkw, show=False, line_color=2, width=2, name=choose(name, f'fid{surface:d}')))

    def get_full_fid(self):
        return self.get_fid(name='full size')

    def draw_fid(self, surface=False, name=None, **dkw):
        self.get_fid(surface, name, **prep_kw(dkw, show=True))

    @staticmethod
    def point_in_polygon(p, poly: TCutG):
        return poly.IsInside(*p)
    # endregion FIDUCIAL
    # ----------------------------------------

    def get_track_events(self):
        return array(self.Ana.F['Tracks']['EvtFrame'])

    def make_ev(self, ev, n=None):
        c = zeros(choose(n, self.Ana.NEvents), '?')
        c[make_list(ev).astype('i')] = True
        return c

    def ev2pl(self, x, pl=None):
        return self.trk2pl(self.ev2trk(x), pl)

    def dut2ev(self, x):
        c = zeros(self.Ana.NEvents, '?')
        c[self.get_track_events()[self.pl2trk(x)]] = True
        return c

    def pl2trk(self, x, pl=None):
        c = zeros(self.Ana.NTracks, '?')
        c[self.make_cluster(pl)] = True if self(x) is ... else self(x)
        return c

    def ev2trk(self, x):
        return x[self.get_track_events()]

    def trk2ev(self, x):
        c = zeros(self.Ana.NEvents, '?')
        c[self.get_track_events()[x]] = True
        return c

    def trk2pl(self, x, pl=None):
        return x[self.make_cluster(pl)]

    @staticmethod
    def to_trk(cut):
        return -1 if type(cut) is bool or cut is ... else cut

    def make_hdf5_path(self, *args, **kwargs):
        return self.Ana.make_hdf5_path(*args, **prep_kw(kwargs, sub_dir=self.MetaSubDir))

    def info(self, *args, **kwargs):
        return self.Ana.info(*args, **kwargs)
