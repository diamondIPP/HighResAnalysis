#!/usr/bin/env python
# --------------------------------------------------------
#       cuts for analysis of a single DUT
# created on March 26th 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from numpy import array, invert, all, zeros, quantile, max, inf, sqrt, where, ndarray, any, append, ones

from plotting.draw import make_box_args, Draw, prep_kw, TCutG, Config
from src.cut import Cuts
from utility.utils import Dir, save_hdf5, parallel, make_list, choose, save_pickle, uarr2n


def save_cut(*pargs, suf_args='[]', field=None, verbose=False, cfg=None, cfield=None, **pkwargs):
    def inner(f):
        def wrapper(*args, **kwargs):
            ana = args[0]
            redo, kw_redo = False, kwargs.pop('_redo') if '_redo' in kwargs else False
            if cfg is not None:
                @save_pickle(cfg.replace(' ', ''), field=cfield, run='', sub_dir='cuts', verbose=verbose)
                def get_config(a, c, _redo=False):
                    return a.get_config(c)
                redo = ana.get_config(cfg) != get_config(ana, cfg)
                redo = any(redo) if type(redo) is ndarray else redo
                get_config(ana, cfg, _redo=redo)
            return save_hdf5(*pargs, suf_args=suf_args, arr=True, dtype='?', field=field, verbose=verbose, **pkwargs)(f)(*args, _redo=redo or kw_redo, **kwargs)
        return wrapper
    return inner


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
        return None if values is None else super().register(name, self(cut=False, data=values), level, description)

    def init_config(self):
        return Config(Dir.joinpath('cuts', f'cut{self.Ana.BeamTest.Tag}.ini'), section=self.Ana.DUT.Name)

    def make(self, redo=False):
        self.register('fid', self.make_fiducial(redo=redo), 10, 'fid cut')
        self.register('mask', self.make_conf_mask(), 20, 'mask pixels')
        self.register('tmask', self.make_cal_thresh_mask(), 21, 'mask pixels with high treshold')
        self.register('cmask', self.make_cal_chi2_mask(), 22, f'mask pixels with calibration fit chi2 > {self.get_config("calibration chi2", default=10.)}')
        self.register('charge', self.Ana.get_phs(cut=False) != 0, 30, 'events with non-zero charge')

    def make_additional(self, redo=False):
        self.register('res', self.make_ref_residual(redo=redo), 60, 'small residuals to REF plane')
        self.register('tp', self.make_trigger_phase(_redo=redo), 40, 'trigger phase')
        self.register('tstart', self.make_start_time(_redo=redo), 35, 'exclude first events')
        self.register('chi2', self.make_chi2(_redo=redo), 50, f'small chi2 < q({self.get_config("chi2 quantile", dtype=float)})')
        self.register('slope', self.make_slope(_redo=redo), 55, f'straight tracks < q({self.get_config("slope quantile", dtype=float)})')

    # ----------------------------------------
    # region GENERATE
    @save_cut('Fid', suf_args='all', cfg='fiducial', cfield='MetaSubDir')
    @parallel('point_in_polygon', 'fiducial cut')
    def make_fiducial_(self, surface=False, _redo=False):
        x, y = self.Ana.get_xy(local=True, cut=False)
        return array([x, y]).T, self.get_fid(surface=surface)

    def make_fiducial(self, surface=False, redo=False):
        return None if self.get_fid(surface=surface) is None else self.make_fiducial_(surface, _redo=redo)

    def make_cluster_mask(self, mx, my, t=.5):
        """ exclude all clusters within half a pixel distance of the masked pixel"""
        x, y = self.Ana.get_xy(local=True, cut=False)
        return all([invert((x >= mx[i] - t) & (x <= mx[i] + t) & (y >= my[i] - t) & (y <= my[i] + t)) for i in range(mx.size)], axis=0)

    @save_cut('Mask', cfg='mask')
    def make_conf_mask(self, _redo=False):
        return self.make_cluster_mask(*self.get_config('mask', default=zeros((0, 2))).T)

    def get_thresh_mask(self):
        return where(self.Ana.Calibration.get_thresholds() > self.Ana.Calibration.Trim * 1.5)

    @save_cut('TMask')
    def make_cal_thresh_mask(self, _redo=False):
        return self.make_cluster_mask(*self.get_thresh_mask())

    def get_cal_chi2_mask(self):
        return where(self.Ana.Calibration.get_chi2s() > self.get_config('calibration chi2', default=10.))

    @save_cut('CMask', cfg='calibration chi2')
    def make_cal_chi2_mask(self, _redo=False):
        return self.make_cluster_mask(*self.get_cal_chi2_mask())

    @save_cut('TP', cfg='trigger phase')
    def make_trigger_phase(self, _redo=False):
        tp = self.Ana.get_trigger_phase(cut=False)
        low, high = self.get_config('trigger phase')
        return (tp >= low) & (tp <= high)

    @save_cut('Clu', suf_args='all')
    def make_cluster(self, pl=None, _redo=False):
        return self.Ana.get_cluster_size(cut=False, pl=pl) > 0

    def make_cs(self, n=1, pl=None):
        return self.Ana.get_cluster_size(cut=0, pl=pl) == n

    @save_cut('Time', cfg='start time')
    def make_start_time(self, _redo=False):
        t = self.Ana.get_time(cut=False)
        return t >= t[0] + self.get_config('start time', dtype=int) * 60

    @save_cut('Chi2', suf_args='all', cfg='chi2 quantile')
    def make_chi2(self, q=None, _redo=False, _save=True):
        x, q = self.Ana.get_chi2(cut=False), choose(q, self.get_config('chi2 quantile', dtype=float))
        return ones(x.size, '?') if q == 1 else x < quantile(x, q)

    @save_cut('Res', suf_args='all', cfg='residuals')
    def make_residual(self, v=None, _redo=False):
        x, y, (mx, my) = self.Ana.Residuals.du(cut=0), self.Ana.Residuals.dv(cut=0), self.Ana.Residuals.get_means(cut=0)
        return sqrt((x - mx.n) ** 2 + (y - my.n) ** 2) < choose(v, self.get_config('residuals', dtype=float))

    def make_trk_residual(self, redo=False):
        return self.pl2trk(self.make_residual(_redo=redo))

    def make_ref_residual(self, pl=None, redo=False):
        return self.trk2pl(self.Ana.REF.Cut.make_trk_residual(redo), pl)

    def make_ph(self, xmax, xmin=None):
        x = self.Ana.get_phs(cut=False)
        return (x >= choose(xmin, -inf)) & (x < xmax)

    def make_correlation(self, pl0, pl1=None):
        return self.make_cluster(pl0) & self.make_cluster(pl1)

    @save_cut('Slope', suf_args='all', cfg='slope quantile')
    def make_slope(self, q=None, _redo=False, _save=True):
        x, y, q = self.Ana.get_slope_x(cut=False), self.Ana.get_slope_y(cut=False), choose(q, self.get_config('slope quantile', dtype=float))
        (xmin, xmax), (ymin, ymax) = quantile(x, [q, 1 - q]), quantile(y, [q, 1 - q])
        return ones(x.size, '?') if q >= .5 else (x > xmin) & (x < xmax) & (y > ymin) & (y < ymax)
    # endregion GENERATE
    # ----------------------------------------

    # ----------------------------------------
    # region SET
    def set(self, name, values=None):
        if values is not None:
            super().set(name, self(cut=False, data=values))

    def set_chi2(self, q):
        self.set('chi2', self.make_chi2(q, _save=False))

    def set_slope(self, q):
        self.set('slope', self.make_slope(q, _save=False))
    # endregion SET
    # ----------------------------------------

    # ----------------------------------------
    # region FIDUCIAL
    def get_fid_config(self, surface=False, name=None):
        p = self.get_config(choose(name, 'surface fiducial' if surface else 'fiducial'), default=self.get_config('full size'))
        if p is None:
            return None, None
        p = make_box_args(*p[[0, 2, 1, 3]]) if p.size == 4 else p  # unpack short box notation
        p[p == max(p, axis=1).reshape((-1, 1))] += 1  # extend one pixel to the top and right
        return p - .5  # pixel centre is at the integer

    def get_fid(self, surface=False, name=None, **dkw):
        x, y = self.get_fid_config(surface, name)
        return None if x is None else Draw.polygon(x, y, **prep_kw(dkw, show=False, line_color=2, width=2, name=choose(name, f'fid{surface:d}')))

    def get_full_fid(self):
        return self.get_fid(name='full size')

    def draw_fid(self, surface=False, name=None, **dkw):
        self.get_fid(surface, name, **prep_kw(dkw, show=True))

    @staticmethod
    def point_in_polygon(p, poly: TCutG):
        return poly.IsInside(*p)
    # endregion FIDUCIAL
    # ----------------------------------------

    def get_res(self, **dkw):
        r, mx, my = append(self.get_config('residuals', dtype=float), uarr2n(self.Ana.Residuals.get_means(cut=0))) * 1e3
        return Draw.circle(r, mx, my, **prep_kw(dkw, show=False))

    def draw_res(self, **dkw):
        self.get_res(**prep_kw(dkw, show=True))

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

    def make_pickle_path(self, *args, **kwargs):
        return self.Ana.make_pickle_path(*args, **prep_kw(kwargs, sub_dir=self.MetaSubDir))

    def info(self, *args, **kwargs):
        return self.Ana.info(*args, **kwargs)
