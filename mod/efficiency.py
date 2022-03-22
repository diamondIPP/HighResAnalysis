#!/usr/bin/env python
# --------------------------------------------------------
#       class for efficiency of a single DUT
# created on March 22nd 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from src.dut_analysis import DUTAnalysis, calc_eff, prep_kw, find_bins, bins, get_2d_hist_vec, update_canvas, partial


class Efficiency(DUTAnalysis):

    def __init__(self, parent):  # noqa
        self.__dict__.update(parent.__dict__)
        self.MetaSubDir = 'Efficiency'

    def __call__(self):
        return self.get_value()

    @property
    def y_args(self):
        return {'y_tit': 'Efficiency [%]', 'y_range': [0, 105]}

    def get_values(self, trk_cut=None):
        return (self.get_cluster_size(trk_cut=trk_cut) > 0).astype('u2') * 100

    def get_segment_values(self, nx=2, ny=3, cut=None):
        return get_2d_hist_vec(self.draw_map(local=True, cut=cut, binning=bins.make2d(*self.get_segments(nx, ny)), show=False), err=False, flat=False)

    def get_value(self, trk_cut=None):
        return calc_eff(values=self.get_values(trk_cut))

    def draw(self, bw=None, **dkw):
        t, e = self.get_time(trk_cut=None), self.get_values() / 100
        return self.Draw.efficiency(t, e, find_bins(t, w=bw), **prep_kw(dkw, **self.t_args(), **self.y_args, stats=False))

    def draw_vs_tp(self, **dkw):
        x, y = [f(trk_cut=self.Tracks.Cuts.exclude('triggerphase')) for f in [self.get_trigger_phase, self.get_values]]
        return self.Draw.profile(x, y, bins.TP, 'Efficiency vs. Trigger Phase', **prep_kw(dkw, x_tit='Trigger Phase', **self.y_args))

    def draw_map(self, res=.5, local=True, eff=True, both=False, fid=False, cut=None, **dkw):
        (x, y), e = [f(trk_cut=self.Tracks.Cuts.get_nofid(cut, fid)) for f in [partial(self.Tracks.get_xy, local=local), self.get_values]]
        p = self.Draw.prof2d(x, y, e, **prep_kw(dkw, title='Efficiency Map', binning=bins.get_xy(local, self.Plane, res), **self.get_ax_tits(local), leg=self.get_fid() if eff else None))
        self.draw_text(self.Surface, cut, eff)
        self.draw_text(not self.Surface, cut, eff and both)
        return p

    def draw_text(self, surface, cut, show=True):
        if show:
            self.activate_surface(surface)
            self.get_fid(show=show)
            x, y = self.Cuts.get_fid_config(surface)
            self.Draw.tlatex(x[0] + (x[2] - x[0]) / 2, y[0] + (y[1] - y[0]) / 2, '{:2.1f}%'.format(self.get_value(cut)[0]), 'eff', 22, size=.04)

    def draw_segments(self, res=.5, local=True, nx=10, ny=15, cut=None, **dkw):
        e = self.get_segment_values(nx, ny, cut)
        self.draw_map(res, local, cut=cut, eff=False, **dkw)
        xs, ys = self.get_segment_centres(nx, ny)
        print(xs.size, ys.size)
        [self.Draw.tlatex(x, y, '{:2.1f}'.format(e[i][j]), str(i * y.size + j), 22, size=.02) for j, x in enumerate(xs) for i, y in enumerate(ys)]
        self.draw_grid(nx, ny)

    def draw_segment_distribution(self, nx=10, ny=15, cut=None, segments=True, **dkw):
        e = self.get_segment_values(nx, ny, cut).flatten() if segments else get_2d_hist_vec(self.draw_map(.5, show=False), err=False)
        self.Draw.distribution(e, **prep_kw(dkw, title='Segment Efficiencies', x_tit='Efficiency [%]'))

    def draw_in_pixel(self, res=.1, cut=None, show=True, cell=False):
        (x, y), e = self.Tracks.get_xy(trk_cut=cut), self.get_values(cut)
        x, y, e = self.expand_inpixel(x, y, e, cell)
        self.Draw.prof2d(x, y, e, bins.get_pixel(self.Plane, res, cell=cell), 'Efficiency Map in {}'.format('3D Cell' if cell else 'Pixel'), show=show, stats=0)
        self.Draw.box(0, 0, 1, 1)
        update_canvas()
