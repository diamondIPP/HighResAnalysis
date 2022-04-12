#!/usr/bin/env python
# --------------------------------------------------------
#       class for efficiency of a single DUT
# created on March 22nd 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from src.dut_analysis import DUTAnalysis, calc_eff, prep_kw, bins, get_2d_hist_vec, partial, mean, ax_range


class Efficiency(DUTAnalysis):

    def __init__(self, parent):  # noqa
        self.__dict__.update(parent.__dict__)
        self.MetaSubDir = 'Eff'
        self.Cut = self.Tracks.Cut

    def __call__(self):
        return self.get_value()

    @property
    def y_args(self):
        return {'y_tit': 'Efficiency [%]', 'y_range': [0, 105]}

    def get_values(self, cut=None):
        return (self.get_cluster_size(cut=cut) > 0).astype('u2') * 100

    def get_segment_values(self, nx=2, ny=3, cut=None):
        return get_2d_hist_vec(self.draw_map(local=True, cut=cut, binning=bins.make2d(*self.get_segments(nx, ny)), show=False), err=False, flat=False)

    def get_value(self, cut=None):
        return calc_eff(values=self.get_values(cut))

    def draw(self, bw=None, **dkw):
        t, e = self.get_time(), self.get_values() / 100
        return self.Draw.efficiency(t, e, **prep_kw(dkw, w=bw, **self.t_args(), **self.y_args, stats=False))

    def draw_vs_tp(self, **dkw):
        x, y = [f(cut=self.Cut.exclude('tp')) for f in [self.get_trigger_phase, self.get_values]]  # noqa
        return self.Draw.efficiency(x, y / 100, bins.TP, 'Efficiency vs. Trigger Phase', **prep_kw(dkw, x_tit='Trigger Phase', **self.y_args))

    def draw_map(self, res=.5, local=True, eff=True, both=False, fid=False, cut=None, **dkw):
        (x, y), e = [f(cut=self.Cut.get_nofid(cut, fid)) for f in [partial(self.Tracks.get_xy, local=local), self.get_values]]
        p = self.Draw.prof2d(x, y, e, **prep_kw(dkw, title='Efficiency Map', binning=bins.get_xy(local, self.Plane, res), **self.ax_tits(local)))
        self.draw_text(self.Surface, cut, eff)
        self.draw_text(not self.Surface, cut, eff and both)
        return p

    def draw_text(self, surface, cut, show=True):
        if show:
            self.activate_surface(surface)
            self.Cut.get_fid(surface, show=show)
            x, y = [.5 if i is None else mean(ax_range(i)) for i in self.Cut.get_fid_config(surface)]
            self.Draw.textbox(f'{self.get_value(cut)[0]:2.1f}%', x, y, s=2, opacity=.4, ndc=False)

    def draw_segments(self, res=.5, local=True, nx=10, ny=15, cut=None, **dkw):
        e = self.get_segment_values(nx, ny, cut)
        self.draw_map(res, local, cut=cut, eff=False, **dkw)
        xs, ys = self.get_segment_centres(nx, ny)
        [self.Draw.tlatex(x, y, '{:2.1f}'.format(e[i][j]), str(i * y.size + j), 22, size=.02) for j, x in enumerate(xs) for i, y in enumerate(ys)]
        self.draw_grid(nx, ny)

    def draw_segment_distribution(self, nx=10, ny=15, cut=None, segments=True, **dkw):
        e = self.get_segment_values(nx, ny, cut).flatten() if segments else get_2d_hist_vec(self.draw_map(.5, show=False), err=False)
        self.Draw.distribution(e, **prep_kw(dkw, title='Segment Efficiencies', x_tit='Efficiency [%]'))

    def draw_in_pixel(self, n=10, ox=0, oy=0, cut=None, **dkw):
        return super().draw_in_pixel(ox, oy, n, cut, self.get_values, tit='Efficiency', **dkw)

    def draw_in_cell(self, n=10, ox=0, oy=0, cut=None, **dkw):
        return super().draw_in_cell(ox, oy, n, cut, self.get_values, tit='Efficiency', **dkw)
