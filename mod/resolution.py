#!/usr/bin/env python
# --------------------------------------------------------
#       Analysis of the reference plane
# created on August 19th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from src.dut_analysis import DUTAnalysis, bins, update_pbar, linspace, ufloat, tile, save_pickle, array
from plotting.draw import mean_sigma, prep_kw, set_statbox
from uncertainties.umath import sqrt as usqrt  # noqa
from mod.reso_cuts import ResCut
from mod.residuals import ResidualAnalysis


class Resolution(DUTAnalysis):

    def __init__(self, parent):  # noqa

        self.Parent = parent
        self.__dict__.update(parent.__dict__)
        self.Cut = ResCut(self.Cut)
        self.Residuals = ResidualAnalysis(self)
        self.MetaSubDir = 'res'

    def __call__(self, redo=False):
        return self.x(_redo=redo), self.y(_redo=redo)

    @save_pickle('X', suf_args='all')
    @update_pbar
    def x(self, chi2=None, slope=None, _redo=False):
        self.Cut.set_chi2(chi2)
        self.Cut.set_slope(slope)
        return usqrt(max(ufloat(1 / 12, 0), mean_sigma(self.Residuals.dx())[1] ** 2) - 1 / 12) * self.Plane.PX * 1e3

    @save_pickle('Y', suf_args='all')
    @update_pbar
    def y(self, chi2=None, slope=None, _redo=False):
        self.Cut.set_chi2(chi2)
        self.Cut.set_slope(slope)
        return usqrt(mean_sigma(self.Residuals.dy())[1] ** 2 - 1 / 12) * self.Plane.PY * 1e3

    def draw_x(self, **dkw):
        return self.Residuals.draw_x(**prep_kw(dkw, lf=.5, rf=.5, stats=set_statbox(all_stat=True, form='.2f')))

    def draw_y(self, **dkw):
        return self.Residuals.draw_y(**prep_kw(dkw, lf=.5, rf=.5, stats=set_statbox(all_stat=True, form='.2f')))

    def draw_cross(self, bw=.25, col=20, row=40, **dkw):
        self.draw_hit_map(bw, cut=self.Cut.make_cross(col, row), **prep_kw(dkw, leg=self.Plane.get_grid()))

    def draw_chi2_cross(self, bw=.25, col=20, row=40, **dkw):
        (x, y), zz = [f(cut=self.Cut.make_cross(col, row)) for f in [self.get_txy, self.get_chi2]]
        self.Draw.prof2d(x, y, zz, bins.get_local(self.Plane, bw), **prep_kw(dkw, **self.ax_tits(), z_tit='#chi^{2}'))

    def draw_chi_x(self, col=20, **dkw):
        x, y = [f(cut=self.Cut.make_col(col)) for f in [self.Residuals.tx, self.get_chi2]]  # noqa
        self.Draw.profile(x, y, **prep_kw(dkw))

    def draw_slope_x(self, col=20, **dkw):
        x, y = [f(cut=self.Cut.make_col(col)) for f in [self.Residuals.tx, self.get_slope_x]]  # noqa
        self.Draw.profile(x, y, **prep_kw(dkw))

    def draw_x_vs_chi2(self, n=20):
        x = linspace(1 / n, 1, n)
        self.PBar.start(n)
        self.Draw.graph(x, [self.x(i) for i in x])

    def draw_x_vs_slope(self, n=20):
        x = linspace(0, .5 - .5 / n, n)
        self.PBar.start(n)
        self.Draw.graph(x, [self.x(slope=i) for i in x])

    def draw_x_vs_pars(self, n=10):
        x, y = linspace(.5 / n, .5, n), linspace(.1, .4, n)
        self.PBar.start(n ** 2)
        d = array([[self.x(i, j) for j in y] for i in x])
        self.Draw.prof2d(x.repeat(10), tile(y, 10), d.flatten(), bins.from_p(x) + bins.from_p(y))
        return d
