#!/usr/bin/env python
# --------------------------------------------------------
#       Residual calculations
# created on August 19th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from numpy import arctan, sqrt, array
from src.dut_analysis import DUTAnalysis, bins, prep_kw, partial, find_bins, FitRes, save_pickle
from plotting.fit import Gauss


class ResidualAnalysis(DUTAnalysis):

    def __init__(self, parent, transformate=False):  # noqa

        self.Parent = parent
        self.__dict__.update(parent.__dict__)
        self.MetaSubDir = 'Res'

    # ----------------------------------------
    # region DATA
    def du(self, cut=None, pl=None):
        return self.get_u(cut, pl) - self.get_tu(cut, pl)

    def dv(self, cut=None, pl=None):
        return self.get_v(cut, pl) - self.get_tv(cut, pl)

    def dx(self, cut=None, pl=None):
        return self.get_x(cut, pl) - self.get_tx(cut, pl)

    def dy(self, cut=None, pl=None):
        return self.get_y(cut, pl) - self.get_ty(cut, pl)

    def __call__(self, cut=None, pl=None):
        return sqrt(self.du(cut, pl) ** 2 + self.dv(cut, pl) ** 2)

    def du_(self, cut=None, pl=None):
        return self.get_data('Tracks', 'dU', cut, pl) * 1e3  # still not clear what this dU is... error on U of the track?

    def dv_(self, cut=None, pl=None):
        return self.get_data('Tracks', 'dV', cut, pl) * 1e3
    # endregion DATA
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    @save_pickle('RM', suf_args='all')
    def get_means(self, local=False, cut=None, pl=None, _redo=False):
        return array([Gauss(f(show=False, cut=cut, pl=pl), thresh=.05).fit()[1] for f in ([self.draw_x, self.draw_y] if local else [self.draw_u, self.draw_v])]) / (1 if local else 1e3)
    # endregion DATA
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw_u(self, cut=None, pl=None, **dkw):
        return self.Draw.distribution(self.du(cut, pl) * 1e3, **prep_kw(dkw, r=[-300, 300], title='X Residuals', x_tit='Residual [#mum]'))

    def draw_v(self, cut=None, pl=None, **dkw):
        return self.Draw.distribution(self.dv(cut, pl) * 1e3, **prep_kw(dkw, r=[-300, 300], title='Y Residuals', x_tit='Residual [#mum]'))

    def draw_x(self, cut=None, pl=None, **dkw):
        return self.Draw.distribution(self.dx(cut, pl), **prep_kw(dkw, title='X Residuals', x_tit='Residual [Columns]'))

    def draw_y(self, cut=None, pl=None, **dkw):
        return self.Draw.distribution(self.dy(cut, pl), **prep_kw(dkw, title='Y Residuals', x_tit='Residual [Rows]'))

    def draw_xy(self, bw=10, cut=None, pl=None, **dkw):
        x, y = array([f(cut=self.Cut.exclude('res', cut), pl=pl) for f in [self.du, self.dv]]) * 1e3  # noqa
        self.Draw.histo_2d(x, y, bins.make(-1000, 1000, bw) * 2, **prep_kw(dkw, title='XY Residual', x_tit='dX [#mum]', y_tit='dY [#mum]'))

    def draw(self, bw=10, pl=None, **dkw):
        self.Draw.distribution(self(cut=False, pl=pl) * 1e3, bins.make(0, 1000, bw), **prep_kw(dkw, title='Residuals', x_tit='Residual [#mum]'))

    def draw_map(self, res=.3, local=True, cut=None, fid=False, pl=None, **dkw):
        (x, y), z_ = [f(cut=self.Cut.get_nofid(fid, cut), pl=pl) for f in [partial(self.get_txy, local=local), self]]
        self.Draw.prof2d(x, y, z_ * 1e3, bins.get_xy(local, self.Plane, res), 'Residuals', **prep_kw(dkw, z_tit='Residuals [#mum]', **self.ax_tits(local)))

    def _draw_angle(self, x, y, prof=False, xb=True, local=False, pl=None, **dkw):
        b = (bins.get_x if xb else bins.get_y)(self.plane(pl), local=local) + find_bins(y)
        return arctan(FitRes((self.Draw.profile if prof else self.Draw.histo_2d)(x, y, b[:2 if prof else 4], graph=True, **dkw).Fit('pol1', 'qs'))[1].n)

    def draw_udv(self, cut=None, prof=True, pl=None, **dkw):
        return self._draw_angle(self.get_u(cut, pl, self.T), self.dv(cut, pl), prof, pl=pl, **prep_kw(dkw, title='X dY', x_tit='X [mm]', y_tit='dY [mm]'))

    def draw_vdu(self, cut=None, prof=True, pl=None, **dkw):
        return self._draw_angle(self.get_v(cut, pl, self.T), self.du(cut, pl), prof, pl=pl, xb=False, **prep_kw(dkw, title='Y dX', x_tit='Y [mm]', y_tit='dX [mm]'))
    # endregion DRAW
    # ----------------------------------------
