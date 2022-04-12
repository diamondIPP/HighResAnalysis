#!/usr/bin/env python
# --------------------------------------------------------
#       Residual calculations
# created on August 19th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from numpy import arctan, sqrt, array, quantile, mean, polyfit, identity, arange
from src.dut_analysis import DUTAnalysis, bins, prep_kw, partial, find_bins, FitRes, save_pickle
from plotting.fit import Gauss
from plotting.draw import np_profile, choose
from utility.affine_transformations import transform, m_transform, matrix, scale_matrix, inv


class ResidualAnalysis(DUTAnalysis):

    def __init__(self, parent, transformate=False):  # noqa

        self.Parent = parent
        self.__dict__.update(parent.__dict__)
        self.MetaSubDir = 'residuals'

        self.Rot, self.Trans = [], []

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

    def get_all(self, cut=None, pl=None):
        return [getattr(self, f'get_{i}')(cut, pl) for i in ['x', 'y', 'tx', 'ty']]
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
        return self._draw_angle(self.get_u(cut, pl), self.dv(cut, pl), prof, pl=pl, **prep_kw(dkw, title='X dY', x_tit='X [mm]', y_tit='dY [mm]'))

    def draw_vdu(self, cut=None, prof=True, pl=None, **dkw):
        return self._draw_angle(self.get_v(cut, pl), self.du(cut, pl), prof, pl=pl, xb=False, **prep_kw(dkw, title='Y dX', x_tit='Y [mm]', y_tit='dX [mm]'))
    # endregion DRAW
    # ----------------------------------------

    # ----------------------------------------
    # region ALIGN
    @property
    def m(self):
        return self.align()

    @save_pickle('AM', sub_dir='alignment', run='', suf_args='[6]')
    def align(self, d=None, m=None, cut=None, pl=None, p=.05, i=0, imax=20, _redo=False, _save=False):
        if d is None:
            self.PBar.start(imax)
        sx, sy = self.plane(pl).PX, self.plane(pl).PY
        x, y = transform(*self.get_xy(local=True, cut=cut, pl=pl), sx, sy) if d is None else d[:2]  # convert to mm
        tx, ty = transform(*self.get_txy(local=True, cut=cut, pl=pl, trans=False), sx, sy) if d is None else d[2:]
        d = self.rotate(x, y, *((tx, ty) if m is None else m_transform(m, *d[2:])), p=p)
        t = self.translate(*d[:4])
        self.Rot.append(d[4])
        self.Trans.append(t)
        m = choose(m, identity(3)) @ matrix(1, 1, *t, rx=d[4], order='str')
        self.PBar.update()
        if i < imax - 1:
            return self.align([x, y, tx, ty], m, p=p, i=i + 1, imax=imax, _redo=_redo)
        s = scale_matrix(sx, sy)
        return inv(s) @ m @ s

    def draw_align(self, imax=20):
        self.align(imax=imax, _redo=True)
        x = arange(imax)
        self.Draw.graph(x, self.Rot, x_tit='Iteration', y_tit='Rotation Angle [rad]')
        self.Draw.multigraph([self.Draw.graph(x, m, x_tit='Iteration', y_tit='Mean Residual [mm]', show=False) for m in array(self.Trans).T], 'Residuals', ['x', 'y'])

    @staticmethod
    def rotate(*d, p):
        x, y, tx, ty = d
        dx, dy = x - tx, y - ty
        qx, qy = quantile([dx, dy], [p, 1 - p], axis=1).T
        c = (dx > qx[0]) & (dx < qx[1]) & (dy > qy[0]) & (dy < qy[1])  # cut out the tails of the residuals in both dimensions
        x, y, tx, ty = x[c], y[c], tx[c], ty[c]
        r = mean([fit_angle(x, dy[c]), -1 * fit_angle(y, dx[c])])
        return x, y, *transform(tx, ty, rx=r), r

    @staticmethod
    def translate(*d):
        x, y, tx, ty = d
        return mean([x - tx, y - ty], axis=1)
    # endregion ALIGN
    # ----------------------------------------


def fit_angle(x, y):
    x, y, ey = np_profile(x, y, u=False)
    return arctan(polyfit(x, y, w=1 / ey, deg=1)[0])
