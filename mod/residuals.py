#!/usr/bin/env python
# --------------------------------------------------------
#       Residual calculations
# created on August 19th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from numpy import arctan, sqrt, array, quantile, mean, polyfit, identity, arange
from src.dut_analysis import DUTAnalysis, bins, prep_kw, partial, find_bins, FitRes, save_pickle, no_trans
from plotting.fit import Gauss
from plotting.draw import np_profile, choose, set_x_range, ax_range, ufloat
from utility.utils import PBAR, uarr2n
from utility.affine_transformations import transform, m_transform, matrix, scale_matrix, inv


def res_analysis(cls):
    class ResidualAnalysis(cls):

        def __init__(self, parent: DUTAnalysis, transformate=False):  # noqa

            self.Parent = parent
            self.__dict__.update(parent.__dict__)
            self.MetaSubDir = 'residuals'
            self.Residuals = self  # for self reference

            self.Rot, self.Trans = [], []

        # ----------------------------------------
        # region DATA
        def du(self, cut=None, pl=None):
            return self.get_u(cut, pl) - self.get_tu(cut, pl)

        def dv(self, cut=None, pl=None):
            return self.get_v(cut, pl) - self.get_tv(cut, pl)

        def dx(self, cut=None, pl=None):
            return self.dxy(cut, pl, local=True)[0]

        def dy(self, cut=None, pl=None):
            return self.dxy(cut, pl, local=True)[1]

        def tx(self, cut=None, pl=None):
            return self.txy(cut, pl, local=True)[0]

        def ty(self, cut=None, pl=None):
            return self.txy(cut, pl, local=True)[1]

        def txy(self, cut=None, pl=None, local=True):
            return self.get_txy(local, cut, pl)

        def dxy(self, cut=None, pl=None, local=True):
            (x, y), (tx, ty) = self.get_xy(local, cut, pl), self.txy(cut, pl, local)
            return x - tx, y - ty

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
        @save_pickle('ResFit', suf_args='all', field='DUT')
        def fit(self, local=False, cut=None, pl=None, _redo=False):
            """Fit the residual distributions with a Gaussian. The distributions are usually not Gaussian shaped though..."""
            return array([Gauss(f(show=False, cut=cut, pl=pl), thresh=.05).fit()[1:] for f in ([self.draw_x, self.draw_y] if local else [self.draw_u, self.draw_v])]) / (1 if local else 1e3)

        @staticmethod
        def mean_std(f, local=False, cut=None, pl=None, thresh=.3):
            """estimate mean and std from the residual distribution"""
            h = f(show=False, cut=cut, pl=pl, w=5, save=False)
            set_x_range(*ax_range(h=h, thresh=thresh * h.GetMaximum()))
            return array([ufloat(h.GetMean(), h.GetMeanError()), ufloat(h.GetStdDev(), h.GetStdDevError())]) / (1 if local else 1e3)

        @save_pickle('Res', suf_args='all', field='DUT')
        def mean_stds(self, local=False, cut=None, pl=None, thresh=.3, _redo=False):
            return array([self.mean_std(f, local, cut, pl, thresh) for f in ([self.draw_x, self.draw_y] if local else [self.draw_u, self.draw_v])])

        def means(self, local=False, cut=None, pl=None, redo=False):
            return self.mean_stds(local, cut, pl, _redo=redo)[:, 0]

        def std(self, local=False, cut=None, pl=None, redo=False):
            return self.mean_stds(local, cut, pl, _redo=redo)[:, 1]

        def show(self, redo=False):
            x, y = self.mean_stds(_redo=redo) * 1e3
            print(f'Res X: mean = {x[0]:.1f}, std = {x[1]:.1f}')
            print(f'Res Y: mean = {y[0]:.1f}, std = {y[1]:.1f}')
        # endregion DATA
        # ----------------------------------------

        # ----------------------------------------
        # region DRAW
        def draw_u(self, cut=None, pl=None, **dkw):
            return self.Draw.distribution(self.du(cut, pl) * 1e3, **prep_kw(dkw, w=1, r=[-300, 300], title='X Residuals', x_tit='Residual [#mum]', file_name='ResX'))

        def draw_v(self, cut=None, pl=None, **dkw):
            return self.Draw.distribution(self.dv(cut, pl) * 1e3, **prep_kw(dkw, w=1, r=[-300, 300], title='Y Residuals', x_tit='Residual [#mum]', file_name='ResY'))

        def draw_x(self, cut=None, pl=None, **dkw):
            return self.Draw.distribution(self.dx(cut, pl), **prep_kw(dkw, title='X Residuals', x_tit='Residual [Columns]'))

        def draw_y(self, cut=None, pl=None, **dkw):
            return self.Draw.distribution(self.dy(cut, pl), **prep_kw(dkw, title='Y Residuals', x_tit='Residual [Rows]'))

        def draw_uv(self, bw=10, cut=None, pl=None, **dkw):
            return self.draw_xy(bw, cut, pl, **dkw)

        def draw_xy(self, bw=10, cut=None, pl=None, **dkw):
            x, y = array([f(cut=self.Cut.exclude('res', cut), pl=pl) for f in [self.du, self.dv]]) * 1e3  # noqa
            return self.Draw.histo_2d(x, y, bins.make(-1000, 1000, bw) * 2, **prep_kw(dkw, title='XY Residual', x_tit='dX [#mum]', y_tit='dY [#mum]', leg=self.cuts(), file_name=f'XYRes{self.Plane}'))

        def draw(self, bw=10, pl=None, **dkw):
            self.Draw.distribution(self(cut=..., pl=pl) * 1e3, bins.make(0, 1000, bw), **prep_kw(dkw, title='Residuals', x_tit='Residual [#mum]'))

        def draw_x_map(self, res=.3, cut=None, fid=False, pl=None, **dkw):
            (x, y), z_ = [f(cut=self.Cut.get_nofid(cut, fid), pl=pl) for f in [self.txy, self.dx]]  # noqa
            self.Draw.prof2d(x, y, z_, bins.get_local(self.Plane, res), 'Residuals', **prep_kw(dkw, z_tit='Residuals [Col]', **self.ax_tits(True)))

        def draw_y_map(self, res=.3, cut=None, fid=False, pl=None, **dkw):
            (x, y), z_ = [f(cut=self.Cut.get_nofid(cut, fid), pl=pl) for f in [self.txy, self.dy]]  # noqa
            self.Draw.prof2d(x, y, z_, bins.get_local(self.Plane, res), 'Residuals', **prep_kw(dkw, z_tit='Residuals [Row]', **self.ax_tits(True)))

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

        def draw_xdy(self, cut=None, prof=True, pl=None, **dkw):
            self._draw_angle(self.tx(cut, pl), self.dy(cut, pl), prof, pl=pl, local=True, **prep_kw(dkw, title='X dY', x_tit='X [Cols]', y_tit='dY [Rows]'))

        def draw_ydx(self, cut=None, prof=True, pl=None, **dkw):
            self._draw_angle(self.ty(cut, pl), self.dx(cut, pl), prof, pl=pl, xb=False, local=True, **prep_kw(dkw, title='Y dX', x_tit='Y [Rows]', y_tit='dX [Cols]'))

        def cuts(self):
            return [self.get_pixel(), self.get_cut()]

        def get_pixel(self):
            xy, oxy = self.Plane.PXY * 1e3, uarr2n(self.means()) * 1e3
            return self.Draw.box(*(xy - 1.5 * xy + oxy), *(xy - .5 * xy + oxy), show=False)

        def get_cut(self, local=False):
            (x, y), n = self.means(local) * 1e3, self.Cut.get_config('residuals', dtype=float)
            return self.Draw.ellipse(*self.Plane.PXY / 2 * n * 1e3, x.n, y.n, show=False)
        # endregion DRAW
        # ----------------------------------------

        # ----------------------------------------
        # region ALIGN
        @property
        def m(self):
            return self.align(_save=False)

        @no_trans
        @save_pickle('AM', sub_dir='alignment', run='', suf_args='[6]')
        def align(self, d=None, m=None, cut=None, pl=None, p=.05, i=0, imax=20, _redo=False, _save=False):
            if d is None:
                PBAR.start(imax)
            sx, sy = self.plane(pl).PX, self.plane(pl).PY
            cut = self.Cut(cut) & self.Cut['res']
            x, y = transform(*self.get_xy(local=True, cut=cut, pl=pl), sx, sy) if d is None else d[:2]  # convert to mm
            tx, ty = transform(*self.get_txy(local=True, cut=cut, pl=pl), sx, sy) if d is None else d[2:]
            d = self.rotate(x, y, *((tx, ty) if m is None else m_transform(m, *d[2:])), p=p)
            t = self.translate(*d[:4])
            self.Rot.append(d[4])
            self.Trans.append(t)
            m = choose(m, identity(3)) @ matrix(1, 1, *t, rx=d[4], order='str')
            PBAR.update()
            if i < imax - 1:
                return self.align([x, y, tx, ty], m, p=p, i=i + 1, imax=imax, _redo=_redo, _save=_save)
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

    return ResidualAnalysis
