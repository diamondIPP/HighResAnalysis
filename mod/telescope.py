#!/usr/bin/env python
# --------------------------------------------------------
#       Telescope analysis of the high resolution tests at DESY and CERN
# created on June 11th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from src.dut_analysis import DUTAnalysis
import src.bins as bins
from src.analysis import *
from src.utils import get_trans_matrix, dot
from numpy.linalg import inv
from numpy import identity
from mod.tel_cuts import TelCut


class TelescopeAnalysis(DUTAnalysis):

    def __init__(self, parent: DUTAnalysis):  # noqa

        self.Parent = parent
        self.__dict__.update(parent.__dict__)
        self.MetaSubDir = 'tel'
        
        self.NPlanes = self.Config.getint('TELESCOPE', 'planes')
        self.Plane = self.Planes[0]
        self.Cut = TelCut(self)

    # ----------------------------------------
    # region MISC
    @property
    def planes(self):
        return arange(self.NPlanes)
    # endregion MISC
    # -----------------------------------------

    # ----------------------------------------
    # region DATA
    def get_time(self, cut=None):
        return self.Cut(cut, super().get_time(...))

    def get_data(self, grp, key=None, cut=None, pl=0):
        self.Data = self.F[str(self.Planes[pl])]
        return super().get_data(grp, key, cut, pl)

    def get_all(self, f, cut=True, **kwargs):
        return [f(pl=pl, cut=self.Cut.make_all_cluster() & cut, **kwargs) for pl in range(self.NPlanes)]

    def get_n_clusters(self, cut=None, pl=0):
        return self.get_data('Clusters', 'N', cut, pl)

    def get_cluster_size(self, cut=None, pl=0):
        return self.get_data('Clusters', 'Size', pl, cut)

    def get_u(self, cut=None, pl=0):
        return self.l2g(self.get_x(cut, pl), pl=pl, centre=True)[0]

    def get_us(self, cut=True):
        return array(self.get_all(self.get_u, cut))

    def get_v(self, pl=0, cut=None):
        return self.l2g(y=self.get_y(cut, pl), pl=pl, centre=True)[1]
    # endregion DATA
    # ----------------------------------------

    # ----------------------------------------
    # region COORDINATE TRANSFORM
    def scale_matrix(self, pl):
        m = self.Plane.M if pl < self.NPlanes else self.Parent.Plane.M
        return array([append(m[0], 0), append(m[1], 0), [0, 0, 1]])

    def s(self, x=None, y=None, pl=0):
        return dot(self.scale_matrix(pl), array([ones(y.size, 'i2') if x is None else x, ones(x.size, 'i2') if y is None else y, ones(choose(x, y).size)]))[:2]

    def rot_matrix(self, pl=0):
        d = self.Converter.get_alignment()['sensors'][pl]
        return array([append(d['unit_u'][:2], 0), append(d['unit_v'][:2], 0), [0, 0, 1]])

    def get_offsets(self, pl=0):
        return self.Converter.get_alignment()['sensors'][pl]['offset'][:2]

    def trans_matrix(self, pl=0, centre=False):
        o = array(self.get_offsets(pl))
        return get_trans_matrix(o - array([self.Plane.W / 2, self.Plane.H / 2]) if centre else o)

    def m(self, pl, scale=True, centre=False):
        s, r, t = self.scale_matrix(pl) if scale else identity(3), self.rot_matrix(pl), self.trans_matrix(pl, centre)
        return dot(t, dot(r, s))

    def l2g(self, x=None, y=None, pl=0, scale=True, centre=False):
        return dot(self.m(pl, scale, centre), array([ones(y.size, 'i2') if x is None else x, ones(x.size, 'i2') if y is None else y, ones(choose(x, y).size)]))[:2]

    def g2l(self, x=None, y=None, pl=0, scale=True):
        return dot(inv(self.m(pl, scale)), array([ones(y.size, 'i2') if x is None else x, ones(x.size, 'i2') if y is None else y, ones(choose(x, y).size)]))[:2]
    # endregion COORDINATE TRANSFORM
    # ----------------------------------------

    # ----------------------------------------
    # region RESIDUAL
    def get_tu(self, pl=0, cut=None):
        return self.l2g(self.get_data('Tracks', 'X', pl, cut), pl=pl, centre=True)[0]

    def get_du_(self, pl=0, cut=None):
        return self.get_data('Tracks', 'dU', pl, cut) * 1e3

    def get_du(self, pl=0, cut=None):
        return self.Cut(cut, pl, self.get_u(pl) - self.get_tu(pl, self.Cut.make_cluster(pl))) * 1e3

    def get_dv(self, pl=0, cut=None):
        return self.get_data('Tracks', 'dV', pl, cut) * 1e3

    def draw_x_residual(self, pl=0, cut=None, bw=10, **dkw):
        self.Draw.distribution(self.get_du(pl, cut), bins.make(-1000, 1000, bw), **prep_kw(dkw, title='XRes', x_tit='dX [#mum]'))

    def draw_y_residual(self, pl=0, cut=None, **dkw):
        self.Draw.distribution(self.get_dv(pl, cut), **prep_kw(dkw, title='XRes', x_tit='dX [#mum]'))
    # endregion RESIDUAL
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw_occupancy(self, pl=0, local=True, bw=10, cut=None, **dkw):
        return super().draw_occupancy(local, bw, cut, pl=pl, **dkw)

    def draw_n_clusters(self, pl=0, cut=None, **dkw):
        return self.Draw.distribution(self.get_n_clusters(pl, cut), **prep_kw(dkw, title='NClusters', x_tit='Number of Clusters', w=1, x0=-.5))

    def draw_x_correlation(self, pl0=0, pl1=1, **dkw):
        c = self.Cut.make_correlation(pl0, pl1)
        x, y = self.get_x(c, pl0), self.get_x(c, pl1)
        return self.Draw.histo_2d(x, y, **prep_kw(dkw, title='XCorr', x_tit=f'Column Plane {pl0}', y_tit=f'Column Plane {pl1}'))

    def draw_y_correlation(self, pl0=0, pl1=1, **dkw):
        c = self.Cut.make_correlation(pl0, pl1)
        x, y = self.get_y(c, pl0), self.get_y(c, pl1)
        return self.Draw.histo_2d(x, y, **prep_kw(dkw, title='YCorr', x_tit=f'Row Plane {pl0}', y_tit=f'Row Plane {pl1}'))

    def draw_correlation_trend(self, pl0=0, pl1=1, thresh=.2, **dkw):
        c = self.Cut.make_correlation(pl0, pl1)
        d0, d1, t = self.get_xy(pl=pl0, cut=c), self.get_xy(pl=pl1, cut=c), self.get_time(c)
        g = [self.Draw.graph(*get_3d_correlations(self.Draw.histo_3d(t, d0[i], d1[i]), thresh=thresh), y_tit='Correlation Factor', show=False) for i in range(2)]
        return self.Draw.multigraph(g, 'CorrFac', ['x', 'y'], draw_opt='pl', **prep_kw(dkw, **self.t_args(), y_range=[0, 1.05]))
    # endregion DRAW
    # ----------------------------------------
