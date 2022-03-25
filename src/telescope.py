#!/usr/bin/env python
# --------------------------------------------------------
#       Telescope analysis of the high resolution tests at DESY and CERN
# created on June 11th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

import src.bins as bins
from src.analysis import *
from src.cut import Cut
from src.dut import Plane
from src.utils import save_hdf5, get_trans_matrix, dot


class TelescopeAnalysis(Analysis):

    def __init__(self, dut_analysis):
        self.Ana = dut_analysis
        Analysis.__init__(self, verbose=self.Ana.Verbose)
        self.MetaSubDir = 'tel'

        self.Run = self.Ana.Run
        self.DUT = self.Ana.DUT
        self.Data = self.Ana.Data
        self.Mask = None
        self.Plane = Plane(0, self.Ana.Config)
        self.NPlanes = self.Config.getint('TELESCOPE', 'planes')

        self.init_cuts()

    # ----------------------------------------
    # region CUT
    def init_cuts(self):
        self.Cuts.register('1cluster', self.make_1cluster_cut(pl=2), 90, 'events with 1 cluster in plane 2')

    @save_hdf5('clu1', arr=True, suf_args='all')
    def make_1cluster_cut(self, pl, _redo=False):
        return self.get_n_clusters(pl, cut=False) == 1

    @save_hdf5('all1', arr=True, suf_args='all')
    def make_all_cluster(self, n=1, _redo=False):
        return all(array([self.get_n_clusters(pl, cut=False) for pl in arange(self.NPlanes)]) == n, axis=0)

    @save_hdf5('all-i', arr=True, suf_args='all')
    def make_all_intercepts(self, n=1, _redo=False):
        return all(array([self.get_n_intercepts(pl, cut=False) for pl in arange(self.NPlanes)]) == n, axis=0)

    def make_corr_cuts(self, p0, p1, n=1):
        n0, n1 = [self.get_n_clusters(p, cut=False) for p in [p0, p1]]
        n = (n0 == n) & (n1 == n)
        return n.repeat(n0), n.repeat(n1)

    def make_correlation(self, pl):
        npl = self.get_n_clusters(pl, cut=False)
        dut_events = zeros(npl.size, dtype='?')
        dut_events[self.Ana.get_events()] = True
        return Cut('correlation', Cut.add(self.make_1cluster_cut(pl), dut_events.repeat(npl)), 90, '1 cluster per event in plane {} and cluster and dut'.format(pl))

    def ev_cut(self, ev):
        c = zeros(self.Ana.NEvents, '?')
        c[make_list(ev).astype('i')] = True
        return c

    def ev2pl(self, x, pl=0):
        n = self.get_n_clusters(pl, cut=False)
        return x.repeat(n) if is_iter(x) else ...
    # endregion CUT
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get(self, plane=0):
        return self.Data['Plane{}'.format(self.Plane(plane).Number)]

    def get_group(self, grp, plane=0):
        return self.get(plane)[grp]

    def get_data(self, grp, key=None, pl=0, cut=None):
        data = self.get_group(grp, pl)
        data = array(data) if key is None else array(data[key])
        return data[self.ev2pl(self.Cuts(cut), pl)] if is_iter(cut) else data

    def get_n(self, name='Hits', plane=0, cut=None):
        return self.get_data(name, f'N{name}', plane, cut)

    def get_n_clusters(self, pl=0, cut=None):
        return self.get_n('Clusters', pl, cut)

    def get_n_intercepts(self, pl=0, cut=None):
        return self.get_n('Intercepts', pl, cut)

    def get_n_tracks(self, pl=0, cut=None):
        return self.get_data('Clusters', 'NTracks', pl, cut)

    def get_all(self, f, cut=None, **kwargs):
        return [f(pl, cut=cut, **kwargs) for pl in range(self.NPlanes)]

    def get_ix(self, pl=0):
        return self.l2g(self.get_data('Intercepts', 'X', pl, cut=False), pl=pl, scale=False)[0]

    def get_ixs(self, cut=None):
        return self.get_all(self.get_ix, cut)

    def get_iy(self, pl=0, cut=None):
        return self.l2g(y=self.get_data('Intercepts', 'Y', pl, cut), pl=pl, scale=False)[1]

    def get_iys(self, cut=None):
        return self.get_all(self.get_iy, cut)

    def get_x(self, pl=0, cluster=True, cut=None):
        return self.get_data('Clusters' if cluster else 'Hits', 'X', pl, cut)

    def get_xs(self, cluster=True, cut=None):
        return [self.get_x(pl, cluster, cut) for pl in range(self.NPlanes)]

    def get_u(self, pl=0, cut=None):
        return self.l2g(self.get_x(pl, cut=cut), pl=pl)[0]

    def get_us(self, cut=None):
        return [self.get_u(pl, cut) for pl in range(self.NPlanes)]

    def get_y(self, plane=0, cluster=True, cut=None):
        return self.get_data('Clusters' if cluster else 'Hits', 'Y', plane, cut)

    def get_v(self, pl=0, cut=None):
        return self.l2g(y=self.get_y(pl, cut=cut), pl=pl)[1]

    def get_xy(self, pl=0, cluster=True, cut=None, local=True):
        x, y = self.get_x(pl, cluster, cut), self.get_y(pl, cluster, cut)
        return (x, y) if local else self.l2g(x, y, pl)

    def ax_tits(self, local=True):
        return self.Ana.get_ax_tits(local)

    @property
    def scale_matrix(self):
        m = self.Plane.M
        return array([append(m[0], 0), append(m[1], 0), [0, 0, 1]])

    def rot_matrix(self, pl=0):
        d = self.Ana.Converter.get_alignment()['sensors'][pl]
        return array([append(d['unit_u'][:2], 0), append(d['unit_v'][:2], 0), [0, 0, 1]])

    def get_offsets(self, pl=0):
        return self.Ana.Converter.get_alignment()['sensors'][pl]['offset'][:2]

    def trans_matrix(self, pl=0, centre=True):
        o = self.get_offsets(pl)
        return get_trans_matrix(o - array([self.Plane.W / 2, self.Plane.H / 2]) if centre else o)

    def ms(self, pl):
        return dot(dot(self.trans_matrix(pl), self.rot_matrix(pl)), self.scale_matrix)

    def m(self, pl):
        return dot(self.trans_matrix(pl, centre=False), self.rot_matrix(pl))

    def l2g(self, x=None, y=None, pl=0, scale=True):
        return dot(self.ms(pl) if scale else self.m(pl), array([ones(y.size, 'i2') if x is None else x, ones(x.size, 'i2') if y is None else y, ones(choose(x, y).size)]))[:2]
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw_n(self, plane=0, name='Hits', show=True):
        n, n_pl = name, self.Plane(plane).get_name()
        self.Draw.distribution(self.get_n(n, plane), bins.make(0, 30), 'Number of {} in {}'.format(n, n_pl), lm=.13, show=show, x_tit='Number of {}'.format(n), y_off=2)

    def draw_occupancy(self, pl=0, local=True, cluster=True, bw=10, cut=None, **dkw):
        x, y = self.get_xy(pl, cluster, cut, local)
        title = f'{"Cluster" if cluster else "Hit"} Occupancy of Plane {pl}'
        self.Draw.histo_2d(x, y, bins.get_xy(local, self.Plane, bw, aspect_ratio=True), **prep_kw(dkw, title=title, **self.ax_tits(local), stats=set_statbox(entries=True, m=True)))

    def draw_n_hits(self, plane=0, show=True):
        self.draw_n(plane, 'Hits', show)

    def draw_n_clusters(self, plane=0, show=True):
        self.draw_n(plane, 'Clusters', show)

    def draw_n_intercepts(self, plane=0, show=True):
        self.draw_n(plane, 'Intercepts', show)

    def draw_x_correlation(self, p0=0, p1=1, n=1, **dkw):
        c = self.make_corr_cuts(p0, p1, n)
        x, y = [self.get_x(p, cut=c[i]) for i, p in enumerate([p0, p1])]
        return self.Draw.histo_2d(x, y, **prep_kw(dkw, title='XCorr', x_tit='Column', y_tit='Row'))

    def draw_y_correlation(self, p0=0, p1=1, n=1, **dkw):
        c = self.make_corr_cuts(p0, p1, n)
        x, y = [self.get_y(p, cut=c[i]) for i, p in enumerate([p0, p1])]
        return self.Draw.histo_2d(x, y, **prep_kw(dkw, title='YCorr', x_tit='Column', y_tit='Row'))
    # endregion DRAW
    # ----------------------------------------
