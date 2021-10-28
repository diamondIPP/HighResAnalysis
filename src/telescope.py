#!/usr/bin/env python
# --------------------------------------------------------
#       Telescope analysis of the high resolution tests at DESY and CERN
# created on June 11th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from src.analysis import *
from src.utils import do_hdf5
from src.dut import Plane
from src.cut import Cut
import src.bins as bins


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
        self.Cuts.register('1cluster', self.make_1cluster_cut(plane=2), 90, 'events with 1 cluster in plane 2')

    def make_1cluster_cut(self, plane, redo=False):
        def f():
            n = self.get_n('Clusters', plane=plane, cut=False)
            return n.repeat(n) == 1
        return array(do_hdf5(self.make_hdf5_path('clu1', plane), f, redo))

    def make_correlation(self, plane):
        npl = self.get_n('Clusters', plane=plane, cut=False)
        dut_events = zeros(npl.size, dtype='?')
        dut_events[self.Ana.get_events()] = True
        return Cut('correlation', Cut.add(self.make_1cluster_cut(plane), dut_events.repeat(npl)), 90, '1 cluster per event in plane {} and cluster and dut'.format(plane))
    # endregion CUT
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get(self, plane=0):
        return self.Data['Plane{}'.format(self.Plane(plane).Number)]

    def get_group(self, grp, plane=0):
        return self.get(plane)[grp]

    def get_data(self, grp, key=None, plane=0, cut=None):
        data = self.get_group(grp, plane)
        data = array(data) if key is None else array(data[key])
        return data[self.Cuts(cut)]

    def get_n(self, name='Hits', plane=0, cut=None):
        return self.get_data(name, 'N{}'.format(name), plane, cut)

    def get_x(self, plane=0, cluster=True, cut=None):
        return self.get_data('Clusters' if cluster else 'Hits', 'X', plane, cut)

    def get_y(self, plane=0, cluster=True, cut=None):
        return self.get_data('Clusters' if cluster else 'Hits', 'Y', plane, cut)

    def get_coods(self, plane=0, cluster=True, cut=None):
        return self.get_x(plane, cluster, cut), self.get_y(plane, cluster, cut)
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw_n(self, plane=0, name='Hits', show=True):
        n, n_pl = name, self.Plane(plane).get_name()
        self.Draw.distribution(self.get_n(n, plane), bins.make(0, 30), 'Number of {} in {}'.format(n, n_pl), lm=.13, show=show, x_tit='Number of {}'.format(n), y_off=2)

    def draw_occupancy(self, plane=0, cluster=True, bin_width=10, cut=None, show=True):
        x, y = self.get_coods(plane, cluster, cut)
        title = '{} Occupancy of Plane {}'.format('Cluster' if cluster else 'Hit', plane)
        self.Draw.histo_2d(x, y, title, bins.get_local(self.Plane, bin_width, aspect_ratio=True), x_tit='Column', y_tit='Row', show=show, stats=set_statbox(entries=True, m=True))

    def draw_n_hits(self, plane=0, show=True):
        self.draw_n(plane, 'Hits', show)

    def draw_n_clusters(self, plane=0, show=True):
        self.draw_n(plane, 'Clusters', show)

    def draw_n_intercepts(self, plane=0, show=True):
        self.draw_n(plane, 'Intercepts', show)

    def draw_correlation(self, mode='y', res=1, plane=2, thresh=.1, show=True):
        v0 = getattr(self, 'get_{}'.format(mode.lower()))(plane, cut=self.make_correlation(plane))
        v1 = getattr(self.Ana, 'get_{}'.format(mode.lower()))(cut=self.Ana.make_correlation(plane))
        tit = '{} Correlation between {} and {}'.format(mode.upper(), self.Plane, self.Ana.Plane)
        h = self.Draw.histo_2d(v0, v1, bins.get_corr(mode, self.Plane, self.Ana.Plane, res), tit, x_tit=mode.upper(), y_tit=mode.upper(), show=show)
        format_histo(h, z_range=[thresh * h.GetMaximum(), h.GetMaximum()])
        update_canvas()
        return h
    # endregion DRAW
    # ----------------------------------------
