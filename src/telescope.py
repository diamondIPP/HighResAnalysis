#!/usr/bin/env python
# --------------------------------------------------------
#       Telescope analysis of the high resolution tests at DESY and CERN
# created on June 11th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from analysis import *
from dut import Plane


class TelescopeAnalysis(Analysis):

    def __init__(self, dut_analysis):
        self.Ana = dut_analysis
        Analysis.__init__(self, verbose=self.Ana.Verbose)

        self.Data = self.Ana.Data
        self.Mask = None
        self.Plane = Plane(0, self.Ana.Config)
        self.NPlanes = self.Config.getint('TELESCOPE', 'planes')

        self.init_cuts()
        
    def init_cuts(self):
        self.Cuts.register('1cluster', self.get_n('Clusters', plane=2, cut=False) == 1, 90, 'events with 1 cluster in plane 2')

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

    def get_1_cluster_cut(self, plane=0):
        return self.get_n('Clusters', plane, cut=False) == 1
    # endregion GET
    # ----------------------------------------

    def draw_n(self, plane=0, name='Hits', show=True):
        self.format_statbox(all_stat=True)
        n, n_pl = name, self.Plane(plane).get_name()
        self.draw_disto(self.get_n(n, plane), 'Number of {} in {}'.format(n, n_pl), bins.make(0, 30), lm=.13, show=show, x_tit='Number of {}'.format(n), y_off=2)

    def draw_occupancy(self, plane=0, cluster=True, bin_width=10, cut=None, show=True):
        self.format_statbox(entries=True, x=.83, m=True)
        x, y = self.get_coods(plane, cluster, cut)
        title = '{} Occupancy of Plane {}'.format('Cluster' if cluster else 'Hit', plane)
        self.draw_histo_2d(x, y, title, bins.get_local(self.Plane, bin_width, aspect_ratio=True), x_tit='Column', y_tit='Row', show=show)

    def draw_n_hits(self, plane=0, show=True):
        self.draw_n(plane, 'Hits', show)

    def draw_n_clusters(self, plane=0, show=True):
        self.draw_n(plane, 'Clusters', show)

    def draw_n_intercepts(self, plane=0, show=True):
        self.draw_n(plane, 'Intercepts', show)
