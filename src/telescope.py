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

    def draw_occupancy(self, plane=0, cluster=True, bin_width=10, show=True):
        return self.Ana.draw_occupancy(self.Plane(plane), cluster, bin_width, show=show)

    def draw_n_hits(self, plane=0, show=True):
        self.Ana.draw_n(self.Plane(plane), 'Hits', show)

    def draw_n_clusters(self, plane=0, show=True):
        self.Ana.draw_n(self.Plane(plane), 'Clusters', show)

    def draw_n_intercepts(self, plane=0, show=True):
        self.Ana.draw_n(self.Plane(plane), 'Intercepts', show)

    def draw_cluster_size(self, plane=0, show=True):
        self.Ana.draw_cluster_size(self.Plane(plane), show)

    def get_1_cluster_cut(self, plane=0):
        n = self.Ana.get_n_clusters(self.Plane(plane))
        return n == 1
