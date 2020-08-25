#!/usr/bin/env python
# --------------------------------------------------------
#       Analysis of the reference plane
# created on August 19th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from analysis import Analysis, warning, bins, do_hdf5
from draw import array
from dut import Plane
from numpy import sqrt, zeros


class RefAnalysis(Analysis):

    def __init__(self, dut_analysis):
        self.Ana = dut_analysis
        Analysis.__init__(self, verbose=self.Ana.Verbose)
        self.MetaSubDir = 'REF'

        self.Run = self.Ana.Run
        self.DUT = self.Ana.DUT
        self.Plane = Plane(self.Config.getint('DUT', 'reference plane'), self.Config, 'DUT')
        self.Data = self.Ana.Data[self.Plane.get_name()]
        self.init_cuts()

    # ----------------------------------------
    # region INIT
    def init_cuts(self):
        self.Cuts.register('cluster', self.make_cluster(), 90, 'tracks with a cluster')
        self.Cuts.register('cluster2', self.make_cluster_corr(), 90, 'tracks a cluster in DUT and REF')

    def make_cluster(self, redo=False):
        def f():
            return self.get('Clusters', 'Size') > 0
        return array(do_hdf5(self.make_hdf5_path('clu'), f, redo))

    def make_cluster_corr(self, redo=False):
        def f():
            return (self.Cuts.get('cluster') + self.Ana.Cuts.get('cluster')).Values
        return array(do_hdf5(self.make_hdf5_path('cluCorr'), f, redo))

    def make_residuals(self, redo=False):
        def f():
            cut = zeros(self.Ana.Tracks.N, bool)
            cut[self.Cuts.get('cluster')()] = self.get_residuals() < self.Cuts.get_config('residuals', dtype=float)
            return cut
        return array(do_hdf5(self.make_hdf5_path('res'), f, redo))

    def make_dut_residuals(self):
        return self.make_residuals()[self.Ana.Cuts.get('cluster')()]
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get(self, grp, key):
        if key not in self.Data[grp].keys():
            warning('key {} not found in {}'.format(key, list(self.Data[grp].keys())))
            return
        return array(self.Data[grp][key])

    def get_du(self):
        return self.get('Clusters', 'U') - self.get('Tracks', 'U')[self.Cuts.get('cluster')()]

    def get_dv(self):
        return self.get('Clusters', 'V') - self.get('Tracks', 'V')[self.Cuts.get('cluster')()]

    def get_residuals(self):
        return sqrt(self.get_du() ** 2 + self.get_dv() ** 2)
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw_x_residuals(self):
        self.format_statbox(all_stat=True)
        self.draw_disto(self.get_du(), bins.make(-3, 3, .01), 'X Residuals', x_tit='Residual [mm]', lm=.12, y_off=1.8)

    def draw_y_residuals(self):
        self.format_statbox(all_stat=True)
        self.draw_disto(self.get_dv(), bins.make(-3, 3, .01), 'Y Residuals', x_tit='Residual [mm]', lm=.12, y_off=1.8)

    def draw_residuals(self):
        self.format_statbox(entries=True)
        self.draw_disto(self.get_residuals(), bins.make(0, 6, .01), 'Residuals', x_tit='Residual [mm]')

    def draw_correlation(self, mode='y', show=True):
        c1, c2 = self.Cuts.get('cluster2')()[self.Cuts.get('cluster')()], self.Cuts.get('cluster2')()[self.Ana.Cuts.get('cluster')()]
        v1, v2 = self.get('Clusters', mode.upper())[c1], self.Ana.get_data('Clusters', mode.upper(), cut=False)[c2]
        self.format_statbox(entries=True, x=.84)
        self.draw_histo_2d(v1, v2, bins.get_corr(mode, self.Plane, self.Ana.Plane), x_tit='{} REF'.format(mode), y_tit='{} DUT'.format(mode), show=show)
    # endregion DRAW
    # ----------------------------------------
