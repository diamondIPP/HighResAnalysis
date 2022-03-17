#!/usr/bin/env python
# --------------------------------------------------------
#       Analysis of the reference plane
# created on August 19th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from src.analysis import Analysis, warning
from src.utils import do_hdf5, ufloat
from plotting.draw import array, format_histo, update_canvas
from src.dut import Plane
from numpy import sqrt, zeros
import src.bins as bins


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
        self.Cuts.set_config(self.BeamTest, self.DUT.Name)
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

    def make_dut_residuals(self, redo=False):
        return self.make_residuals(redo)[self.Ana.Cuts.get('cluster')()]
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get(self, grp, key):
        if key not in self.Data[grp].keys():
            warning('key {} not found in {}'.format(key, list(self.Data[grp].keys())))
            return
        return array(self.Data[grp][key])

    def get_du(self, centred=False):
        u = self.get('Clusters', 'U') - self.get('Tracks', 'U')[self.Cuts.get('cluster')()]
        return u - self.get_md('x').n if centred else u

    def get_dv(self, centred=False):
        v = self.get('Clusters', 'V') - self.get('Tracks', 'V')[self.Cuts.get('cluster')()]
        return v - self.get_md('y').n if centred else v

    def get_md(self, mode):
        h = getattr(self, 'draw_{}_residuals'.format(mode))(show=False)
        fit = h.Fit('gaus', 'sq0', '', -1, 1)
        return ufloat(fit.Parameter(1), fit.ParError(1))

    def get_residuals(self, centred=False):
        return sqrt(self.get_du(centred) ** 2 + self.get_dv(centred) ** 2)
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw_x_residuals(self, show=True, centred=False):
        return self.Draw.distribution(self.get_du(centred), bins.make(-2, 2, .01), 'X Residuals', x_tit='Residual [mm]', lm=.12, y_off=1.8, show=show)

    def draw_y_residuals(self, show=True, centred=False):
        return self.Draw.distribution(self.get_dv(centred), bins.make(-2, 2, .01), 'Y Residuals', x_tit='Residual [mm]', lm=.12, y_off=1.8, show=show)

    def draw_residuals_map(self, centred=False):
        x, y = self.get_du(centred), self.get_dv(centred)
        self.Draw.histo_2d(x, y, bins.make(-1, 1, .01) * 2)

    def draw_residuals(self, centred=False):
        self.Draw.distribution(self.get_residuals(centred), bins.make(0, 6, .01), 'Residuals', x_tit='Residual [mm]', lm=.13, y_off=2)

    def draw_correlation(self, mode='y', thresh=.01, show=True):
        c1, c2 = self.Cuts.get('cluster2')()[self.Cuts.get('cluster')()], self.Cuts.get('cluster2')()[self.Ana.Cuts.get('cluster')()]
        v1, v2 = self.get('Clusters', mode.upper())[c1], self.Ana.get_data('Clusters', mode.upper(), cut=False)[c2]
        h = self.Draw.histo_2d(v1, v2, bins.get_corr(mode, self.Plane, self.Ana.Plane), x_tit='{} REF'.format(mode), y_tit='{} DUT'.format(mode), show=show)
        format_histo(h, z_range=[thresh * h.GetMaximum(), h.GetMaximum()])
        update_canvas()
    # endregion DRAW
    # ----------------------------------------
