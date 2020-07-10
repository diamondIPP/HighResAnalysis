#!/usr/bin/env python
# --------------------------------------------------------
#       Track analysis of the high resolution tests at DESY and CERN
# created on June 10th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from analysis import Analysis, warning
from draw import fill_hist, format_histo, array
from ROOT import TH2F
from numpy import rad2deg, arange, sqrt
import bins


class TrackAnalysis(Analysis):

    def __init__(self, dut_analysis):
        self.Ana = dut_analysis
        Analysis.__init__(self, verbose=self.Ana.Verbose)

        self.Data = self.Ana.Data['Tracks']
        self.N = self.Data['X'].size
        self.init_cuts()

    # ----------------------------------------
    # region INIT
    def init_cuts(self):
        dof = self.get_dof(cut=False)
        self.Cuts.register('dof', (dof > 10) & (dof < 20), 10, 'degress of freedom > 10')
        matched = self.get('Matched', cut=False)
        self.Cuts.register('matched', matched, 10, '{} matched tracks'.format(matched.nonzero()[0].size))
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get(self, key, cut=None):
        if key not in self.Data.keys():
            warning('key {} not found in {}'.format(key, list(self.Data.keys())))
            return
        return self.get_values(array(self.Data[key]), self.Cuts(cut))

    def get_n(self):
        return self.get('NTracks')

    def get_x(self, cut=None):
        return self.get('X', cut)

    def get_y(self, cut=None):
        return self.get('Y', cut)

    def get_coods(self, cut):
        return self.get_x(cut), self.get_y(cut)

    def get_dof(self, cut=None):
        return self.get('Dof', cut)

    def get_chi2(self, cut=None):
        return self.get('Chi2', cut)
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw_n(self, show=True):
        self.format_statbox(all_stat=True)
        return self.draw_disto(self.get_n(), 'Number of Tracks', bins.make(0, 10, 1), show=show, x_tit='Number of Tracks', y_off=2.1, lm=.14)

    def draw_dof(self, show=True, cut=None):
        self.format_statbox(all_stat=True)
        return self.draw_disto(self.get_dof(cut), 'Track Degrees of Freedom', bins.make(0, 20, 1), show=show, x_tit='Degrees of Freedom', y_off=2.1, lm=.14)

    def draw_dof_trend(self, cut=False, show=True):
        self.format_statbox(entries=True, x=.9)
        self.draw_prof(arange(self.N), self.get_dof(cut), bins.make(0, self.N, sqrt(self.N)), 'DOF Trend', x_tit='Track Number', y_tit='Degrees of Freedom', show=show, rm=.08)

    def draw_occupancy(self, scale=4, cut=None):
        h = TH2F('hto', 'Track Occupancy', *bins.get_global(self.Ana.Telescope.Plane, scale))
        fill_hist(h, *self.get_coods(cut))
        format_histo(h, x_tit='Track X [mm]', y_tit='Track Y [mm]', y_off=1.3, z_tit='Number of Entries', z_off=1.2)
        self.format_statbox(all_stat=True, x=.83)
        self.draw_histo(h, lm=.12, draw_opt='colz', rm=.15)

    def draw_chi2(self, cut=None):
        self.format_statbox(all_stat=True)
        values = self.get_chi2(cut) / self.get_dof(cut)
        self.draw_disto(values, 'Track #chi^{2}', bins.make(0, 100, 1.), x_tit='#chi^{2}')

    def draw_slope(self, mode='x', bin_width=.01):
        values = rad2deg(array(self.Data['Slope{}'.format(mode.title())]))
        self.format_statbox(all_stat=True)
        self.draw_disto(values, 'Track Slope {}'.format(mode.title()), bins.make(-1, 1, bin_width), x_tit='Track Slope [deg]', y_off=1.8)

    def draw_slope_x(self, bin_width=.01):
        self.draw_slope('x', bin_width)

    def draw_slope_y(self, bin_width=.01):
        self.draw_slope('y', bin_width)
    # endregion DRAW
    # ----------------------------------------
