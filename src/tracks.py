#!/usr/bin/env python
# --------------------------------------------------------
#       Track analysis of the high resolution tests at DESY and CERN
# created on June 10th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from analysis import Analysis
from draw import fill_hist, format_histo, array
from ROOT import TH2F
from binning import make_bins
from numpy import rad2deg


class TrackAnalysis(Analysis):

    def __init__(self, dut_analysis):
        self.Ana = dut_analysis
        Analysis.__init__(self, verbose=self.Ana.Verbose)

        self.Data = self.Ana.Data
        self.Bins = self.Ana.Bins

    def draw_occupancy(self, scale=4):
        h = TH2F('hto', 'Track Occupancy', *self.Bins.get_global(scale))
        fill_hist(h, self.Data['X'], self.Data['Y'])
        format_histo(h, x_tit='Track X [mm]', y_tit='Track Y [mm]', y_off=1.3, z_tit='Number of Entries', z_off=1.2)
        self.format_statbox(all_stat=True, x=.83)
        self.draw_histo(h, lm=.12, draw_opt='colz', rm=.15)

    def draw_chi2(self):
        values = array(self.Data['Chi2']) / array(self.Data['Dof'])
        self.draw_disto(values, 'Track #chi^{2}', make_bins(0, 100, 1.), x_tit='#chi^{2}')

    def draw_slope(self, mode='x', bin_width=.01):
        values = rad2deg(array(self.Data['Slope{}'.format(mode.title())]))
        self.format_statbox(all_stat=True)
        self.draw_disto(values, 'Track Slope {}'.format(mode.title()), make_bins(-1, 1, bin_width), x_tit='Track Slope [deg]', y_off=1.8)

    def draw_slope_x(self, bin_width=.01):
        self.draw_slope('x', bin_width)

    def draw_slope_y(self, bin_width=.01):
        self.draw_slope('y', bin_width)