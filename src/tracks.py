#!/usr/bin/env python
# --------------------------------------------------------
#       Track analysis of the high resolution tests at DESY and CERN
# created on June 10th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from analysis import Analysis
from draw import fill_hist, format_histo, array
from ROOT import TH2F
from numpy import rad2deg
import bins


class TrackAnalysis(Analysis):

    def __init__(self, dut_analysis):
        self.Ana = dut_analysis
        Analysis.__init__(self, verbose=self.Ana.Verbose)

        self.Data = self.Ana.Data

    def draw_n(self, show=True):
        self.format_statbox(all_stat=True)
        return self.draw_disto(self.Data['NTracks'], 'Number of Tracks', bins.make(0, 10, 1), show=show, x_tit='Number of Tracks', y_off=2.1, lm=.14)

    def draw_dof(self, show=True):
        self.format_statbox(all_stat=True)
        return self.draw_disto(self.Data['Dof'], 'Track Degrees of Freedom', bins.make(0, 20, 1), show=show, x_tit='Degrees of Freedom', y_off=2.1, lm=.14)

    def draw_occupancy(self, scale=4):
        h = TH2F('hto', 'Track Occupancy', *bins.get_global(self.Ana.Telescope.Plane, scale))
        fill_hist(h, self.Data['X'], self.Data['Y'])
        format_histo(h, x_tit='Track X [mm]', y_tit='Track Y [mm]', y_off=1.3, z_tit='Number of Entries', z_off=1.2)
        self.format_statbox(all_stat=True, x=.83)
        self.draw_histo(h, lm=.12, draw_opt='colz', rm=.15)

    def draw_chi2(self):
        values = array(self.Data['Chi2']) / array(self.Data['Dof'])
        self.draw_disto(values, 'Track #chi^{2}', bins.make(0, 100, 1.), x_tit='#chi^{2}')

    def draw_slope(self, mode='x', bin_width=.01):
        values = rad2deg(array(self.Data['Slope{}'.format(mode.title())]))
        self.format_statbox(all_stat=True)
        self.draw_disto(values, 'Track Slope {}'.format(mode.title()), bins.make(-1, 1, bin_width), x_tit='Track Slope [deg]', y_off=1.8)

    def draw_slope_x(self, bin_width=.01):
        self.draw_slope('x', bin_width)

    def draw_slope_y(self, bin_width=.01):
        self.draw_slope('y', bin_width)
