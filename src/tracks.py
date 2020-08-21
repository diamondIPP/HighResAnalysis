#!/usr/bin/env python
# --------------------------------------------------------
#       Track analysis of the high resolution tests at DESY and CERN
# created on June 10th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from analysis import Analysis, warning
from draw import fill_hist, format_histo, array
from ROOT import TH2F
from numpy import rad2deg
import bins


class TrackAnalysis(Analysis):

    def __init__(self, dut_analysis):
        self.Ana = dut_analysis
        Analysis.__init__(self, verbose=self.Ana.Verbose)

        self.Data = self.Ana.Data['Tracks']
        self.PlaneData = self.Ana.get_group('Tracks')
        self.N = self.Data['X'].size

    # ----------------------------------------
    # region GET
    def get(self, key, cut=None):
        if key not in self.Data.keys():
            warning('key {} not found in {}'.format(key, list(self.Data.keys())))
            return
        return array(self.Data[key])[self.Cuts(cut)]

    def get_n(self, cut=None):
        """ returns: number of tracks per event. """
        return self.get('NTracks', cut)

    def get_x(self, cut=None, trk_cut=-1, local=True):
        return self.Ana.get_track_data('Tracks', 'X' if local else 'U', cut, trk_cut)

    def get_y(self, cut=None, trk_cut=-1, local=True):
        return self.Ana.get_track_data('Tracks', 'Y' if local else 'V', cut, trk_cut)

    def get_coods(self, local=True, cut=None, trk_cut=-1):
        return (self.get_x(cut, trk_cut), self.get_y(cut, trk_cut)) if local else (self.get_u(cut, trk_cut), self.get_v(cut, trk_cut))

    def get_dof(self, cut=None):
        return self.get('Dof', cut)

    def get_chi2(self, cut=None):
        return self.get('Chi2', cut) / self.get_dof(cut)

    def get_u(self, cut=None, trk_cut=-1):
        return self.Ana.get_track_data('Tracks', 'U', cut, trk_cut)

    def get_v(self, cut=None, trk_cut=-1):
        return self.Ana.get_track_data('Tracks', 'V', cut, trk_cut)

    def get_events(self, cut=False):
        return self.get('EvtFrame', cut=cut)
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw_n(self, show=True):
        self.format_statbox(all_stat=True)
        return self.draw_disto(self.get_n(), 'Number of Tracks', bins.make(0, 10, 1), show=show, x_tit='Number of Tracks', y_off=2.1, lm=.14)

    def draw_dof(self, show=True, cut=None):
        self.format_statbox(all_stat=True)
        return self.draw_disto(self.get_dof(cut), bins.make(0, 20, 1), 'Track Degrees of Freedom', show=show, x_tit='Degrees of Freedom', y_off=2.1, lm=.14)

    def draw_occupancy0(self, scale=4, cut=None):
        h = TH2F('hto', 'Track Occupancy', *bins.get_global(self.Ana.Telescope.Plane, scale))
        fill_hist(h, *self.get_coods(cut))
        format_histo(h, x_tit='Track X [mm]', y_tit='Track Y [mm]', y_off=1.3, z_tit='Number of Entries', z_off=1.2)
        self.format_statbox(all_stat=True, x=.83)
        self.draw_histo(h, lm=.12, draw_opt='colz', rm=.15)

    def draw_occupancy(self, scale=4, cut=None, raw=False, show=True):
        self.format_statbox(all_stat=True, x=.83)
        x, y = self.get_coods(cut) if raw else (self.get_u(trk_cut=True), self.get_v(cut, trk_cut=True))
        self.draw_histo_2d(x, y, bins.get_global(self.Ana.Telescope.Plane, scale), 'Track Occupancy', x_tit='Track X [mm]', y_tit='Track Y [mm]', show=show)

    def draw_map(self, bin_width=.1, cut=None, dut_plane=True, show=True):
        self.format_statbox(all_stat=True, x=.83)
        binning = bins.get_global(self.Ana.Plane if dut_plane else self.Ana.Telescope.Plane, bin_width)
        self.draw_histo_2d(self.get_u(cut), self.get_v(cut), binning, 'Track Map', x_tit='Track X [mm]', y_tit='Track Y [mm]', show=show)

    def draw_chi2(self, cut=None):
        self.format_statbox(all_stat=True)
        values = self.get_chi2(cut)
        self.draw_disto(values, bins.make(0, 100, 1.), 'Track #chi^{2}', x_tit='#chi^{2}')

    def draw_slope(self, mode='x', bin_width=.01):
        values = rad2deg(array(self.Data['Slope{}'.format(mode.title())]))
        self.format_statbox(all_stat=True)
        self.draw_disto(values, bins.make(-1, 1, bin_width), 'Track Slope {}'.format(mode.title()), x_tit='Track Slope [deg]', y_off=1.8, lm=.13)

    def draw_slope_x(self, bin_width=.01):
        self.draw_slope('x', bin_width)

    def draw_slope_y(self, bin_width=.01):
        self.draw_slope('y', bin_width)
    # endregion DRAW
    # ----------------------------------------
