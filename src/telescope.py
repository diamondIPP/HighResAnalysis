#!/usr/bin/env python
# --------------------------------------------------------
#       Telescope analysis of the high resolution tests at DESY and CERN
# created on June 11th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from analysis import Analysis
from draw import fill_hist, format_histo, array
from ROOT import TH2F
import toml
from numpy import in1d, where


class TelescopeAnalysis(Analysis):

    def __init__(self, dut_analysis):
        self.Ana = dut_analysis
        Analysis.__init__(self, verbose=self.Ana.Verbose)

        self.Data = self.Ana.Data
        self.Bins = self.Ana.Bins
        self.Mask = None
        self.NPlanes = self.Config.getint('TELESCOPE', 'planes')

    def load_mask(self):
        if self.Mask is None:
            d = toml.load(self.Ana.Converter.FileNames[1])['sensors']
            self.Mask = [array(d[i]['masked_pixels']) for i in range(self.NPlanes)]
        return self.Mask

    def get_data(self, plane, grp, key, cut=None):
        data = array(self.Data['Plane{}'.format(plane)][grp][key])
        return data if cut is None else data[cut]

    def get_x(self, plane=0, cluster=False, cut=None):
        return self.get_data(plane, 'Clusters' if cluster else 'Hits', 'X')

    def get_y(self, plane=0, cluster=False, cut=None):
        return self.get_data(plane, 'Clusters' if cluster else 'Hits', 'Y')

    def get_hits(self, plane=0, cut=None):
        return array([self.get_x(plane, cluster=False, cut=cut), self.get_y(plane, cluster=False, cut=cut)])

    def get_clusters(self, plane=0, cut=None):
        return array([self.get_x(plane, cluster=True, cut=cut), self.get_y(plane, cluster=True, cut=cut)])

    def get_mask(self, plane=0):
        x, y = self.get_hits(plane)
        data = x.astype('i') * 10000 + y  # make unique number out of the tuple... Is there a way to compare tuples?
        mx, my = self.load_mask()[plane].T
        mask = mx * 10000 + my
        return where(in1d(data, mask, invert=True))[0]

    def draw_mask(self, plane=0, show=True):
        h = TH2F('htm', 'Masked Pixels in Plane {}'.format(plane), *self.Bins.get_pixel())
        m = self.load_mask()[plane].T
        fill_hist(h, m[0], m[1])
        format_histo(h, x_tit='Column', y_tit='Row', y_off=1.3, fill_color=1)
        self.format_statbox(entries=True)
        self.draw_histo(h, show=show, lm=.12, draw_opt='box')
