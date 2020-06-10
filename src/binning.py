#!/usr/bin/env python
# --------------------------------------------------------
#       class for creating the histogram bins for the analysis
# created on June 10th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from numpy import arange
from dut import Plane


class Bins:
    def __init__(self, config):

        self.Plane = Plane(config)

        self.NCols, self.NRows = self.Plane.NCols, self.Plane.NRows
        self.PX, self.PY = self.Plane.PX, self.Plane.PY

        # Pixel
        self.MaxADC = 2**8 - 1
        self.MinPH = -5000
        self.MaxPH = 100000
        self.PHBinWidth = 200
        self.MinVcal = -100
        self.MaxVcal = 1250
        self.VcalToEl = 47.

        # Miscellaneous
        self.GlobalCoods = [-.5025, .5175, -.505, .515]  # range in x and y in telescope coordinates [cm]

    def get_pixel(self):
        return self.get_pixel_x() + self.get_pixel_y()

    def get_pixel_x(self):
        return [self.NCols, arange(-.5, self.NCols)]

    def get_pixel_y(self):
        return [self.NRows, arange(-.5, self.NRows)]

    def get_global(self, scale=1):
        return self.get_global_x(scale) + self.get_global_y(scale)

    def get_global_cood(self, mode, scale=1):
        return self.get_global_x(scale) if mode == 'x' else self.get_global_y(scale)

    def get_global_x(self, scale=1):
        """ calculates the global telescope bins
        :param scale: telescope resolution in um
        :return: [nbins, bin_array] """
        xmax = max(self.PX * self.NCols, self.PY * self.NRows) * .6  # to keep the aspect ratio
        bins = arange(-xmax, xmax, self.PX * scale)
        return [bins.size - 1, bins]

    def get_global_y(self, scale=1):
        ymax = max(self.PX * self.NCols, self.PY * self.NRows) * .6
        bins = arange(-ymax, ymax, self.PY * scale)
        return [bins.size - 1, bins]

    # ----------------------------------------
    # region PIXEL
    def get_adc(self):
        return make_bins(0, self.MaxADC, bin_width=1)

    def get_vcal(self):
        return make_bins(self.MinVcal, self.MaxVcal, int(self.PHBinWidth / self.VcalToEl))

    def get_electrons(self, bin_width=None):
        return make_bins(self.MinPH, self.MaxPH, self.PHBinWidth if bin_width is None else bin_width)

    def get_ph(self, vcal=False, adc=False, bin_width=None):
        return self.get_vcal() if vcal else self.get_adc() if adc else self.get_electrons(bin_width)
    # endregion PIXEL
    # ----------------------------------------


def make_bins(min_val, max_val, bin_width):
    bins = arange(min_val, max_val + bin_width / 100., bin_width, dtype='d')
    return [bins.size - 1, bins]
