#!/usr/bin/env python
# --------------------------------------------------------
#       class for creating the histogram bins for the analysis
# created on June 10th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from numpy import arange


# Pixel
MaxADC = 2**8 - 1
MinPH = -5000
MaxPH = 100000
PHBinWidth = 200
MinVcal = -100
MaxVcal = 1250
VcalToEl = 47.

# Miscellaneous
GlobalCoods = [-.5025, .5175, -.505, .515]  # range in x and y in telescope coordinates [cm]


def get_pixel(self, plane):
    return self.get_pixel_x(plane) + self.get_pixel_y(plane)


def get_pixel_x(plane):
    extra_pixel = (plane.get_max_width() - plane.get_x_width()) / plane.PX / 2  # keep aspect ratio
    return make(-extra_pixel - .5, plane.NCols + extra_pixel - .5)


def get_pixel_y(plane):
    extra_pixel = (plane.get_max_width() - plane.get_y_width()) / plane.PY / 2  # keep aspect ratio
    return make(-extra_pixel - .5, plane.NRows + extra_pixel - .5)


# ----------------------------------------
# region PIXEL
def get_adc():
    return make(0, MaxADC)


def get_vcal():
    return make(MinVcal, MaxVcal, int(PHBinWidth / VcalToEl))


def get_electrons(bin_width=None):
    return make(MinPH, MaxPH, PHBinWidth if bin_width is None else bin_width)


def get_ph(vcal=False, adc=False, bin_width=None):
    return get_vcal() if vcal else get_adc() if adc else get_electrons(bin_width)
# endregion PIXEL
# ----------------------------------------


def make(min_val, max_val, bin_width=1):
    bins = arange(min_val, max_val + bin_width / 100., bin_width, dtype='d')
    return [bins.size - 1, bins]
