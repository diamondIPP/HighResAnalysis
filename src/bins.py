#!/usr/bin/env python
# --------------------------------------------------------
#       class for creating the histogram bins for the analysis
# created on June 10th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from numpy import arange, concatenate


# Pixel
MaxADC = 2**8 - 1
MinPH = -5000
MaxPH = 50000
PHBinWidth = 200
MinVcal = -100
MaxVcal = 1250


def get_time(t_vec, bin_width, last=False):
    return make(t_vec[0], t_vec[-1], bin_width, last)


# ----------------------------------------
# region PIXEL
def get_local(plane, bin_width=1):
    return get_local_x(plane, bin_width) + get_local_y(plane, bin_width)


def get_local_x(plane, bin_width=1, aspect_ratio=False):
    extra_pixel = round((plane.get_max_width() - plane.get_x_width()) / plane.PX / 2) if aspect_ratio else 0  # keep aspect ratio
    return make(-extra_pixel - .5, plane.NCols + extra_pixel - .5, bin_width)


def get_local_y(plane, bin_width=1, aspect_ratio=False):
    extra_pixel = round((plane.get_max_width() - plane.get_y_width()) / plane.PY / 2) if aspect_ratio else 0  # keep aspect ratio
    return make(-extra_pixel - .5, plane.NRows + extra_pixel - .5, bin_width)


def get_corr(mode, plane0, plane1):
    return list(concatenate([get_local_x(pl) if mode.lower() == 'x' else get_local_y(pl) for pl in [plane0, plane1]]))


def get_global(plane, res=1):
    return get_global_x(plane, res) + get_global_y(plane, res)


def get_global_x(plane, res=1):
    """ calculates the global telescope bins
    :return: [nbins, bin_array] """
    xmax = plane.get_max_width() * .6  # to keep the aspect ratio
    return make(-xmax, xmax, res * plane.PX)


def get_global_y(plane, res=1):
    ymax = plane.get_max_width() * .6
    return make(-ymax, ymax, res * plane.PY)


def get_adc():
    return make(0, MaxADC)


def get_vcal(bin_width=1):
    return make(MinVcal, MaxVcal, bin_width)


def get_electrons(bin_width=200):
    return make(MinPH, MaxPH, bin_width)


def get_ph(vcal=False, adc=False, bin_width=None):
    return get_vcal() if vcal else get_adc() if adc else get_electrons(bin_width)
# endregion PIXEL
# ----------------------------------------


def make(min_val, max_val, bin_width=1, last=False):
    bins = arange(min_val, max_val + (bin_width if last else 0), bin_width, dtype='d')
    return [bins.size - 1, bins]
