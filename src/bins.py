#!/usr/bin/env python
# --------------------------------------------------------
#       class for creating the histogram bins for the analysis
# created on June 10th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from numpy import arange, append, diff


# Pixel
MaxADC = 2**8 - 1
MinPH = -5000
MaxPH = 65000
PHBinWidth = 200
MinVcal = -100
MaxVcal = 1250


def get_time(t_vec, bin_width, last=False):
    return make(t_vec[0], t_vec[-1], bin_width, last)


# ----------------------------------------
# region PIXEL
def get_local(plane, bin_width=1, aspect_ratio=False):
    return get_local_x(plane, bin_width, aspect_ratio) + get_local_y(plane, bin_width, aspect_ratio)


def get_local_x(plane, bin_width=1, aspect_ratio=False):
    extra_pixel = round((plane.get_max_width() - plane.get_x_width()) / plane.PX / 2) if aspect_ratio else 0  # keep aspect ratio
    return make(-extra_pixel - .5, plane.NCols + extra_pixel - .5, bin_width)


def get_local_y(plane, bin_width=1, aspect_ratio=False):
    extra_pixel = round((plane.get_max_width() - plane.get_y_width()) / plane.PY / 2) if aspect_ratio else 0  # keep aspect ratio
    return make(-extra_pixel - .5, plane.NRows + extra_pixel - .5, bin_width)


def get_corr(mode, pl0, pl1, bw=1):
    return sum([get_local_x(pl, bw) if mode.lower() == 'x' else get_local_y(pl, bw) for pl in [pl0, pl1]], start=[])


def get_global(plane, res=1):
    return get_global_x(plane, res) + get_global_y(plane, res)


def get_xy(local, plane, bin_width=1, aspect_ratio=False):
    return get_local(plane, bin_width, aspect_ratio) if local else get_global(plane, bin_width)


def get_global_x(plane, res=1):
    """ calculates the global telescope bins
    :return: [nbins, bin_array] """
    xmax = plane.get_max_width() * .6  # to keep the aspect ratio
    return make(-xmax, xmax, res * plane.PX)


def get_x(plane, bw=1, res=1, local=True, aspect_ratio=False):
    return get_local_x(plane, bw, aspect_ratio) if local else get_global_x(plane, res)


def get_y(plane, bw=1, res=1, local=True, aspect_ratio=False):
    return get_local_y(plane, bw, aspect_ratio) if local else get_global_y(plane, res)


def get_global_y(plane, res=1):
    ymax = plane.get_max_width() * .6
    return make(-ymax, ymax, res * plane.PY)


def get_pixel(plane, res, outer=.5, cell=False):
    x0 = outer if plane.PX > plane.PY or cell else (2 * outer + 1) * (plane.R - 1) / 2 + outer
    y0 = outer if plane.PY > plane.PX or cell else (2 * outer + 1) * (plane.R - 1) / 2 + outer
    return make(-x0, 1 + x0, res, last=True) + make(-y0, 1 + y0, res, last=True)


def get_adc():
    return make(0, MaxADC)


def get_vcal(bin_width=1):
    return make(MinVcal, MaxVcal, bin_width)


def get_electrons(bin_width=200):
    return make(MinPH, MaxPH, bin_width)


def get_ph(vcal=False, adc=False, bin_width=None):
    return get_vcal() if vcal else get_adc() if adc else get_electrons(bin_width)


def get_triggerphase():
    return make(-.5, 10.5)
# endregion PIXEL
# ----------------------------------------


def make(min_val, max_val, bin_width=1, last=False):
    bins = arange(min_val, max_val + (bin_width if last else 0), bin_width, dtype='d')
    return [bins.size - 1, bins]


def make2d(x, y):
    return make(x[0], x[-1], x[1] - x[0], last=True) + make(y[0], y[-1], y[1] - y[0], last=True)


def from_vec(x):
    return [x.size, append([i.n - i.s for i in x], x[-1].n + x[-1].s).astype('d')]


def from_p(x):
    d = x[-1] - x[-2]
    x = append(x, [x[-1] + d, x[-1] + 2 * d])
    return [x.size - 2, x[:-1] - diff(x) / 2]


# Binning
TP = get_triggerphase()
