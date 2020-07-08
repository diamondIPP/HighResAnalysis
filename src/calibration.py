#!/usr/bin/env python
# --------------------------------------------------------
#       pulse height calibration of the cms pixel planes
# created on July 2nd 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import TF1, TGraph
from utils import *
from analysis import glob, join
from run import Run
from numpy import genfromtxt, split, all, delete, concatenate, round, zeros
from draw import Draw, format_histo, update_canvas
from copy import deepcopy


class Calibration:

    def __init__(self, run: Run, n=None):

        # INFO
        self.Run = run
        self.Plane = run.DUT.Plane
        self.Draw = Draw(config=self.Run.Config)

        # TODO add calibration in hdf5 file (no charge for MIMOSA anyway)
        # Calibration
        self.HighRangeFactor = 7
        self.Trim, self.Number = self.get_trim_number(n)
        self.Fit = TF1('ErFit', '[3] * (TMath::Erf((x - [0]) / [1]) + [2])', -500, 255 * 7)
        self.Points = None
        self.Fits = None

        self.PBar = PBar()
        self.correct_file()

    # ----------------------------------------
    # region INIT
    def get_trim_number(self, n=None):
        trim, number = [int(v) for v in self.Run.Logs['trim{}'.format(self.Run.DUT.Number)].split('-')]
        return trim, choose(n, number)

    def get_file_name(self):
        f = glob(join(self.Run.TCDir, 'calibrations', self.Run.DUT.Name, 'phCalibration{}*-{}.dat'.format(self.Trim, self.Number)))
        if len(f):
            return f[0]
        critical('Pulse height calibration file does not exist...')

    def get_fit_file(self):
        f = glob(join(self.Run.TCDir, 'calibrations', self.Run.DUT.Name, 'phCalibrationFitErr{}*-{}.dat'.format(self.Trim, self.Number)))
        if len(f):
            return f[0]
        warning('Pulse height calibration fit file does not exist...')

    def correct_file(self):
        """ correct corrupted files from pxar if there are more than 100 points... """
        data = genfromtxt(self.get_file_name(), skip_header=3, max_rows=1, dtype=object)[:-3].astype('u2')
        if data.size < 100:
            return True
        hr_size = self.get_vcals()['high range'].size
        if all(data[100:100 + hr_size] == data[-hr_size:]):
            warning('Found corrupted file from pxar. Correcting ...')
            low_range = self.get_vcals()['low range']
            lines = []
            data = self.read()
            i_corr = arange(100, 100 + hr_size)  # indices of the corrupted data
            with open(self.get_file_name(), 'r+') as f:
                for line in f.readlines():
                    if line.startswith('Low'):
                        line = 'Low range: {} \n'.format(' '.join(delete(low_range, i_corr).astype(str)))
                    lines.append(line)
                    if line.startswith('\n'):
                        break
                for col, rows in enumerate(data):
                    for row, lst in enumerate(rows):
                        # interpolate straight line from the good values around the corrupted one
                        missing = round([interpolate(99, 101 + hr_size, lst[99], lst[101 + hr_size], x) for x in i_corr]).astype('u2')
                        lst[-hr_size:] -= missing
                        lines.append('{} Pix {} {} \n'.format(' '.join(delete(lst, i_corr).astype(str)), col, row))
                f.seek(0)
                f.truncate(0)
                f.writelines(lines)
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def read(self):
        cols = arange(sum(a.size for a in self.get_vcals().values()))
        return array(split(genfromtxt(self.get_file_name(), skip_header=3, usecols=cols, dtype='u2'), self.get_splits()))

    def read_fit_pars(self):
        return array(split(genfromtxt(self.get_fit_file(), skip_header=3, usecols=arange(4), dtype='f2'), self.get_splits()))

    def get(self, col, row):
        if self.Points is None:
            self.Points = self.read()
        return self.Points[col][row]

    def get_all(self):
        if self.Points is None:
            self.Points = self.read()
        return self.Points

    def get_splits(self):
        """ :returns: array of indices when to split the calibration data -> at every new column and after NRows of the chip. """
        return arange(self.Plane.NRows, self.Plane.NCols * self.Plane.NRows, self.Plane.NRows)

    def get_vcals(self):
        """ :returns: the vcal dacs used for the calibration, which are saved in the header of the files. """
        f = self.get_file_name()
        return {'low range': genfromtxt(f, 'u2', skip_header=1, max_rows=1)[2:], 'high range': genfromtxt(f, 'u2', skip_header=2, max_rows=1)[2:]}

    def get_vcal_vec(self):
        return concatenate([v * f for v, f in zip(self.get_vcals().values(), [1, self.HighRangeFactor])])

    def get_charge(self, col, row, adc):
        fit = self.fit_erf(self.get_vcal_vec(), self.get(col, row)) if self.Fits is None else self.Fits[col][row]
        return -1 if fit is None else fit.GetX(adc)

    def get_chi2s(self):
        chi2s = zeros(self.Fits.shape, dtype='f2')
        for col, rows in enumerate(self.Fits):
            for row, fit in enumerate(rows):
                chi2s[col][row] = 1000 if fit is None else fit.GetChisquare() / fit.GetNDF()
        return chi2s
    # endregion GET
    # ----------------------------------------

    def load_fits(self, pbar=True):
        info('Fitting calibration points ...', prnt=pbar)
        do(self.PBar.start, self.Plane.NPixels, pbar)
        x = self.get_vcal_vec()
        pars = zeros((self.Plane.NCols, self.Plane.NRows), dtype=object)
        self.Fit.SetParameters(309.2062, 112.8961, 1.022439, 35.89524)
        for col, rows in enumerate(self.get_all()):
            for row, y in enumerate(rows):
                pars[col][row] = self.fit_erf(x, y)
                if pbar:
                    self.PBar.update()
        self.Fits = array(pars)

    def fit_erf(self, x, y):
        x, y = x[y > 0], y[y > 0]  # take only non zero values
        if x.size < 5:
            return
        g = TGraph(x.size, x.astype('d'), y.astype('d'))
        g.Fit(self.Fit, 'q0', '', 0, 255 * 7)
        return deepcopy(self.Fit)

    def draw_fit(self, col=14, row=14, show=True):
        self.Fit.SetParameters(309.2062, 112.8961, 1.022439, 35.89524)
        self.fit_erf(self.get_vcal_vec(), self.get(col, row))
        self.draw(col, row, show).SetTitle('Calibration Fit for Pix {} {}'.format(col, row))
        self.Fit.Draw('same')
        update_canvas()
        return self.Fit

    # ----------------------------------------
    # region DRAW
    def draw(self, col=14, row=14, show=True):
        g = self.Draw.make_tgrapherrors('gcal{}{}'.format(col, row), 'Calibration Points for Pixel {} {}'.format(col, row), x=self.get_vcal_vec(), y=self.get(col, row))
        format_histo(g, x_tit='vcal', y_tit='adc', y_off=1.4, markersize=.4)
        self.Draw.draw_histo(g, show, .12, draw_opt='ap')
        return g

    def draw_calibration_fit(self, col=14, row=14, show=True):
        self.Fit.SetParameters(*self.read_fit_pars()[col][row])
        self.draw(col, row, show).SetTitle('Calibration Fit for Pixel {} {}'.format(col, row))
        self.Fit.Draw('same')
    # endregion DRAW
    # ----------------------------------------
