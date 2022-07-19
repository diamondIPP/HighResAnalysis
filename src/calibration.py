#!/usr/bin/env python
# --------------------------------------------------------
#       pulse height calibration of the CMS pixel planes
# created on July 2nd 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import TGraph
from numpy import genfromtxt, all, delete, round, argmax, savetxt

from plotting.draw import Draw, FitRes
from plotting.fit import Erf
from src.run import Run
from src.analysis import Analysis
import src.bins as bins
from utility.utils import *


class Calibration:

    def __init__(self, run: Run, n=None):

        # INFO
        self.Run = run
        self.DUT = run.DUT
        self.Plane = run.DUT.Plane
        self.NX, self.NY, self.NPix = self.Plane.NCols, self.Plane.NRows, self.Plane.NPixels
        self.Draw = Draw(config=Analysis.Config)
        self.Dir = self.Run.TCDir.joinpath('calibrations', str(self.Run.DUT))

        # Calibration
        self.HighRangeFactor = 7
        self.Trim, self.Number = self.get_trim_number(n)
        self.Fit = Erf(fit_range=[0, 255 * 7]).Fit
        self.Points = None
        self.Fits = None

        self.FileName = join(self.Dir, f'{self.Trim}-{self.Number}.h5py')
        self.FitFileName = self.Dir.joinpath(f'fitpars-{self.Trim}.txt')
        self.CalPath = self.Dir.parent.joinpath('fitpars-.txt')  # enter trim and DUT in eudaq

        self.PBar = PBar()
        self.correct_file()

    def __call__(self, x, y, adc):
        return self.get_lookup_table()[x, y, adc]

    def __repr__(self):
        return f'ADC Calibration of {self.Run.DUT}'

    # ----------------------------------------
    # region INIT
    def get_trim_number(self, n=None):
        if f'trim{self.DUT.Number}' not in self.Run.Logs:
            return None, None
        trim, number = [int(v) for v in self.Run.Logs[f'trim{self.Run.DUT.Number}'].split('-')]
        return trim, choose(n, number)

    def get_file_name(self):
        trim, n = choose('', self.Trim), '' if self.Number is None else f'-{self.Number}'
        f = list(self.Dir.glob(f'phCalibration{trim}*{n}.dat'))
        info(f'reading calibration file {f[0].name} for {self.DUT}') if f else do_nothing()
        return f[0] if f else critical(f'Pulse height calibration file {"" if trim is None else f"{trim}-{n} "}does not exist in {self.Dir} ...')

    def get_fit_file(self):
        f = list(self.Dir.glob(f'phCalibrationFitErr{self.Trim}*-{self.Number}.dat'))
        return f[0] if len(f) else warning('Pulse height calibration fit file does not exist...')

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
            data = self.get_all_points()
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
                        lines.append(f'{" ".join(delete(lst, i_corr).astype(str))} Pix {col} {row} \n')
                f.seek(0)
                f.truncate(0)
                f.writelines(lines)
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def make_hdf5_path(self, name='', suf=''):
        d = ensure_dir(Dir.joinpath('metadata', 'calibration'))
        return d.joinpath(f'{"_".join(str(v) for v in [name, self.Run.DUT, self.Trim, self.Number, suf] if v)}.hdf5')

    @save_hdf5('Points', arr=True)
    def get_all_points(self):
        return genfromtxt(self.get_file_name(), skip_header=3, dtype='u2')[:, :-3].reshape((self.NX, self.NY, -1))   # last three entries are pixel info

    def get_thresholds(self):
        return self.vcals[argmax(self.get_all_points() > 0, axis=2)]

    def read_fit_pars(self):
        return genfromtxt(self.get_fit_file(), skip_header=3, usecols=arange(4), dtype='f2').reshape((self.NX, self.NY, -1))

    def get_formula(self):
        return str(genfromtxt(self.get_fit_file(), skip_header=1, max_rows=1, dtype=str)).replace('par', '').replace('x[0]', 'x')

    def get_points(self, col, row):
        return self.get_all_points()[col, row]

    def get_vcals(self):
        """ :returns: the vcal dacs used for the calibration, which are saved in the header of the files. """
        f = self.get_file_name()
        return {'low range': genfromtxt(f, 'u2', skip_header=1, max_rows=1)[2:], 'high range': genfromtxt(f, 'u2', skip_header=2, max_rows=1)[2:]}

    @property
    def vcals(self):
        return concatenate([v * f for v, f in zip(self.get_vcals().values(), [1, self.HighRangeFactor])])
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region FIT
    @update_pbar
    def fit(self, x, y):
        x, y = x[y > 0], y[y > 0]  # take only non zero values
        if x.size < 5:
            return None
        self.Fit.SetParameters(255 / 2, 255 / 2, 400, 500)
        TGraph(x.size, x.astype('d'), y.astype('d')).Fit(self.Fit, 'q0', '', 0, 255 * 7)
        return FitRes(self.Fit)

    def fit_all(self):
        info('fit calibration points ...')
        x, y = self.vcals, self.get_all_points()
        PBAR.start(self.NPix)
        self.Fits = [[self.fit(x, iy) for iy in lst] for lst in y]
        return self.Fits

    def get_fits(self):
        return self.fit_all() if self.Fits is None else self.Fits

    @save_hdf5('Chi2', arr=True, dtype='f4')
    def get_chi2s(self, _redo=False):
        return array([[1000 if fit is None else fit.get_chi2() for fit in lst] for lst in self.get_fits()])

    @update_pbar
    def get_vcal(self, f, adc):
        return f.GetX(adc) if adc > f.GetMinimum() else 0

    @save_hdf5('LUT', arr=True, dtype='f2')
    def get_lookup_table(self, _redo=False):
        fits = self.get_fits()
        info('creating calibration LUT ... ')
        PBAR.start(256 * self.NPix)
        return array([[[0 if f is None else self.get_vcal(f.Fit, i) for i in range(256)] for f in lst] for lst in fits])

    def save_fit_pars(self):
        """ [0] + [1] * TMath::Erf((x - [2]) / [3]) -> [3] * (TMath::Erf((x - [0]) / [1]) + [2])"""
        info(f'generating calibration fit parameters for {self.Run.DUT} ...')
        pars = array([[[0, 0, 0, 1] if fit is None else fit.get_pars(err=False)[[2, 3, 0, 1]] for fit in lst] for lst in self.get_fits()])  # fix order
        pars[:, :, 2] /= pars[:, :, 3]  # fix offset
        lines = [' '.join(f'{p:.6e}' for p in ps) + f' Pix {col:2d} {row:2d}' for col, lst in enumerate(pars) for row, ps in enumerate(lst)]
        header = list(genfromtxt(self.get_fit_file(), dtype=str, max_rows=2, delimiter='abc')) + ['']
        savetxt(self.FitFileName, header + lines, fmt='%s')  # noqa
    # endregion FIT
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw(self, col=14, row=14, **dkw):
        x, y = self.vcals, self.get_points(col, row)
        return self.Draw.graph(x, y, **prep_kw(dkw, title=f'Calibration Points for Pixel {col} {row}', draw_opt='ap', x_tit='VCAL', y_tit='ADC', markersize=.7))

    def draw_fit(self, col=14, row=14, **dkw):
        f = self.fit(self.vcals, self.get_points(col, row))
        return self.draw(col, row, **prep_kw(dkw, title=f'Calibration Fit for Pix {col} {row}', leg=f.Fit))

    def draw_pxar_fit(self, col=14, row=14, **dkw):
        """ draws the Erf fit from pXar """
        f = Draw.make_f(None, self.get_formula(), 0, 255 * 7, pars=self.read_fit_pars()[col, row])
        return self.draw(col, row, **prep_kw(dkw, title=f'Calibration Fit for Pixel {col} {row}', leg=f))

    def draw_thresholds(self, **dkw):
        return self.Draw.prof2d(self.get_thresholds(), binning=bins.get_local(self.Plane), **prep_kw(dkw, x_tit='Column', y_tit='Row', z_tit='Threshold [vcal]'))

    def draw_chi2_map(self, **dkw):
        return self.Draw.prof2d(self.get_chi2s(), binning=bins.get_local(self.Plane), **prep_kw(dkw, x_tit='Column', y_tit='Row', z_tit='#chi^{2}'))

    def draw_chi2(self, **dkw):
        return self.Draw.distribution(self.get_chi2s().flatten(), **prep_kw(dkw, rf=1, x0=0, x_tit='#chi^{2}'))
    # endregion DRAW
    # ----------------------------------------
