#!/usr/bin/env python
# --------------------------------------------------------
#       pulse height calibration of the cms pixel planes
# created on July 2nd 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import TF1
from utils import *
from analysis import basename, glob, join
from pickle import load as pload
from run import Run


class Calibration:

    def __init__(self, run: Run):

        self.Run = run

        # TODO add calibration in hdf5 file (no charge for MIMOSA anyway)
        # Calibration
        self.File = self.load_file()
        self.Fit = TF1('ErFit', '[3] * (TMath::Erf((x - [0]) / [1]) + [2])', -500, 255 * 7)
        # self.FitParameters = None
        # self.get_calibration_data()


    # ----------------------------------------
    # region INIT
    def load_file(self):
        pass

    # endregion INIT
    # ----------------------------------------

    def get_calibration_number(self):
        numbers = sorted(int(remove_letters(basename(name))) for name in glob(join(self.TCDir, self.DUT.Name, 'calibrations', 'phCal*.dat')) if basename(name)[5].isdigit())
        first_run = int(self.Run.Logs.keys()[0])
        if first_run > numbers[-1]:
            return numbers[-1]
        next_number = next(nr for nr in numbers if nr >= first_run)
        return numbers[numbers.index(next_number) - 1]

    def get_calibration_data(self):
        pickle_name = join(self.TCDir, self.DUT.Name, 'calibrations', 'fitpars{}.pickle'.format(self.get_calibration_number()))
        with open(pickle_name, 'r') as f:
            self.FitParameters = pload(f)