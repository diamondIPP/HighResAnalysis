#!/usr/bin/env python
# --------------------------------------------------------
#       raw file conversion & analysis
# created on April 22nd 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from subprocess import check_call
import h5py
from numpy import all, array

from plotting.draw import critical, Draw, info, remove_file, Path
from src.converter import Converter


class Raw:

    def __init__(self, c: Converter, load_file=False):

        self.Parent = c
        self.RunNumber = c.RunNumber
        self.NT = c.NTelPlanes
        self.ND = c.NDUTPlanes
        self.P = range(self.NT)

        self.SoftDir = c.SoftDir.joinpath(c.Config.get('SOFTWARE', 'eudaq2'))
        self.DataDir = c.DataDir
        self.SaveDir = c.SaveDir

        self.FilePath = self.load_file_path()
        self.OutFilePath = c.SaveDir.joinpath(f'run{self.RunNumber:06d}.root')

        if load_file:
            self.F = c.F if c.F is not None else h5py.File(c.OutFileName, 'r')

        self.Steps = [(self.convert, self.OutFilePath)]
        self.Draw = Draw(c.Config.FilePath)

    def __repr__(self):
        return f'Raw file analysis run {self.RunNumber} ({self.FilePath.name})'

    # ----------------------------------------
    # region CONVERT
    def load_file_path(self):
        n = list(self.DataDir.joinpath('raw').glob(f'run{self.RunNumber:06d}*.raw'))
        return n[0] if len(n) else None

    def generate_fit_files(self):
        c = None
        for dut in range(self.Parent.NDUTPlanes):
            c = self.Parent.get_calibration(dut)
            if not c.FitFileName.exists():
                c.save_fit_pars()
        return c.CalPath

    def convert(self):
        """ convert binary raw file to root file with eudaq"""
        cal_path = self.generate_fit_files()
        if self.FilePath is None:
            critical('raw file does not exist for run: {}'.format(self.RunNumber))
        self.OutFilePath.parent.mkdir(exist_ok=True)
        cmd = f'{self.SoftDir.joinpath("bin", "euCliConverter")} -i {self.FilePath} -o {self.OutFilePath} -c {cal_path}'
        info(f'Convert {self.FilePath.name} to {self.OutFilePath.name} using EUDAQ-2\n')
        info(f'{cmd}\n')
        check_call(cmd, shell=True)
        for f in Path().glob('AutoDict_vector*'):
            remove_file(f)
    # endregion CONVERT
    # ----------------------------------------

    # ----------------------------------------
    # region ANALYSIS
    def get(self, p, g, k):
        return array(self.F[f'Plane{p}'][g][k])

    def z(self):
        return self.Parent.Proteus.get_z_positions(raw=True)[:self.NT]

    def get_all(self, g, k):
        return array([self.get(p, g, k)[c] for p, c in zip(self.P, self.cuts())])

    def cuts(self):
        n = [self.get(p, 'Clusters', 'Size') > 0 for p in self.P]
        c = all(n, axis=0)
        return [c[i] for i in n]

    def draw_track_x(self, i=0):
        self.Draw.graph(self.z(), self.get_all('Clusters', 'U').T[i], x_tit='Z [mm]', y_tit='U [mm]')

    def draw_track_y(self, i=0):
        self.Draw.graph(self.z(), self.get_all('Clusters', 'V').T[i], x_tit='Z [mm]', y_tit='V [mm]')
    # ----------------------------------------
    # endregion ANALYSIS


if __name__ == '__main__':
    from src.analysis import Analysis

    a = Analysis()
    c_ = Converter(a.BeamTest.Path, 11, a.Config)
    z = Raw(c_, load_file=True)
