#!/usr/bin/env python
# --------------------------------------------------------
#       raw file conversion & analysis
# created on April 22nd 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from src.converter import Converter
from plotting.draw import critical, Draw, info, remove_file, Path
from utility.utils import get_tree_vec
from ROOT import TFile
from subprocess import check_call
from numpy import all


class Raw:

    def __init__(self, c: Converter, load_file=False):

        self.Parent = c
        self.RunNumber = c.RunNumber

        self.SoftDir = c.SoftDir.joinpath(c.Config.get('SOFTWARE', 'eudaq2'))
        self.DataDir = c.DataDir
        self.SaveDir = c.SaveDir

        self.FilePath = self.load_file_path()
        self.OutFilePath = c.SaveDir.joinpath(f'run{self.RunNumber:06d}.root')

        if load_file:
            self.F = TFile(str(self.OutFilePath)) if load_file else None
            self.P = [self.F.Get(key.GetName()).Get('Hits') for key in self.F.GetListOfKeys() if key.GetName().startswith('Plane')]
        self.Draw = Draw(c.Config.FilePath)

    def __repr__(self):
        return f'Raw file analysis run {self.RunNumber} ({self.FilePath.name})'

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

    def make_cuts(self):
        [t.SetEstimate(t.GetEntries() * 10) for t in self.P]
        n = [get_tree_vec(t, 'NHits', dtype='i2') for t in self.P[:6]]
        c = all([i == 1 for i in n], axis=0)
        return [c.repeat(i) for i in n]


if __name__ == '__main__':
    from src.analysis import Analysis

    a = Analysis()
    c0 = Converter(a.BeamTest.Path, 17, a.Config)
    z = Raw(c0, load_file=True)
