#!/usr/bin/env python
# --------------------------------------------------------
#       raw file sub class for CERN
# created on April 22nd 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from src.raw import Raw
from src.converter import Converter, Analysis


class CERNRaw(Raw):

    def __init__(self, c: Converter, load_file=False, step=-1):
        super().__init__(c, load_file, step)

    def load_raw_file_path(self):
        return self.DataDir.joinpath('tel', f'acq{self.Run.Info["telescope run"]:03d}.bin')

    def load_out_file_path(self):
        return self.RawFilePath.with_name(f'tel-run{self.Run:04d}.root')

    def load_soft_dir(self):
        return self.Parent.SoftDir.joinpath(Analysis.Config.get('SOFTWARE', 'judith'))

    @property
    def soft(self):
        return self.SoftDir.joinpath('Judith')

    @property
    def options(self):
        return f'-c convert -m {self.SoftDir.joinpath("configs", "readout", "CERN.cfg")}'

    def convert(self):
        """convert binary raw file (from KARTEL telescope) to root file with judith."""
        super().convert()
