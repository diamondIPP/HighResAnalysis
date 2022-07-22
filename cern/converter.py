#!/usr/bin/env python
# --------------------------------------------------------
#       adds clustering and charge to trees created with pXar
# created on August 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from src.converter import Converter
from subprocess import check_call
from plotting.utils import info as pinfo


class CERNConverter(Converter):
    """
    Converts the raw data taken at CERN in several steps to hdf5 files.
    The raw data consists of one root file for the DUTs saved with pXar and a binary file from the KARTEL telescope.
    STEP -2: adc -> vcal for DUTs         (python)
    STEP -1: raw -> root for tel          (judith)
    STEP  0: merge tel and dut root files (python)\n"""
    __doc__ += '\n'.join(Converter.__doc__.split('\n')[3:])

    def __init__(self, data_dir, run_number):

        Converter.__init__(self, data_dir, run_number)
        self.Adc2Vcal = self.init_adc2vcal()

    def proteus_raw_file_path(self):
        return self.DataDir.joinpath('data', f'run{self.Run:04d}.root')

    def trigger_info_file(self):
        return self.Adc2Vcal.OutFilePath

    @property
    def first_steps(self):
        return self.Adc2Vcal.Steps + self.Raw.Steps + [(self.merge_root_files, self.proteus_raw_file_path())]

    @property
    def raw_files(self):
        return [self.Raw.OutFilePath, self.Adc2Vcal.OutFilePath] + super().raw_files

    def init_raw(self):
        from cern.raw import CERNRaw
        return CERNRaw(self)

    def init_adc2vcal(self):
        from cern.adc import Adc2Vcal
        return Adc2Vcal(self)

    def merge_root_files(self, force=False):
        """merge the telescope and DUT root files"""
        cmd = f'hadd {"-f" if force else ""} {self.proteus_raw_file_path()} {self.Raw.OutFilePath} {self.Adc2Vcal.OutFilePath}'
        pinfo(cmd)
        check_call(cmd, shell=True)


if __name__ == '__main__':

    from argparse import ArgumentParser
    import uproot, awkward as aw  # noqa
    from numpy import *

    p = ArgumentParser()
    p.add_argument('run', nargs='?', default=128, type=int)
    p.add_argument('dut', nargs='?', default=1, type=int)
    args = p.parse_args()

    z = CERNConverter.from_ana(args.run, args.dut, single_mode=True)
    run = z.Run
    r = z.Raw
    c = z.load_calibration()
    p = z.Proteus
    draw = c.Draw
    adc = z.Adc2Vcal
