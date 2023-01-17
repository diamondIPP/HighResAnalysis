#!/usr/bin/env python
# --------------------------------------------------------
#       adds clustering and charge to trees created with pXar
# created on August 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from subprocess import check_call
import uproot

from plotting.utils import array, info as pinfo
from src.converter import Converter, DUT


class CERNConverter(Converter):
    """
    Converts the raw data taken at CERN in several steps to hdf5 files.
    The raw data consists of one root file for the DUTs saved with pXar and a binary file from the KARTEL telescope.
    STEP -2: raw -> root for tel          (judith)
    STEP -1: adc -> vcal for DUTs         (python)
    STEP  0: merge tel and dut root files (python)\n"""
    __doc__ += '\n'.join(Converter.__doc__.split('\n')[3:])

    def __init__(self, data_dir, run_number, dut_name=None):

        self.DUTName = dut_name
        Converter.__init__(self, data_dir, run_number)

        self.Adc2Vcal = self.init_adc2vcal()
        self.EventAlignment = self.init_event_alignment()
        self.Ref = self.init_ref()

    def proteus_raw_file_path(self):
        return self.DataDir.joinpath('data', f'run{self.Run:04d}.root')

    def trigger_info_file(self):
        return self.Adc2Vcal.OutFilePath

    @property
    def first_steps(self):
        return self.Raw.Steps + self.Adc2Vcal.Steps + self.Ref.Steps + [(self.merge_root_files, self.proteus_raw_file_path())]

    @property
    def aux_files(self):
        return [self.Raw.OutFilePath, self.Adc2Vcal.OutFilePath, self.Ref.OutFilePath] + super().aux_files

    @property
    def raw_files(self):
        return [self.EventAlignment.RefPath, self.Adc2Vcal.RawFilePath] + super().raw_files

    def init_raw(self):
        from cern.raw import CERNRaw
        return CERNRaw(self)

    def init_event_alignment(self):
        from cern.event_alignment import EventAlignment
        return EventAlignment(self.Raw)

    def init_adc2vcal(self):
        from cern.adc import Adc2Vcal
        return Adc2Vcal(self)

    def init_ref(self):
        from cern.ref import RefConverter
        return RefConverter(self)

    def init_duts(self):
        return super().init_duts() if self.DUTName is None else [DUT.from_name(self.DUTName, self.Run.Info)]

    def merge_root_files(self, force=False):
        """merge the telescope and DUT root files"""
        self.OutFilePath.parent.mkdir(exist_ok=True)
        cmd = f'hadd{" -f" if force else ""} {self.proteus_raw_file_path()} {self.Raw.OutFilePath} {self.Ref.OutFilePath} {self.Adc2Vcal.OutFilePath}'
        pinfo(cmd)
        check_call(cmd, shell=True)

    @property
    def time_stamp_file(self):
        return self.Adc2Vcal.OutFilePath

    def get_time_stamp(self):
        return (array(uproot.open(self.time_stamp_file)['Event']['TimeStamp']) / 1000).astype('d')

    @property
    def calibration(self):
        from cern.calibration import CERNCalibration
        return CERNCalibration
