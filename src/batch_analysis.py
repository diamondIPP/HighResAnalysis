#!/usr/bin/env python
# --------------------------------------------------------
#       class for analysis of a batch of runs for a single DUT
# created on October 27th 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from src.dut_analysis import DUTAnalysis, Analysis, ev2str, Path, datetime
from src.run import Batch
from src.converter import batch_converter


class BatchAnalysis(DUTAnalysis):

    def __init__(self, batch_name, dut_number, test_campaign, verbose=True, test=False):

        self.Batch = batch_name if type(batch_name) is Batch else Batch(batch_name, dut_number, Analysis(test_campaign).BeamTest.Path)
        self.DUT = self.Batch.DUT
        super().__init__(self.Batch.min_run.Number, dut_number, test_campaign, verbose, test)

    @classmethod
    def from_batch(cls, batch: Batch, verbose=True, test=False):
        return cls(batch, batch.DUT.Number, test_campaign=batch.DataDir.stem, verbose=verbose, test=test)

    @property
    def server_save_dir(self):
        return Path('duts', str(self.DUT), self.BeamTest.Tag, f'b-{self.Batch}')

    @property
    def ev_str(self):
        return f'{ev2str(self.NEvents if hasattr(self, "NEvents") else self.Batch.n_ev)} ev'

    @property
    def suffix(self):
        return f'{self.DUT}-{self.Batch}-{self.BeamTest.Location}'.lower().replace('ii6-', '')

    @property
    def unit_str(self):
        return f'batch {self.Batch}'

    @property
    def run_str(self):
        return f'b-{self.Batch}'

    def init_converter(self):
        return self.converter.from_batch(self.Batch)

    @property
    def converter(self):
        return batch_converter(super(BatchAnalysis, self).converter)

    @property
    def file_name(self):
        return self.Batch.FileName

    def get_end_time(self):
        return datetime.fromtimestamp(self.Batch.Runs[-1].EndTime)
