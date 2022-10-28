#!/usr/bin/env python
# --------------------------------------------------------
#       class for analysis of a batch of runs for a single DUT
# created on October 27th 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from src.dut_analysis import DUTAnalysis, Analysis, ev2str
from src.run import Batch
from src.converter import batch_converter


class BatchAnalysis(DUTAnalysis):

    def __init__(self, batch_name, dut_number, test_campaign, verbose=True, test=False):

        self.Batch = Batch(batch_name, dut_number, Analysis(test_campaign).BeamTest.Path)
        super().__init__(self.Batch.min_run.Number, dut_number, test_campaign, verbose, test)

    def __repr__(self):
        return f'{self} of batch {self.Batch} ({self.BeamTest}), {self.ev_str}'

    @property
    def ev_str(self):
        return f'{ev2str(self.NEvents if hasattr(self, "NEvents") else self.Batch.n_ev)} ev'

    def init_converter(self):
        return self.converter.from_batch(self.Batch)

    @property
    def converter(self):
        return batch_converter(super(BatchAnalysis, self).converter)

    @property
    def file_name(self):
        return self.Batch.FileName
