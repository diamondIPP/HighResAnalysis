#!/usr/bin/env python
# --------------------------------------------------------
#       Analysis of the reference plane
# created on August 19th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from src.dut_analysis import DUTAnalysis, Run
from mod.residuals import ResidualAnalysis
from mod.ref_cuts import RefCut
from src.calibration import Calibration


class RefAnalysis(DUTAnalysis):

    def __init__(self, parent):  # noqa

        self.IsRef = True
        self.Parent = parent
        self.REF = self.Parent  # make DUTAna ref of the reference
        self.__dict__.update(parent.__dict__)
        self.Plane = self.Planes[self.Config.getint('DUT', 'reference plane')]
        self.N = self.n
        self.MetaSubDir = 'REF'

        self.Run = Run(parent.Run.Number, self.get_dut_number(), self.BeamTest.Path, parent.Run.SingleMode)
        self.DUT = self.Run.DUT
        self.Calibration = self.Converter.load_calibration(self.Run.DUT.Number)
        self.Cut = RefCut(self)

        self.Residuals = ResidualAnalysis(self)
        self.Tracks = self.init_tracks()
        self.Efficiency = self.init_eff()
        self.Cut.make_additional()

    def get_dut_number(self):
        default = next(i for i, dut in enumerate(self.Parent.Run.Logs['duts']) if dut != self.Parent.DUT.Name)
        return next((i for i, dut in enumerate(self.Parent.Run.Logs['duts']) if dut.startswith('Si') or dut.startswith('D')), default)
