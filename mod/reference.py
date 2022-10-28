#!/usr/bin/env python
# --------------------------------------------------------
#       Analysis of the reference plane
# created on August 19th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from src.dut_analysis import Run
from src.dut import REF
from mod.residuals import ResidualAnalysis
from mod.ref_cuts import RefCut


def ref_analysis(cls):
    class RefAnalysis(cls):

        def __init__(self, parent):  # noqa

            self.IsRef = True
            self.Parent = parent
            self.REF = self.Parent  # make DUTAna ref of the reference
            self.__dict__.update(parent.__dict__)
            self.Plane = self.Planes[self.Config.getint('DUT', 'reference plane')]
            self.N = self.n
            self.MetaSubDir = 'REF'

            self.Run = Run(parent.Run.Number, self.dut_nr, self.BeamTest.Path)
            self.DUT = REF(self.dut_nr) if self.Parent.Proteus.NRefPlanes else self.Run.DUT
            self.Calibration = self.Converter.load_calibration(self.Run.DUT.Number)
            self.Cut = RefCut(self)

            self.Residuals = ResidualAnalysis(self)
            self.Tracks = self.init_tracks()
            self.Efficiency = self.init_eff()
            self.Cut.make_additional()

        @property
        def dut_nr(self):
            default = next((i for i, dut in enumerate(self.Parent.Run.Info['duts']) if dut != self.Parent.DUT.Name), 0)
            return 0 if self.Parent.Proteus.NRefPlanes else next((i for i, dut in enumerate(self.Parent.Run.Info['duts']) if dut.startswith('Si') or dut.startswith('D')), default)

    return RefAnalysis
