#!/usr/bin/env python
# --------------------------------------------------------
#       Analysis of the reference plane
# created on August 19th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from src.dut_analysis import Run
from src.dut import REF
from mod.residuals import res_analysis
from mod.ref_cuts import RefCut


def ref_analysis(cls):
    class RefAnalysis(cls):

        def __init__(self, parent):  # noqa

            self.IsRef = True
            self.Parent = parent
            self.REF = self.Parent  # make DUTAna ref of the reference
            self.__dict__.update(parent.__dict__)
            self.N = self.n
            self.MetaSubDir = 'REF'

            self.Run = Run(parent.Run.Number, self.dut_nr, self.BeamTest.Path)
            self.DUT = REF(self.dut_nr) if self.Parent.Proteus.NRefPlanes else self.Run.DUT
            self.Plane = self.Planes[self.Config.getint('DUT', 'reference plane')] if self.Parent.Proteus.NRefPlanes > 0 else self.DUT.Plane
            self.Calibration = None if self.Parent.Proteus.NRefPlanes else self.Converter.load_calibration(self.Run.DUT.Number)
            self.Cut = RefCut(self)

            self.Residuals = res_analysis(cls)(self)
            self.Tracks = self.init_tracks()
            self.Efficiency = self.init_eff()
            self.Cut.make_additional()

        @property
        def dut_nr(self):
            default = next((i for i, dut in enumerate(self.Parent.Run.DUTs) if dut != self.Parent.DUT.Name), 0)
            is_ref = lambda dut: dut in self.Config.get('DUT', 'reference detectors') and dut != self.Parent.DUT.Name
            return 0 if self.Parent.Proteus.NRefPlanes > 0 else next((i for i, dut in enumerate(self.Parent.Run.DUTs) if is_ref(dut)), default)

    return RefAnalysis
