#!/usr/bin/env python
# --------------------------------------------------------
#       Analysis of the reference plane
# created on August 19th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from src.dut_analysis import DUTAnalysis
from mod.residuals import ResidualAnalysis
from mod.ref_cuts import RefCut


class RefAnalysis(DUTAnalysis):

    def __init__(self, parent):  # noqa

        self.Parent = parent
        self.__dict__.update(parent.__dict__)
        self.MetaSubDir = 'REF'

        self.Plane = self.Planes[self.Config.getint('DUT', 'reference plane')]
        self.N = self.n
        self.DUTCut = self.Cut
        self.Cut = RefCut(self)
        self.Residuals = ResidualAnalysis(self)
