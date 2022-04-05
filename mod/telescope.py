#!/usr/bin/env python
# --------------------------------------------------------
#       Telescope analysis of the high resolution tests at DESY and CERN
# created on June 11th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from src.dut_analysis import DUTAnalysis
from src.analysis import *
from mod.tel_cuts import TelCut
from mod.residuals import ResidualAnalysis


class TelescopeAnalysis(DUTAnalysis):

    def __init__(self, parent: DUTAnalysis):  # noqa

        self.Parent = parent
        self.__dict__.update(parent.__dict__)
        self.MetaSubDir = 'tel'
        
        self.NPlanes = self.Config.getint('TELESCOPE', 'planes')
        self.Plane = self.Planes[0]
        self.Cut = TelCut(self)
        self.Residuals = ResidualAnalysis(self)

    @property
    def planes(self):
        return arange(self.NPlanes)

    # ----------------------------------------
    # region DATA
    def get_all(self, f, cut=True, **kwargs):
        return [f(pl=pl, cut=self.Cut.make_all_cluster() & cut, **kwargs) for pl in range(self.NPlanes)]

    def get_n_clusters(self, cut=None, pl=0):
        return self.get_data('Clusters', 'N', cut, pl)

    def get_us(self, cut=True):
        return array(self.get_all(self.get_u, cut, centre=True))

    def get_vs(self, cut=True):
        return array(self.get_all(self.get_u, cut, centre=True))
    # endregion DATA
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw_occupancy(self, local=True, bw=10, cut=None, pl=None, **dkw):
        return super().draw_occupancy(local, bw, cut, fid=False, pl=pl, **dkw)

    def draw_n_clusters(self, pl=0, cut=None, **dkw):
        return self.Draw.distribution(self.get_n_clusters(cut, pl), **prep_kw(dkw, title='NClusters', x_tit='Number of Clusters', w=1, x0=-.5))
    # endregion DRAW
    # ----------------------------------------
