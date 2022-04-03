#!/usr/bin/env python
# --------------------------------------------------------
#       cuts for the telescope planes
# created on March 30th 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from numpy import ndarray, all
from mod.dut_cuts import DUTCut


class TelCut(DUTCut):
    def __init__(self, ana):
        super().__init__(ana, meta_sub_dir='tel_cuts')

    def __call__(self, cut=None, data=None, pl=0):
        cut = super().__call__(cut, pl=pl)
        if data is None:
            return cut if type(cut) is ndarray and cut.size != self.Ana.NTracks else self.trk2pl(cut, pl)
        return data if cut is ... else data[cut] if cut.size == data.size else data[self.trk2pl(cut, pl)]

    def make(self, redo=False):
        pass

    # ----------------------------------------
    # region GENERATE
    def make_cluster(self, pl=0):
        return self.Ana.get_data('Clusters', 'Size', cut=False, pl=pl) > 0

    def make_all_cluster(self):
        return all([self.make_cluster(pl) for pl in self.Ana.planes], axis=0)

    def make_correlation(self, pl0, pl1):
        return self.make_cluster(pl0) & self.make_cluster(pl1)

    def trk2pl(self, x, pl=0):
        return x[self.make_cluster(pl)]
    # endregion GENERATE
    # ----------------------------------------
