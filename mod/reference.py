#!/usr/bin/env python
# --------------------------------------------------------
#       Analysis of the reference plane
# created on August 19th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from src.dut_analysis import DUTAnalysis
from src.utils import save_pickle
from plotting.draw import prep_kw, fit_fwhm, remove_low_stat_bins
from src.dut import Plane
import src.bins as bins


class RefAnalysis(DUTAnalysis):

    def __init__(self, parent):  # noqa

        self.Parent = parent
        self.__dict__.update(parent.__dict__)
        self.MetaSubDir = 'REF'

        self.Plane = Plane(self.Config.getint('DUT', 'reference plane'), self.Config, 'DUT')
        self.Data = self.F[str(self.Plane)]
        self.N = self.Data['Clusters']['X'].size

        from mod.ref_cuts import RefCut
        self.DUTCut = self.Cut
        self.Cut = RefCut(self)

    # ----------------------------------------
    # region GET
    def get_du(self, cut=False, **k):
        return (self.get_data('Clusters', 'U') - self.get_data('Tracks', 'U')[self.Cut['cluster']])[self.Cut(cut)]

    def get_dv(self, cut=False, **k):
        return (self.get_data('Clusters', 'V') - self.get_data('Tracks', 'V')[self.Cut['cluster']])[self.Cut(cut)]

    @save_pickle('RM')
    def get_res_means(self):
        return [fit_fwhm(f(show=False))[1] / 1e3 for f in [self.draw_x_residuals, self.draw_y_residuals]]
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw_correlation(self, m='y', cut=True, **dkw):
        c1, c2 = self.Cut.make_corrs(cut)
        v1, v2 = self.get_data('Clusters', m.upper()), self.Parent.get_data('Clusters', m.upper(), cut=False)
        return remove_low_stat_bins(self.Draw.histo_2d(v1[c1], v2[c2], bins.get_corr(m, self.Plane, self.Parent.Plane), **prep_kw(dkw, x_tit=f'{m} REF', y_tit=f'{m} DUT')), .1, of_max=True)
    # endregion DRAW
    # ----------------------------------------
