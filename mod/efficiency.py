#!/usr/bin/env python
# --------------------------------------------------------
#       class for efficiency of a single DUT
# created on March 22nd 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from src.dut_analysis import DUTAnalysis, calc_eff, prep_kw, bins, get_2d_hist_vec, partial, mean, ax_range


def eff_analysis(cls):
    class Efficiency(cls):

        YArgs = {'y_tit': 'Efficiency [%]', 'y_range': [0, 105]}

        def __init__(self, parent: DUTAnalysis):  # noqa
            self.__dict__.update(parent.__dict__)
            self.MetaSubDir = 'Eff'
            self.Cut = self.Tracks.Cut

        def __call__(self):
            return self.value()

        def value(self, cut=None):
            return calc_eff(values=self.values(cut))

        def values(self, cut=None):
            return self.get_cluster_size(cut=cut).astype('?')

        def segment_values(self, nx=2, ny=3, cut=None):
            return get_2d_hist_vec(self.draw_map(local=True, cut=cut, binning=bins.make2d(*self.segments(nx, ny)), save=False), err=False, flat=False)

        # ----------------------------------------
        # region DRAW
        def draw(self, bw=None, **dkw):
            t, e = self.time(), self.values()
            return self.Draw.efficiency(t, e, **prep_kw(dkw, w=bw, **self.t_args(), **self.YArgs, stats=False, file_name='EfficiencyTrend'))

        def draw_vs_tp(self, **dkw):
            x, y = [f(cut=self.Cut.exclude('tp')) for f in [self.get_trigger_phase, self.values]]  # noqa
            return self.Draw.efficiency(x, y, bins.TP, 'Efficiency vs. Trigger Phase', **prep_kw(dkw, x_tit='Trigger Phase', **self.YArgs, file_name='EffieciencyTP'))

        def draw_map(self, res=.5, local=True, eff=True, both=False, fid=False, cut=None, **dkw):
            (x, y), e = [f(cut=self.Cut.get_nofid(cut, fid)) for f in [partial(self.Tracks.get_xy, local=local), self.values]]
            t = [*self.draw_text(self.Surface, cut, eff), *self.draw_text(not self.Surface, cut, eff and both)]
            b = bins.get_xy(local, self.Plane, res)
            return self.Draw.prof2d(x, y, e * 100, **prep_kw(dkw, title='Efficiency Map', binning=b, **self.ax_tits(local), leg=t, file_name='EfficiencyMap'))

        def draw_text(self, surface, cut, show=True):
            self.activate_surface(surface)
            x, y = [.5 if i is None else mean(ax_range(i)) for i in self.Cut.get_fid_config(surface)]
            return [self.Cut.get_fid(surface), self.Draw.textbox(f'{self.value(cut)[0]:2.1f}%', x, y, s=2, opacity=.4, ndc=False, show=False)] if show else []

        def draw_segments(self, res=.5, local=True, nx=10, ny=15, cut=None, **dkw):
            e = self.segment_values(nx, ny, cut)
            xs, ys = self.segment_centres(nx, ny)
            t = self.draw_grid(nx, ny) + [self.Draw.tlatex(x, y, '{:2.1f}'.format(e[i][j]), str(i * y.size + j), 22, size=.02) for j, x in enumerate(xs) for i, y in enumerate(ys)]
            return self.draw_map(res, local, cut=cut, eff=False, **prep_kw(dkw, leg=t, file_name=f'EffMap{nx}x{ny}'))

        def draw_segment_distribution(self, nx=10, ny=15, cut=None, segments=True, **dkw):
            e = self.segment_values(nx, ny, cut).flatten() if segments else get_2d_hist_vec(self.draw_map(.5, show=False), err=False)
            self.Draw.distribution(e, **prep_kw(dkw, title='Segment Efficiencies', x_tit='Efficiency [%]'))

        def draw_in_pixel(self, n=10, ox=0, oy=0, cut=None, **dkw):
            return super().draw_in_pixel(ox, oy, n, cut, self.values, tit='Efficiency', **dkw)

        def draw_in_cell(self, n=10, ox=0, oy=0, cut=None, **dkw):
            return super().draw_in_cell(ox, oy, n, cut, self.values, tit='Efficiency', **dkw)
        # endregion DRAW
        # ----------------------------------------

    return Efficiency
