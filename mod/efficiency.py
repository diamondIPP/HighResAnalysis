#!/usr/bin/env python
# --------------------------------------------------------
#       class for efficiency of a single DUT
# created on March 22nd 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from src.dut_analysis import DUTAnalysis, calc_eff, prep_kw, bins, hist_values_2d, partial, mean, ax_range
from numpy import rad2deg


def eff_analysis(cls):
    class Efficiency(cls):

        YArgs = {'y_tit': 'Efficiency [%]', 'y_range': [0, 105]}
        ZArgs = {'z_tit': 'Efficiency [%]', 'z_range': [0, 100]}

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

        def pvalues(self, cut=None):
            return (self.get_cluster_size(cut=cut).astype('?') * 100).astype('d')

        def segment_values(self, nx=2, ny=3, cut=None):
            return hist_values_2d(self.draw_map(local=True, cut=cut, binning=bins.make2d(*self.segments(nx, ny)), save=False), err=False, flat=False)

        # ----------------------------------------
        # region DRAW
        def draw(self, bw=None, **dkw):
            t, e = self.time(), self.values()
            return self.Draw.efficiency(t, e, **prep_kw(dkw, w=bw, **self.t_args(), **self.YArgs, stats=False, file_name='EfficiencyTrend'))

        def draw_vs_tp(self, **dkw):
            x, y = [f(cut=self.Cut.exclude('tp')) for f in [self.get_trigger_phase, self.values]]  # noqa
            return self.Draw.efficiency(x, y, bins.TP, 'Efficiency vs. Trigger Phase', **prep_kw(dkw, x_tit='Trigger Phase', **self.YArgs, file_name='EffieciencyTP'))

        def draw_map(self, res=.5, local=True, eff=True, both=False, fid=False, cut=None, **dkw):
            (x, y), e = [f(cut=self.Cut.get_nofid(cut, fid)) for f in [partial(self.Tracks.get_xy, local=local), self.pvalues]]
            e[e < 1] = 1e-4  # zero values are not shown in the map ...
            t = [*self.draw_text(not self.Surface, cut, local, eff and both), *self.draw_text(self.Surface, cut, local, eff)]
            b = bins.get_xy(local, self.Plane, res)
            return self.Draw.prof2d(x, y, e, **prep_kw(dkw, title='Efficiency Map', binning=b, **self.ZArgs, **self.ax_tits(local), leg=t, file_name='EfficiencyMap'))

        def draw_text(self, surface, cut, local=True, show=True):
            if show:
                self.activate_surface(surface)
                x, y = [.5 if i is None else mean(ax_range(i)) for i in self.Cut.get_fid_config(surface, local=local)]
                return [self.Cut.get_fid(surface, local=local), self.Draw.textbox(f'{self.value(cut)[0]:2.1f}%', x, y, s=2 if local else .3, opacity=.4, ndc=False, show=False)]
            return []

        def draw_segments(self, res=.5, local=True, nx=10, ny=15, cut=None, **dkw):
            e = self.segment_values(nx, ny, cut)
            xs, ys = self.segment_centres(nx, ny)
            t = self.draw_grid(nx, ny) + [self.Draw.tlatex(x, y, '{:2.1f}'.format(e[i][j]), str(i * y.size + j), 22, size=.02) for j, x in enumerate(xs) for i, y in enumerate(ys)]
            return self.draw_map(res, local, cut=cut, eff=False, **prep_kw(dkw, leg=t, file_name=f'EffMap{nx}x{ny}'))

        def draw_segment_distribution(self, nx=10, ny=15, cut=None, segments=True, **dkw):
            e = self.segment_values(nx, ny, cut).flatten() if segments else hist_values_2d(self.draw_map(.5, show=False), err=False)
            self.Draw.distribution(e, **prep_kw(dkw, title='Segment Efficiencies', x_tit='Efficiency [%]'))

        def draw_in_pixel(self, n=10, ox=0, oy=0, cut=None, **dkw):
            return super().draw_in_pixel(ox, oy, n, cut, self.pvalues, tit='Efficiency', **prep_kw(dkw, **self.ZArgs))

        def draw_in_cell(self, n=10, ox=0, oy=0, cut=None, **dkw):
            return super().draw_in_cell(ox, oy, n, cut, self.pvalues, tit='Efficiency', **prep_kw(dkw, **self.ZArgs))

        def draw_col_vs_angle(self, r=10, bias=True, **dkw):
            cut = self.Cut.add(self.Cut.make_cell_fiducial(0, r, *self.DUT.PXYu / 2 if bias else (0, 0)))
            x, y, zz = rad2deg(self.get_slope_x(cut=cut)), rad2deg(self.get_slope_y(cut=cut)), self.pvalues(cut=cut)
            return self.Draw.prof2d(x, y, zz, **prep_kw(dkw, x_tit='Slope X [deg]', y_tit='Slope Y [deg]', **self.ZArgs, file_name='EffColvsAngle'))
        # endregion DRAW
        # ----------------------------------------

    return Efficiency
