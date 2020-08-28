#!/usr/bin/env python
# --------------------------------------------------------
#       small script to read simple text files written by pXar
# created on August 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from ROOT import TH1F
from numpy import cumsum, split, histogram, corrcoef, diff, all, invert, histogram2d, max, quantile, rad2deg, arctan

from analysis import *
from calibration import Calibration
from cern_run import CERNRun
from converter import Converter
from currents import Currents
from cut import Cut
from desy_converter import DESYConverter
from desy_run import DESYRun
from dummy import Dummy
from fit import *
from telescope import TelescopeAnalysis
from tracks import TrackAnalysis
from reference import RefAnalysis
from typing import Any


class DUTAnalysis(Analysis):

    def __init__(self, run_number, dut_number, test_campaign, single_mode=True, verbose=True, test=False):

        Analysis.__init__(self, test_campaign, verbose)
        self.print_start(run_number)

        # MAIN
        self.Run = self.init_run()(run_number, dut_number, self.TCDir, self.Config, single_mode)
        self.DUT = self.Run.DUT
        self.Plane = self.DUT.Plane
        self.MetaSubDir = 'DUT'

        # DATA
        self.Converter = self.init_converter()(self.TCDir, self.Run.Number, self.Config)
        self.Dummy = Dummy(self.TCDir, self.Converter.NTelPlanes, self.Converter.NDUTPlanes, self.Config)
        self.Data = self.load_data(test)
        self.init_cuts()

        # INFO
        self.NEvents = self.get_entries()
        self.StartTime = self.get_start_time()
        self.EndTime = self.get_end_time()
        self.Duration = (self.EndTime - self.StartTime).seconds
        self.Surface = False

        # SUBCLASSES
        self.Calibration = Calibration(self.Run)
        self.Telescope = TelescopeAnalysis(self)
        self.Tracks = TrackAnalysis(self)
        self.REF = RefAnalysis(self)
        self.Currents = Currents(self)
        self.add_cuts()
        self.add_track_cuts()

    # ----------------------------------------
    # region INIT
    def init_run(self):
        return DESYRun if self.Location == 'DESY' else CERNRun

    def init_converter(self):
        return DESYConverter if self.Location == 'DESY' else Converter

    def get_entries(self):
        return self.Data['Tracks']['NTracks'].size

    def load_data(self, test=False):
        if not test:
            self.Converter.run()
            try:
                f = self.load_file()
                _ = f['Tracks']
                _ = f[self.Plane.get_name()]
                return f
            except (KeyError, OSError):
                warning('could not load data file {} -> start with dummy'.format(self.Run.FileName))
        return self.Dummy.load_file()

    def load_file(self):
        return h5py.File(self.Run.FileName, 'r')

    def reload_data(self):
        self.Data = self.load_data()
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region CUTS
    def init_cuts(self, redo=False):
        self.Cuts.set_config(self.TestCampaign, self.DUT.Name)
        self.Cuts.register('fid', self.make_fiducial(redo=redo), 10, 'fid cut')
        self.Cuts.register('mask', self.make_mask(), 20, 'mask pixels')
        self.Cuts.register('charge', self.get_charges(cut=False) != 0, 30, 'events with non-zero charge')
        self.Cuts.register('cluster', self.get_data('Clusters', 'Size', cut=False) > 0, 90, 'tracks with a cluster')

    def add_cuts(self, redo=False):
        self.Cuts.register('res', self.REF.make_dut_residuals(redo), 69, 'small residuals to REF plane')
        self.Cuts.register('triggerphase', self.make_trigger_phase(redo=redo), 61, 'trigger phase')
        self.Cuts.register('tstart', self.make_start_time(redo=redo), 40, 'exclude first events')
        self.Cuts.register('chi2', self.make_chi2(redo=redo), 50, 'small chi2')

    def add_track_cuts(self, redo=False):
        self.Tracks.Cuts.register('triggerphase', self.make_trigger_phase(tracks=True, redo=redo), 10, 'track trigger phase')
        self.Tracks.Cuts.register('res', self.REF.make_residuals(redo=redo), 20, 'tracks with a small residual in the REF')
        self.Tracks.Cuts.register('fid', self.make_fiducial(tracks=True, redo=redo), 30, 'tracks in fiducial area')
        self.Tracks.Cuts.register('tstart', self.make_start_time(tracks=True, redo=redo), 40, 'exclude first events')
        self.Tracks.Cuts.register('chi2', self.make_chi2(tracks=True, redo=redo), 50, 'small chi2')

    def reload_cuts(self, redo=False):
        self.init_cuts(redo)
        self.add_cuts(redo)
        self.add_track_cuts(redo)

    def activate_surface(self):
        self.Cuts.register('fid', self.make_fiducial(option='surface fiducial'), 10)
        self.Tracks.Cuts.register('fid', self.make_fiducial(tracks=True, option='surface fiducial'), 30)
        self.Surface = True

    def deactivate_surface(self):
        self.Cuts.register('fid', self.make_fiducial(), 10, 'fid cut')
        self.Tracks.Cuts.register('fid', self.make_fiducial(tracks=True), 30, 'tracks in fiducial area')
        self.Surface = False

    def make_start_time(self, tracks=False, redo=False):
        def f():
            t0 = self.Cuts.get_config('start time', dtype=int) * 60
            t = self.get_time(cut=False, trk_cut=-1 if not tracks else False)
            return t > t[0] + t0
        return array(do_hdf5(self.make_hdf5_path('time', sub_dir='tracks' if tracks else None), f, redo))

    def make_res(self):
        x, y = self.get_du(cut=0), self.get_dv(cut=0)
        x0, y0 = self.get_du(), self.get_dv()
        mx, my = mean(x0[abs(x0) < .3]), mean(y0[abs(y0) < .3])
        return sqrt((x-mx) ** 2 + (y-my) ** 2) < .1

    def make_fiducial(self, tracks=False, redo=False, option='fiducial'):
        def f():
            x, y = self.Tracks.get_coods(local=True, trk_cut=False) if tracks else self.get_coods(local=True, cut=False)
            x0, x1, y0, y1 = self.Cuts.get_config(option, lst=True)
            return (x >= x0) & (x <= x1) & (y >= y0) & (y <= y1)
        return array(do_hdf5(self.make_hdf5_path(option[:3], sub_dir='tracks' if tracks else None), f, redo))

    def make_chi2(self, tracks=False, redo=False):
        def f():
            chi2 = self.Tracks.get_chi2(cut=False, trk_cut=-1)
            chi_max = quantile(chi2, self.Cuts.get_config('chi2 quantile', dtype=float))
            return self.Tracks.get_chi2(cut=False, trk_cut=False if tracks else -1) < chi_max
        return array(do_hdf5(self.make_hdf5_path('chi2', sub_dir='tracks' if tracks else None), f, redo))

    def make_mask(self):
        x, y = self.get_coods(local=True, cut=False)
        mx, my = array(self.Cuts.get_config('mask', lst=True)).T
        return all([invert((x >= mx[i] - .5) & (x <= mx[i] + .5) & (y >= my[i] - .5) & (y <= my[i] + .5)) for i in range(mx.size)], axis=0)

    def draw_fid_area(self, show=True, off=-.5):
        x1, x2, y1, y2 = self.Cuts.get_fid_config(self.Surface) + off
        self.draw_box(x1, y1, x2 + 1, y2 + 1, color=2, width=2, name='fid', show=show)

    def make_trigger_phase(self, tracks=False, redo=False):
        def f():
            tp = self.get_trigger_phase(cut=False, trk_cut=False if tracks else -1)
            low, high = self.Cuts.get_config('trigger phase', lst=True)
            return (tp >= low) & (tp <= high)
        return array(do_hdf5(self.make_hdf5_path('tp', sub_dir='tracks' if tracks else None), f, redo))

    def make_correlation(self, plane=2):
        n = self.Telescope.get_n('Clusters', plane, cut=False)
        return (n[self.Tracks.get_events()] == 1)[self.Cuts.get('cluster')()]

    def remove_metadata(self):
        if self.make_run_str():
            for f in glob(join(self.MetaDir, '*', '*_{}_{}*'.format(self.TestCampaign.strftime('%Y%m'), self.make_run_str()))):
                remove_file(f, join(basename(dirname(f)), basename(f)))
        self.Cuts.set_config(self.TestCampaign, self.DUT.Name)
        self.init_cuts()
        self.add_cuts()
        self.add_track_cuts()
    # endregion CUTS
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_start_time(self):
        return datetime.fromtimestamp(self.Run.StartTime)

    def get_end_time(self):
        return datetime.fromtimestamp(self.Run.EndTime)

    def get_time(self, cut=None, trk_cut: Any = -1):
        t = array(self.Data['Event']['Time']).astype('f8') + self.Run.StartTime
        t = t[self.Tracks.get_events()]  # map the event times to the track times
        return t[self.Cuts.get('cluster')()][self.Cuts(cut)] if trk_cut == -1 else t[self.Tracks.Cuts(trk_cut)]

    def get(self):
        return self.Data['Plane{}'.format(self.Plane.Number)]

    def get_group(self, grp):
        return self.get()[grp]

    def get_data(self, grp, key=None, cut=None):
        data = self.get_group(grp)
        data = array(data) if key is None else array(data[key])
        return data[self.Cuts(cut)]

    def get_charges(self, e=False, cut=None):
        return self.get_data('Clusters', 'Charge', cut) * (self.DUT.VcalToEl if e else 1)

    def get_charge(self, cut=None):
        values = self.get_charges(cut=cut)
        m, s = mean_sigma(values)
        return ufloat(m, s / sqrt(values.size))

    def get_x(self, cut=None):
        return self.get_data('Clusters', 'X', cut)

    def get_y(self, cut=None):
        return self.get_data('Clusters', 'Y', cut)

    def get_u(self, cut=None):
        return self.get_data('Clusters', 'U', cut=cut)

    def get_v(self, cut=None):
        return self.get_data('Clusters', 'V', cut=cut)

    def get_du(self, cut=None):
        return self.get_u(cut) - self.Tracks.get_u(cut)

    def get_dx(self, cut=None):
        return self.get_x(cut) - self.Tracks.get_x(cut)

    def get_dv(self, cut=None):
        return self.get_v(cut=cut) - self.Tracks.get_v(cut)

    def get_dy(self, cut=None):
        return self.get_y(cut) - self.Tracks.get_y(cut)

    def get_residuals(self, cut=None):
        return sqrt(self.get_du(cut) ** 2 + self.get_dv(cut) ** 2)

    def get_coods(self, local=True, cut=None):
        return array([self.get_x(cut=cut) if local else self.get_u(cut), self.get_y(cut=cut) if local else self.get_v(cut)])

    def get_mask(self):
        return self.get_data('Mask', cut=False)

    def get_time_args(self, rel_t=False):
        return {'x_tit': 'Time [hh:mm]', 't_ax_off': self.Run.StartTime if rel_t else 0}

    @staticmethod
    def get_ax_tits(local):
        return {'x_tit': 'Column' if local else 'X [mm]', 'y_tit': 'Row' if local else 'Y [mm]'}

    def get_track_data(self, grp, key=None, cut=None, trk_cut=-1):
        data = self.get_data(grp, key, cut=False)
        return data[self.Tracks.Cuts(trk_cut)] if type(trk_cut) == ndarray or trk_cut != -1 else data[self.Cuts.get('cluster').Values][self.Cuts(cut)]

    def get_cluster_size(self, cut=None, trk_cut: Any = -1):
        return self.get_track_data('Clusters', 'Size', cut=cut, trk_cut=trk_cut)

    def get_efficiencies(self, trk_cut=None):
        return (self.get_cluster_size(trk_cut=trk_cut) > 0).astype('u2') * 100

    def get_segment_efficiecies(self, nx=2, ny=3, cut=None):
        return get_2d_hist_vec(self.draw_efficiency_map(local=True, cut=cut, binning=bins.make2d(*self.get_segments(nx, ny)), show=False), err=False, flat=False)

    def get_efficiency(self, trk_cut=None, prnt=True):
        eff = calc_eff(values=self.get_efficiencies(trk_cut))
        self.info('{:.1f}(-{:.1f}+{:.1f})%'.format(*eff), prnt=prnt)
        return eff

    def get_trigger_phase(self, cut=None, trk_cut: Any = -1):
        return self.get_track_data('TriggerPhase', cut=cut, trk_cut=trk_cut)

    def get_events(self):
        """ returns: event numbers with clusters. """
        return self.Tracks.get_events(self.Cuts.get('cluster'))

    def get_segments(self, nx, ny):
        x0, x1, y0, y1 = self.Cuts.get_config('full size', lst=True)
        return arange(x0, x1 + bool((x1 - x0 + 1) % nx) * nx + 1.1, nx, dtype='u2'), arange(y0, y1 + bool((y1 - y0 + 1) % ny) * ny + 1.1, ny, dtype='u2')

    def expand_inpixel(self, x, y, e=None, cell=False):
        cx, cy = [self.DUT.CellSize / 1000. / self.Plane.PX, self.DUT.CellSize / 1000. / self.Plane.PY]
        x, y = [(x % cx) / cx,  (y % cy) / cy] if cell else [x % 1, y % 1]
        # add edges of the neighbouring pixels
        xf = concatenate([x, x[x < .5] + 1, x[x >= .5] - 1, x[y < .5], x[y > .5]])
        yf = concatenate([y, y[x < .5], y[x >= .5], y[y < .5] + 1, y[y >= .5] - 1])
        ef = concatenate([e, e[x < .5], e[x >= .5], e[y < .5], e[y >= .5]]) if e is not None else None
        # add corners
        xf = concatenate([xf, x[(x < .5) & (y < .5)] + 1, x[(x >= .5) & (y >= .5)] - 1, x[(x < .5) & (y >= .5)] + 1, x[(x >= .5) & (y < .5)] - 1])
        yf = concatenate([yf, y[(x < .5) & (y < .5)] + 1, y[(x >= .5) & (y >= .5)] - 1, y[(x < .5) & (y >= .5)] - 1, y[(x >= .5) & (y < .5)] + 1])
        ef = concatenate([ef, e[(x < .5) & (y < .5)], e[(x >= .5) & (y >= .5)], e[(x < .5) & (y >= .5)], e[(x >= .5) & (y < .5)]]) if e is not None else None
        return [xf, yf] if e is None else [xf, yf, ef]
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw_mask(self, show=True):
        self.format_statbox(entries=True)
        x, y = self.get_mask().T
        title = 'Masked Pixels in {}'.format(self.DUT.Name)
        self.draw_histo_2d(x, y, bins.get_local(self.Plane), title, x_tit='Column', y_tit='Row', fill_color=1, draw_opt='box', rm=.03, show=show)

    def draw_occupancy(self, local=True, bin_width=1, cut=None, fid=False, show=True):
        self.format_statbox(entries=True, x=.83, m=True)
        cut = self.Cuts(cut) if fid else self.Cuts.exclude('fid', cut)
        x, y = self.get_coods(local, cut)
        title = '{} Cluster Occupancy'.format('Local' if local else 'Global')
        self.draw_histo_2d(x, y, bins.get_coods(local, self.Plane, bin_width), title, x_tit='Column', y_tit='Row', show=show)

    def draw_hit_map(self, res=.3, local=True, cut=None, fid=False, show=True):
        self.format_statbox(entries=True, x=.84)
        cut = self.Cuts(cut) if fid else self.Cuts.exclude('fid', cut)
        x, y = self.Tracks.get_coods(local, cut)
        self.draw_histo_2d(x, y, bins.get_coods(local, self.Plane, res), 'Hit Map', show=show, **self.get_ax_tits(local))

    def draw_cluster_size(self, show=True, cut=None, trk_cut=-1):
        self.format_statbox(all_stat=True)
        v = self.get_cluster_size(cut, trk_cut)
        self.draw_disto(v, bins.make(0, 10), 'Cluster Size in {}'.format(self.Plane), show=show, x_tit='Cluster Size', lm=.14, y_off=2)

    def draw_cluster_size_map(self, res=.3, local=True, cut=None, fid=False, show=True):
        self.format_statbox(entries=True, x=.84)
        cut = self.Cuts(cut) if fid else self.Cuts.exclude('fid')
        x, y = self.Tracks.get_coods(local, cut)
        self.draw_prof2d(x, y, self.get_cluster_size(cut), bins.get_coods(local, self.Plane, res), 'Cluster Size', show=show, z_tit='Cluster Size', **self.get_ax_tits(local))

    def draw_trigger_phase(self, cut=None, trk_cut=-1):
        self.format_statbox(entries=True)
        cut, trk_cut = self.Cuts.exclude('triggerphase', cut), self.Tracks.Cuts.exclude('triggerphase', trk_cut) if trk_cut != -1 else trk_cut
        h = self.draw_disto(self.get_trigger_phase(cut, trk_cut), bins.get_triggerphase(), 'Trigger Phase', x_tit='Trigger Phase', y_off=1.8, lm=.13)
        format_histo(h, y_range=[0, h.GetMaximum() * 1.1])
        update_canvas()

    def draw_time(self, show=True):
        t = self.get_time()
        g = self.make_graph_from_profile(self.draw_prof(arange(t.size), t, bins.make(0, t.size, sqrt(t.size)), show=False), err=0)
        format_histo(g, x_tit='Event Number', y_tit='Time [hh:mm]', y_off=1.8)
        set_time_axis(g, axis='Y')
        self.draw_histo(g, show, .13, draw_opt='al')

    def draw_segments(self, nx=2, ny=3, w=1):
        x, y = array(self.get_segments(nx, ny)) - .5
        self.draw_grid(x, y, w)

    def draw_x(self, res=.1, cluster=False, cut=None):
        self.format_statbox(all_stat=True)
        x = self.get_x(cut) if cluster else self.Tracks.get_x(cut)
        self.draw_disto(x, bins.make(0, self.Plane.NCols, res, last=True))

    def draw_inpixel_map(self, res=.1, cut=None, cell=False, show=True):
        x, y = self.expand_inpixel(cell=cell, *self.Tracks.get_coods(cut=cut))
        self.draw_histo_2d(x, y, bins.get_pixel(self.Plane, res, cell=cell), 'Hit Map in {}'.format('3D Cell' if cell else 'Pixel'), show=show, stats=0)
        self.draw_box(0, 0, 1, 1)
        update_canvas()
    # endregion DRAW
    # ----------------------------------------

    # ----------------------------------------
    # region RESIDUALS
    def draw_x_residuals(self, cut=None):
        self.format_statbox(all_stat=True)
        self.draw_disto(self.get_du(cut), bins.make(-3, 3, .01), 'X Residuals', x_tit='Residual [mm]', lm=.12, y_off=1.8)

    def draw_udv(self, cut=None, tv=None, show=True):
        self.format_statbox(entries=True, x=.84)
        u, dv = self.get_u(cut), (self.get_dv(cut) if tv is None else self.get_v() - tv)
        h = self.draw_histo_2d(u, dv, bins.get_global_x(self.Plane) + bins.make(-1, 1, .02), 'X dY', x_tit='X [col]', y_tit='dY [row]', show=show)
        # h = self.draw_prof(u, dv, bins.get_global_x(self.Plane), 'X dY', x_tit='X [col]', y_tit='dY [row]', show=show)
        fit = h.Fit('pol1', 'qs')
        return arctan(fit.Parameter(1))

    def draw_vdu(self, cut=None, tu=None, show=True):
        self.format_statbox(entries=True, x=.84)
        v, du = self.get_v(cut), (self.get_du(cut) if tu is None else self.get_u() - tu)
        h = self.draw_histo_2d(v, du, bins.get_global_y(self.Plane) + bins.make(-1, 1, .02), 'X dY', x_tit='X [col]', y_tit='dY [row]', show=show)
        # h = self.draw_prof(v, du, bins.get_global_y(self.Plane), 'X dY', x_tit='X [col]', y_tit='dY [row]', show=show)
        fit = h.Fit('pol1', 'qs')
        return arctan(fit.Parameter(1))

    def draw_xdy(self, cut=None, ty=None):
        self.format_statbox(entries=True, x=.84)
        x, dy = self.get_x(cut), (self.get_dy(cut) if ty is None else self.get_y() - ty)
        h = self.draw_histo_2d(x, dy, bins.get_local_x(self.Plane) + bins.make(-3, 3, .1), 'X dY', x_tit='X [col]', y_tit='dY [row]')
        fit = h.Fit('pol1', 'qs')
        return arctan(fit.Parameter(1))

    def draw_ydx(self, cut=None, tx=None):
        self.format_statbox(entries=True, x=.84)
        y, dx = self.get_y(cut), (self.get_dx(cut) if tx is None else self.get_x() - tx)
        h = self.draw_histo_2d(y, dx, bins.get_local_y(self.Plane) + bins.make(-3, 3, .1), 'Y dX', x_tit='Y [row]', y_tit='dX [col]')
        fit = h.Fit('pol1', 'qs')
        return rad2deg(arctan(fit.Parameter(1)))

    def draw_y_residuals(self, cut=None):
        self.format_statbox(all_stat=True)
        self.draw_disto(self.get_dv(cut), bins.make(-3, 3, .01), 'Y Residuals', x_tit='Residual [mm]', lm=.12, y_off=1.8)

    def draw_residuals(self, show=True):
        self.draw_disto(self.get_residuals(), bins.make(0, 6, .01), 'Residuals', x_tit='Residual [mm]', show=show)

    def draw_residuals_map(self, res=.3, local=True, cut=None, fid=False, show=True):
        self.format_statbox(entries=True, x=.84)
        cut = self.Cuts(cut) if fid else self.Cuts.exclude('fid')
        x, y = self.Tracks.get_coods(local, cut)
        self.draw_prof2d(x, y, self.get_residuals(cut), bins.get_coods(local, self.Plane, res), 'Residuals', z_tit='Residuals', show=show, **self.get_ax_tits(local))
    # endregion RESIDUALS
    # ----------------------------------------

    # ----------------------------------------
    # region CORRELATION
    def draw_correlation(self, mode='y', res=1, plane=2, thresh=.1, show=True):
        return self.Telescope.draw_correlation(mode, res, plane, thresh, show)

    def draw_x_correlation(self, res=1, plane=2, show=True):
        self.draw_correlation('x', res, plane, show)

    def draw_y_correlation(self, res=1, plane=2, show=True):
        self.draw_correlation('y', res, plane, show)

    def draw_correlation_trend(self, plane=2, bin_width=120, thresh=.5, show=True):
        v0, v1 = self.get_y(cut=self.make_correlation(plane)), self.Telescope.get_y(plane, cut=self.Telescope.make_correlation(plane))
        splits, t = histogram(self.get_time(self.make_correlation(plane)), bins.get_time(self.get_time(), min(int(self.Duration / 3), bin_width))[1])
        hs = [histogram2d(v[0], v[1], [bins.get_local_y(self.Plane, 2)[1], bins.get_local_y(self.Telescope.Plane, 4)[1]])[0] for v in split(array([v0, v1]), cumsum(splits)[:-1], axis=1)]
        # take only bins with >25% of max entries (noise filter)
        values = [corrcoef(array(where(h > max(h) * thresh)).repeat(h[h > max(h) * thresh].astype('i'), axis=1).T, rowvar=False)[0][1] for h in hs]
        g = self.make_tgrapherrors('gal', 'Event Correlation', x=t[:-1] + diff(t) / 2, y=values)
        format_histo(g, y_tit='Correlation Factor', y_off=1.6, **self.get_time_args(), y_range=[0, 1.05])
        self.draw_histo(g, show, .13, draw_opt='ap')
        return g

    def draw_alignment(self, plane=2, bin_width=120, thresh=.5, show=True):
        x, c = get_graph_vecs(self.draw_correlation_trend(plane, bin_width, show=False), err=False)
        r = [1 if ic < thresh else 2 for ic in c]
        x, y = x.repeat(r), ones(sum(r))
        binning = bins.get_time(self.get_time(), min(int(self.Duration / 3), bin_width)) + [3, 0, 3]
        gStyle.SetPalette(3, array([1, 2, 3], 'i'))
        self.draw_histo_2d(x, y, binning, 'Event Alignment', show=show, t_ax_off=0, x_tit='Time [hh:mm]', y_tit='Alignment', stats=False, l_off_y=99, center_y=True, draw_opt='col', rm=.03)
    # endregion CORRELATION
    # ----------------------------------------

    # ----------------------------------------
    # region SIGNAL
    def draw_ped_map(self, c_max, c_min=None):
        charge = self.get_charges(cut=False)
        cut = Cut('charge', charge < c_max if c_min is None else (charge > c_min) & (charge < c_max)) + self.Cuts.exclude('fid')
        self.draw_hit_map(cut=cut)

    def draw_charge_map(self, res=.3, cluster=False, fid=False, cut=None, show=True):
        cut = self.Cuts(cut) if fid else self.Cuts.exclude('fid', cut)
        res = 1 if cluster else res
        x, y = (self.get_x(cut), self.get_y(cut)) if cluster else (self.Tracks.get_x(cut), self.Tracks.get_y(cut))
        self.format_statbox(entries=True, x=.84)
        self.draw_prof2d(x, y, self.get_charges(cut=cut), bins.get_local(self.Plane, res), 'Charge Map', x_tit='Column', y_tit='Row', z_tit='Charge [vcal]', show=show)
        self.draw_fid_area(not fid and show)

    def draw_charge_occupancy(self, fid=False, cut=None, show=True):
        self.draw_charge_map(1, cluster=True, fid=fid, cut=cut, show=show)

    def draw_charge_distribution(self, bin_width=4, cut=None, x_range=None, show=True):
        self.format_statbox(all_stat=True)
        self.draw_disto(self.get_charges(cut=cut), bins.get_vcal(bin_width), 'Cluster Charge', x_tit='Charge [vcal]', x_range=x_range, show=show, y_off=1.8, lm=.12)
        self.draw_vertical_line(self.Calibration.Trim, 0, 1e5, w=2)
        self.draw_tlatex(.2, .5, 'Threshold = {} vcal'.format(self.Calibration.Trim), ndc=True, angle=90, size=.04)

    def draw_signal_distribution(self, bin_width=200, x_range=None, cut=None, thresh=False, show=True):
        self.format_statbox(all_stat=True)
        values = self.get_charges(cut=cut) * self.DUT.VcalToEl
        h = self.draw_disto(values, bins.get_electrons(bin_width), 'Pulse Height', x_tit='Pulse Height [e]', x_range=x_range, show=show)
        if thresh:
            self.draw_vertical_line(self.Calibration.Trim * self.DUT.VcalToEl, 0, h.GetMaximum(), w=2)
            self.draw_tlatex(.14, .5, 'Pixel Threshold #approx {:1.0f} e'.format(round(self.Calibration.Trim * self.DUT.VcalToEl, -2)), ndc=True, angle=90, size=.04, align=22)
        self.draw_preliminary()
        update_canvas()
        return h

    def draw_charge_vs_trigger_phase(self, cut=None, show=True):
        self.format_statbox(entries=True)
        cut = self.Cuts.exclude('triggerphase', cut)
        x, y = self.get_trigger_phase(cut=cut), self.get_charges(cut=cut)
        self.draw_prof(x, y, bins.get_triggerphase(), 'Charge vs. Trigger Phase', show=show, x_tit='Trigger Phase', y_tit='Charge [vcal]')

    def draw_charge_trend(self, bin_width=30, e=False, y_range=None, cut=None, show=True, stats=True):
        self.format_statbox(entries=True, exe=stats)
        t, charges = self.get_time(cut), self.get_charges(e=e, cut=cut)
        p = self.draw_prof(t, charges, bins.get_time(t, bin_width), x_tit='Time [hh:mm}', y_tit='Charge [{}]'.format('e' if e else 'vcal'), show=show, t_ax_off=0)
        values = get_hist_vec(p)
        # noinspection PyUnresolvedReferences
        format_histo(p, y_range=choose(y_range, increased_range([min(values).n, max(values).n], 1, 2)))
        update_canvas()
        return p

    def fit_charge(self, bin_width=30, y_range=None, show=False):
        self.format_statbox(only_fit=True)
        h = self.draw_charge_trend(bin_width=bin_width, y_range=y_range, show=show, stats=False)
        format_histo(h, stats=1, name='Fit Result')
        fit = h.Fit('pol0', 'sq')
        return fit

    def draw_inpixel_charge(self, res=.1, cut=None, show=True, cell=False):
        (x, y), c = self.Tracks.get_coods(cut), self.get_charges(cut=cut)
        x, y, c = self.expand_inpixel(x, y, c, cell)
        self.draw_prof2d(x, y, c, bins.get_pixel(self.Plane, res, cell=cell), 'Charge Map in {}'.format('3D Cell' if cell else 'Pixel'), show=show, stats=0)
        self.draw_box(0, 0, 1, 1)
        update_canvas()
    # endregion SIGNAL
    # ----------------------------------------

    # ----------------------------------------
    # region EFFICIENCY
    def draw_efficiency(self, bin_width=30, show=True):
        t, e = self.get_time(trk_cut=None), self.get_efficiencies()
        p = self.draw_eff(t, e, bins.get_time(t, bin_width), 'Efficiency', x_tit='Time [hh:mm]', y_tit='Efficiency [%]', t_ax_off=0, y_range=[0, 105], show=show, stats=0)
        return p

    def draw_efficiency_vs_trigger_phase(self, show=True):
        self.format_statbox(entries=True, y=.2)
        cut = self.Tracks.Cuts.exclude('triggerphase')
        x, y = self.get_trigger_phase(trk_cut=cut), self.get_efficiencies(cut)
        return self.draw_prof(x, y, bins.get_triggerphase(), 'Efficiency vs. Trigger Phase', x_tit='Trigger Phase', y_tit='Efficiency [%]', y_range=[0, 105], show=show)

    def draw_efficiency_map(self, res=.25, local=True, eff=True, fid=False, cut=None, binning=None, show=True):
        mcut = self.Tracks.Cuts(cut) if fid else self.Tracks.Cuts.exclude('fid')
        x, y = self.Tracks.get_x(trk_cut=mcut, local=local), self.Tracks.get_y(trk_cut=mcut, local=local)
        self.format_statbox(entries=True, x=.84)
        binning = choose(binning, bins.get_coods, 'None', local, self.Plane, res)
        p = self.draw_prof2d(x, y, self.get_efficiencies(mcut), binning, 'Efficiency Map', show=show, draw_opt='colz', **self.get_ax_tits(local))
        self.draw_fid_area(show=not fid and show)
        if eff:
            x0, x1, y0, y1 = self.Cuts.get_fid_config(self.Surface)
            self.draw_tlatex(x0 + (x1 - x0) / 2, y0 + (y1 - y0) / 2, '{:2.1f}%'.format(self.get_efficiency(cut, False)[0]), 'eff', 22, size=.04)
        return p

    def draw_segment_efficiencies(self, res=.25, local=True, nx=2, ny=3, cut=None, show=True):
        e = self.get_segment_efficiecies(nx, ny, cut)
        self.draw_efficiency_map(res, local, cut=cut, show=show, eff=False)
        x, y = self.get_segments(nx, ny)
        for i, ix in enumerate(x[:-1]):
            for j, iy in enumerate(y[:-1]):
                self.draw_tlatex(ix + nx / 2 - .5, iy + ny / 2 - .5, '{:2.1f}'.format(e[i][j]), str(i * y.size + j), 22, size=.02)
        self.draw_segments(nx, ny)

    def draw_segment_distribution(self, nx=2, ny=3, full=False, cut=None, bin_width=None, segments=True, show=True):
        e = self.get_segment_efficiecies(nx, ny, cut).flatten() if segments else get_2d_hist_vec(self.draw_efficiency_map(.5, show=False), err=False)
        self.format_statbox(all_stat=True)
        binning = bins.make(0, 101.5, choose(bin_width, 1)) if full else bins.make(95, 100.5, choose(bin_width, 5.5 / sqrt(e.size) / 2), last=True)
        self.draw_disto(e, binning, 'Segment Efficiencies', x_tit='Efficiency [%]', show=show)

    def draw_inpixel_eff(self, res=.1, cut=None, show=True, cell=False):
        (x, y), e = self.Tracks.get_coods(trk_cut=cut), self.get_efficiencies(cut)
        x, y, e = self.expand_inpixel(x, y, e, cell)
        self.draw_prof2d(x, y, e, bins.get_pixel(self.Plane, res, cell=cell), 'Efficiency Map in {}'.format('3D Cell' if cell else 'Pixel'), show=show, stats=0)
        self.draw_box(0, 0, 1, 1)
        update_canvas()
    # endregion EFFICIENCY
    # ----------------------------------------

    def fit_langau(self, h=None, nconv=30, show=True, chi_thresh=8, fit_range=None):
        h = self.draw_signal_distribution(show=show) if h is None and hasattr(self, 'draw_signal_distribution') else h
        fit = Langau(h, nconv, fit_range)
        fit.get_parameters()
        fit(show=show)
        get_last_canvas().Modified()
        get_last_canvas().Update()
        if fit.get_chi2() > chi_thresh and nconv < 80:
            self.Count += 5
            self.info('Chi2 too large ({c:2.2f}) -> increasing number of convolutions by 5'.format(c=fit.get_chi2()))
            fit = self.fit_langau(h, nconv + self.Count, chi_thresh=chi_thresh, show=show)
        print('MPV: {:1.1f}'.format(fit.get_mpv()))
        self.Count = 0
        self.add(fit)
        return fit

    def estimate_mean(self, n=10, fit_range=None):
        fit_range = [8000, 30000] if fit_range is None else fit_range
        self.fit_langau(show=False, fit_range=fit_range)
        values = []
        self.PBar.start(n)
        for i in range(n):
            h = TH1F('hph{}'.format(i), 'Pulse Height Distribution', 200, -5000, 70000)
            h.FillRandom('Fitfcn_hph', 25000)
            values.append(h.GetMean())
            self.PBar.update(i)
        return mean_sigma(values)


if __name__ == '__main__':
    t_start = time()
    aparser = ArgumentParser()
    aparser.add_argument('run', nargs='?', default=11, type=int)
    aparser.add_argument('dut', nargs='?', default=1, type=int)
    aparser.add_argument('--testcampaign', '-tc', nargs='?', default=None)
    aparser.add_argument('--verbose', '-v', action='store_true')
    aparser.add_argument('--single_mode', '-s', action='store_false')
    aparser.add_argument('--test', '-t', action='store_true')
    args = aparser.parse_args()
    z = DUTAnalysis(args.run, args.dut, test_campaign=args.testcampaign, single_mode=args.single_mode, verbose=args.verbose, test=args.test)
    z.add_info(t_start, prnt=True)
