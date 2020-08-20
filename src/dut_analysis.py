#!/usr/bin/env python
# --------------------------------------------------------
#       small script to read simple text files written by pXar
# created on August 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from ROOT import TH1F
from numpy import cumsum, split, histogram, corrcoef, diff, all, invert

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

        # DATA
        self.Converter = self.init_converter()(self.TCDir, self.Run.Number, self.Config)
        self.Dummy = Dummy(self.TCDir, self.Converter.NTelPlanes, self.Converter.NDUTPlanes, self.Config)
        self.Data = self.load_file(test)
        self.init_cuts()

        # SUBCLASSES
        self.Calibration = Calibration(self.Run)
        self.Telescope = TelescopeAnalysis(self)
        self.Tracks = TrackAnalysis(self)
        self.REF = RefAnalysis(self)
        self.Currents = Currents(self)
        self.add_cuts()
        self.add_track_cuts()

        # INFO
        self.NEvents = self.get_entries()
        self.StartTime = self.get_start_time()
        self.EndTime = self.get_end_time()

    # ----------------------------------------
    # region INIT
    def init_run(self):
        return DESYRun if self.Location == 'DESY' else CERNRun

    def init_converter(self):
        return DESYConverter if self.Location == 'DESY' else Converter

    def get_entries(self):
        return self.Data['Tracks']['NTracks'].size

    def load_file(self, test=False):
        if not test:
            self.Converter.run()
            try:
                f = h5py.File(self.Run.FileName, 'r')
                _ = f['Tracks']
                _ = f[self.Plane.get_name()]
                return f
            except (KeyError, OSError):
                warning('could not load data file {} -> start with dummy'.format(self.Run.FileName))
        return self.Dummy.load_file()

    def reload_file(self):
        self.Data = self.load_file()
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region CUTS
    def init_cuts(self):
        self.Cuts.register('fid', self.make_fiducial(), 10, 'fid cut')
        self.Cuts.register('mask', self.make_mask(), 20, 'mask pixels')
        self.Cuts.register('charge', self.get_charges(cut=False) != 0, 60, 'events with non-zero charge')
        self.Cuts.register('cluster', self.get_data('Clusters', 'Size', cut=False) > 0, 90, 'tracks with a cluster')

    def make_fiducial(self, tracks=False):
        x, y = self.Tracks.get_coods(local=True, trk_cut=False) if tracks else self.get_coods(local=True, cut=False)
        x0, x1, y0, y1 = self.Cuts.get_config('fiducial', lst=True)
        return (x >= x0) & (x <= x1) & (y >= y0) & (y <= y1)

    def make_mask(self):
        x, y = self.get_coods(local=True, cut=False)
        mx, my = array(self.Cuts.get_config('mask', lst=True)).T
        return all([invert((x > mx[i] - .5) & (x < mx[i] + .5) & (y > my[i] - .5) & (y < my[i] + .5)) for i in range(mx.size)], axis=0)

    def draw_fid_area(self):
        x1, x2, y1, y2 = self.Cuts.get_config('fiducial', lst=True)
        self.draw_box(x1, y1, x2, y2, color=2, width=2, name='fid')

    def make_trigger_phase(self, tracks=False):
        tp = self.get_trigger_phase(cut=False, trk_cut=False if tracks else -1)
        low, high = self.Cuts.get_config('trigger phase', lst=True)
        return (tp >= low) & (tp <= high)

    def add_cuts(self):
        self.Cuts.register('res', self.REF.make_dut_residuals(), 69, 'small residuals to REF plane')
        self.Cuts.register('triggerphase', self.make_trigger_phase(), 61, 'trigger phase')

    def add_track_cuts(self):
        self.Tracks.Cuts.register('triggerphase', self.make_trigger_phase(tracks=True), 10, 'track trigger phase')
        self.Tracks.Cuts.register('res', self.REF.make_residuals(), 20, 'tracks with a small residual in the REF')
        self.Tracks.Cuts.register('fid', self.make_fiducial(tracks=True), 30, 'tracks in fiducial area')
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
        return self.get_u(cut=cut) - self.Tracks.get_u(cut)

    def get_dv(self, cut=None):
        return self.get_v(cut=cut) - self.Tracks.get_v(cut)

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

    def get_efficiency(self, trk_cut=None):
        return (self.get_cluster_size(trk_cut=trk_cut) > 0).astype('u2') * 100

    def get_trigger_phase(self, cut=None, trk_cut: Any = -1):
        return self.get_track_data('TriggerPhase', cut=cut, trk_cut=trk_cut)
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw_mask(self, show=True):
        self.format_statbox(entries=True)
        x, y = self.get_mask().T
        title = 'Masked Pixels in {}'.format(self.DUT.Name)
        self.draw_histo_2d(x, y, bins.get_local(self.Plane), title, x_tit='Column', y_tit='Row', fill_color=1, draw_opt='box', rm=.03, show=show)

    def draw_occupancy(self, local=True, bin_width=1, cut=None, show=True):
        self.format_statbox(entries=True, x=.83, m=True)
        x, y = self.get_coods(local, cut)
        title = '{} Cluster Occupancy'.format('Local' if local else 'Global')
        self.draw_histo_2d(x, y, title, bins.get_coods(local, self.Plane, bin_width), x_tit='Column', y_tit='Row', show=show)

    def draw_x_residuals(self, cut=None):
        self.format_statbox(all_stat=True)
        self.draw_disto(self.get_du(cut), bins.make(-3, 3, .01), 'X Residuals', x_tit='Residual [mm]', lm=.12, y_off=1.8)

    def draw_y_residuals(self, cut=None):
        self.format_statbox(all_stat=True)
        self.draw_disto(self.get_dv(cut), bins.make(-3, 3, .01), 'Y Residuals', x_tit='Residual [mm]', lm=.12, y_off=1.8)

    def draw_residuals(self):
        self.draw_disto(self.get_residuals(), bins.make(0, 6, .01), 'Residuals', x_tit='Residual [mm]')

    def draw_correlation(self, mode, plane=2, show=True):
        tel_plane = self.Telescope.Plane(plane)
        evnt_cut, dut_cut, tel_cut = self.get_corr_cuts(plane)
        x, y = (self.get_x(cut=dut_cut), self.get_x(tel_plane, cut=tel_cut)) if mode.lower() == 'x' else (self.get_y(cut=dut_cut), self.get_y(tel_plane, cut=tel_cut))
        self.format_statbox(entries=True, x=.83)
        self.draw_histo_2d(x, y, 'Cluster Correlation in {} with Plane {}'.format(mode.upper(), plane), bins.get_corr(mode, self.Plane, tel_plane), x_tit='Col', y_tit='Row', show=show)

    def draw_x_correlation(self, plane=2, show=True):
        self.draw_correlation('x', plane, show)

    def draw_y_correlation(self, plane=2, show=True):
        self.draw_correlation('y', plane, show)

    def draw_alignment(self, mode='y', plane=2, bin_width=30, show=True):
        tel_plane = self.Telescope.Plane(plane)
        evnt_cut, dut_cut, tel_cut = self.get_corr_cuts(plane)
        splits, x = histogram(self.Time[evnt_cut], bins.get_time(self.Time, bin_width)[1])
        y, z_ = (self.get_x(cut=dut_cut), self.get_x(tel_plane, cut=tel_cut)) if mode.lower() == 'x' else (self.get_y(cut=dut_cut), self.get_y(tel_plane, cut=tel_cut))
        values = array([corrcoef(arr, rowvar=False)[0][1] for arr in split(array([y, z_]).T, cumsum(splits)[:-1])])
        g = self.make_tgrapherrors('gal', 'Event Alignment', x=x[:-1] + diff(x) / 2, y=values)
        format_histo(g, y_tit='Correlation Factor', y_off=1.6, **self.get_time_args())
        self.draw_histo(g, show, .13, draw_opt='ap')

    def draw_cluster_size(self, show=True, cut=None, trk_cut=-1):
        self.format_statbox(all_stat=True)
        v = self.get_cluster_size(cut, trk_cut)
        self.draw_disto(v, bins.make(0, 10), 'Cluster Size in {}'.format(self.Plane), show=show, x_tit='Cluster Size', lm=.14, y_off=2)

    def draw_cluster_size_map(self, res=.3, local=True, cut=None, fid=False, show=True):
        self.format_statbox(entries=True, x=.84)
        cut = self.Cuts(cut) if fid else self.Cuts.get_special('fid')
        x, y = self.Tracks.get_coods(local, cut)
        self.draw_prof2d(x, y, self.get_cluster_size(cut), bins.get_coods(local, self.Plane, res), 'Cluster Size', show=show, z_tit='Cluster Size', **self.get_ax_tits(local))

    def draw_hit_map(self, res=.3, local=True, cut=None, show=True):
        self.format_statbox(entries=True, x=.84)
        x, y = self.Tracks.get_coods(local, cut)
        self.draw_histo_2d(x, y, bins.get_coods(local, self.Plane, res), 'Hit Map', show=show, **self.get_ax_tits(local))

    def draw_ped_map(self, c_max, c_min=None, cut=None):
        charge = self.get_charges(cut=False)
        cut = Cut('charge', charge < c_max if c_min is None else (charge > c_min) & (charge < c_max)) + self.Cuts(cut)
        self.draw_hit_map(cut=cut)

    def draw_charge_map(self, res=.3, cluster=False, cut=None):
        res = 1 if cluster else res
        x, y = (self.get_x(cut), self.get_y(cut)) if cluster else (self.Tracks.get_x(cut), self.Tracks.get_y(cut))
        self.format_statbox(entries=True, x=.84)
        self.draw_prof2d(x, y, self.get_charges(cut=cut), bins.get_local(self.Plane, res), 'Charge Map', x_tit='Column', y_tit='Row', z_tit='Charge [vcal]')

    def draw_charge_distribution(self, bin_width=4, cut=None, x_range=None, show=True):
        self.format_statbox(all_stat=True)
        self.draw_disto(self.get_charges(cut=cut), bins.get_vcal(bin_width), 'Cluster Charge', x_tit='Charge [vcal]', x_range=x_range, show=show, y_off=1.8, lm=.12)

    def draw_signal_distribution(self, bin_width=200, x_range=None, cut=None, show=True):
        self.format_statbox(all_stat=True)
        values = self.get_charges(cut=cut) * self.DUT.VcalToEl
        return self.draw_disto(values, bins.get_electrons(bin_width), 'Pulse Height', x_tit='Pulse Height [e]', x_range=x_range, show=show)

    def draw_trigger_phase(self, cut=None, raw=False):
        self.format_statbox(entries=True)
        h = self.draw_disto(self.get_trigger_phase(raw, cut), bins.make(0, 11), 'Trigger Phase', x_tit='Trigger Phase', y_off=1.8, lm=.13)
        format_histo(h, y_range=[0, h.GetMaximum() * 1.1])
        update_canvas()

    def draw_charge_vs_trigger_phase(self, cut=None, show=True):
        self.format_statbox(entries=True)
        x, y = self.get_trigger_phase(cut=cut), self.get_charges(cut=cut)
        self.draw_prof(x, y, bins.make(0, 10, last=True), 'Charge vs. Trigger Phase', show=show, x_tit='Trigger Phase', y_tit='Charge [vcal]')

    def draw_charge_trend(self, bin_width=30, e=False, y_range=None, cut=None, show=True, stats=True):
        self.format_statbox(entries=True, exe=stats)
        t, charges = self.get_time(cut), self.get_charges(e=e, cut=cut)
        p = self.draw_prof(t, charges, bins.get_time(t, bin_width), x_tit='Time [hh:mm}', y_tit='Charge [{}]'.format('e' if e else 'vcal'), show=show, t_ax_off=0)
        values = get_hist_vec(p)
        format_histo(p, y_range=choose(y_range, increased_range([min(values).n, max(values).n], 1, 2)))
        update_canvas()
        return p

    def fit_charge(self, bin_width=30, y_range=None, show=False):
        self.format_statbox(only_fit=True)
        h = self.draw_charge_trend(bin_width=bin_width, y_range=y_range, show=show, stats=False)
        format_histo(h, stats=1, name='Fit Result')
        fit = h.Fit('pol0', 'sq')
        return fit

    def draw_efficiency(self, bin_width=30, show=True):
        t, e = self.get_time(trk_cut=None), self.get_efficiency()
        return self.draw_prof(t, e, bins.get_time(t, bin_width), x_tit='Time [hh:mm]', y_tit='Efficiency [%]', t_ax_off=0, y_range=[0, 105], show=show, stats=0)

    def fit_efficiency(self, bin_width=30):
        self.format_statbox(only_fit=True)
        h = self.draw_efficiency(bin_width=bin_width)
        format_histo(h, stats=1, name='Fit Result')
        fit = h.Fit('pol0', 'sq')
        return fit

    def draw_efficiency_vs_trigger_phase(self, show=True):
        self.format_statbox(entries=True)
        cut = self.Tracks.Cuts.get_special('triggerphase')
        x, y = self.get_trigger_phase(trk_cut=cut), self.get_efficiency(cut)
        return self.draw_prof(x, y, bins.make(0, 11), x_tit='Trigger Phase', y_tit='Efficiency [%]', y_range=[0, 105], show=show)

    def draw_efficiency_map(self, res=.3, local=True):
        x, y = self.Tracks.get_x(trk_cut=None), self.Tracks.get_y(trk_cut=None)
        self.format_statbox(entries=True, x=.84)
        self.draw_prof2d(x, y, self.get_efficiency(), bins.get_coods(local, self.Plane, res), **self.get_ax_tits(local))

    def draw_time(self, show=True):
        t = self.get_time()
        g = self.make_graph_from_profile(self.draw_prof(arange(t.size), t, bins.make(0, t.size, sqrt(t.size)), show=False), err=0)
        format_histo(g, x_tit='Event Number', y_tit='Time [hh:mm]', y_off=1.8)
        set_time_axis(g, axis='Y')
        self.draw_histo(g, show, .13, draw_opt='al')

    # endregion DRAW
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
