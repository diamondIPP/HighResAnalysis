#!/usr/bin/env python
# --------------------------------------------------------
#       small script to read simple text files written by pXar
# created on August 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from ROOT import TProfile2D, TH1F, TH2F
from cern_run import CERNRun
from desy_run import DESYRun
from currents import Currents
from desy_converter import DESYConverter
from converter import Converter
from fit import *
import bins
from tracks import TrackAnalysis
from telescope import TelescopeAnalysis
from analysis import *
from numpy import in1d, cumsum, split, concatenate, unique
from calibration import Calibration


class DUTAnalysis(Analysis):

    def __init__(self, run_number, dut_number, test_campaign, single_mode=True, verbose=True):

        Analysis.__init__(self, test_campaign, verbose)
        self.print_start(run_number)

        # MAIN
        self.Run = self.init_run()(run_number, dut_number, self.TCDir, self.Config, single_mode)
        self.DUT = self.Run.DUT
        self.Plane = self.DUT.Plane

        # DATA
        self.Converter = self.init_converter()(self.TCDir, self.Run.Number, self.Config)
        self.Data = self.load_file()

        # INFO
        self.NEvents = self.get_entries()
        self.StartTime = self.get_start_time()
        self.EndTime = self.get_end_time()

        # SUBCLASSES
        # TODO: add cut class
        self.Calibration = Calibration(self.Run)
        self.Telescope = TelescopeAnalysis(self)
        self.Tracks = TrackAnalysis(self)
        self.Currents = Currents(self)

    # ----------------------------------------
    # region INIT
    def init_run(self):
        return DESYRun if self.Location == 'DESY' else CERNRun

    def init_converter(self):
        return DESYConverter if self.Location == 'DESY' else Converter

    def get_entries(self):
        return self.Data['Tracks']['NTracks'].size

    def load_file(self):
        self.Converter.run()
        try:
            f = h5py.File(self.Run.FileName, 'r')
            _ = f['Tracks']
            return f
        except (KeyError, OSError):
            warning('could not load data file {} -> start with dummy'.format(self.Run.FileName))
            return {'Tracks': array([])}  # dummy

    def reload_file(self):
        self.Data = self.load_file()
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_cluster_cut(self):
        return where(self.get_data('Clusters', 'NClusters') == 1)

    def get_start_time(self):
        return datetime.fromtimestamp(self.Run.StartTime)

    def get_end_time(self):
        return datetime.fromtimestamp(self.Run.EndTime)

    def get_time(self, cut=None):
        t = array(self.Data['Event']['Time']).astype('f8') + self.Run.StartTime
        return t if cut is None else t[cut]

    def get(self, plane=None):
        return self.Data['Plane{}'.format(self.get_plane(plane).Number)]

    def get_group(self, grp, plane=None):
        return self.Data['Plane{}'.format(self.get_plane(plane).Number)][grp]

    def get_data(self, grp, key=None, plane=None, cut=None):
        data = self.get_group(grp, plane)
        data = array(data) if key is None else array(data[key])
        return data if cut is None else data[cut]

    def get_n(self, name='Hits', cut=None):
        return self.get_data(name, 'N{}'.format(name), cut=cut)

    def get_split(self, cluster=True):
        n = 'Clusters' if cluster else 'Hits'
        return cumsum(self.get_n(n))[:-1]

    def get_charges(self, cluster=True, e=False, cut=None, flat=False):
        values = self.get_data('Clusters' if cluster else 'Hits', 'Charge', self.Plane, cut if flat else None) * (self.DUT.VcalToEl if e else 1)
        if flat:
            return values
        return concatenate(array(split(values, self.get_split(cluster)))[cut])

    def get_charge(self, cut=None):
        values = self.get_charges(cut=cut, flat=True)
        m, s = mean_sigma(values)
        return ufloat(m, s / sqrt(values.size))

    def get_x(self, plane=None, cluster=True, cut=None):
        return self.get_data('Clusters' if cluster else 'Hits', 'X', plane, cut)

    def get_y(self, plane=None, cluster=True, cut=None):
        return self.get_data('Clusters' if cluster else 'Hits', 'Y', plane, cut)

    def get_hits(self, plane=None, cut=None):
        return array([self.get_x(plane, cluster=False, cut=cut), self.get_y(plane, cluster=False, cut=cut)])

    def get_clusters(self, plane=None, cut=None):
        return array([self.get_x(plane, cluster=True, cut=cut), self.get_y(plane, cluster=True, cut=cut)])

    def get_mask(self, plane=None):
        return self.get_data('Mask', plane=plane)

    def get_mask_cut(self, plane=None):
        x, y = self.get_hits(plane)
        data = x.astype('i') * 10000 + y  # make unique number out of the tuple... Is there a way to compare tuples?
        mx, my = self.get_mask(plane).T.astype('i')
        mask = mx * 10000 + my
        return where(in1d(data, mask, invert=True))[0]

    def get_plane(self, plane):
        return choose(plane, self.Plane)
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw_n(self, plane, name, show=True):
        self.format_statbox(all_stat=True)
        n, pl = name, self.get_plane(plane)
        self.draw_disto(self.get_data(n, 'N{}'.format(n), pl), 'Number of {} in {}'.format(n, pl), bins.make(0, 30), lm=.13, show=show, x_tit='Number of {}'.format(n), y_off=2)

    def draw_mask(self, plane=None, show=True):
        plane = self.get_plane(plane)
        h = TH2F('htm', 'Masked Pixels in {}'.format(plane), *bins.get_local(self.Plane))
        fill_hist(h, *self.get_mask(plane).T)
        format_histo(h, x_tit='Column', y_tit='Row', y_off=1.3, fill_color=1)
        self.format_statbox(entries=True)
        self.draw_histo(h, show=show, lm=.12, draw_opt='box')

    def draw_occupancy(self, plane=None, cluster=True, cut=None, show=True):
        plane = self.get_plane(plane)
        h = TH2F('hto', '{} Occupancy in {}'.format('Cluster' if cluster else 'Hit', plane), *bins.get_local(plane))
        cut = unique(concatenate([self.get_mask_cut(plane), [] if cut is None else cut])) if not cluster else cut
        fill_hist(h, *(self.get_clusters(plane, cut) if cluster else self.get_hits(plane, cut)))
        format_histo(h, x_tit='Column', y_tit='Row', y_off=1.5, z_tit='Number of Entries', z_off=1.2)
        self.format_statbox(entries=True, x=.83, m=True)
        self.draw_histo(h, show=show, lm=.12, draw_opt='colz', rm=.15)

    def draw_n_hits(self, plane=None, show=True):
        self.draw_n(plane, 'Hits', show)

    def draw_n_clusters(self, plane=None, show=True):
        self.draw_n(plane, 'Clusters', show)

    def draw_cluster_size(self, plane=None, show=True):
        self.format_statbox(all_stat=True)
        v = self.get_data('Clusters', 'Size', plane)
        self.draw_disto(v[v > 0], 'Cluster Size in {}'.format(self.get_plane(plane)), bins.make(0, 10), show=show, x_tit='Cluster Size', lm=.14, y_off=2)

    def draw_n_intercepts(self, plane=None, show=True):
        self.draw_n(plane, 'Intercepts', show)

    def draw_charge_map(self, res=1, plane=None, cut=None):
        plane = self.get_plane(plane)
        p = TProfile2D('pam', 'Charge Map', *bins.get_local(plane, res))
        fill_hist(p, self.get_x(plane, cut=cut), self.get_y(plane, cut=cut), self.get_charges(cut=cut, flat=True))
        self.format_statbox(entries=True, x=.78)
        format_histo(p, x_tit='Cluster X', y_tit='Cluster Y', z_tit='Charge [vcal]', y_off=1.2, z_off=1.5)
        self.draw_histo(p, draw_opt='colz', rm=.18)

    def draw_charge_distribution(self, bin_width=4, cut=None, x_range=None, show=True):
        self.format_statbox(all_stat=True)
        self.draw_disto(self.get_charges(cut=cut, flat=True), 'Cluster Charge', bins.get_vcal(bin_width), x_tit='Charge [vcal]', x_range=x_range, show=show)

    def draw_signal_distribution(self, bin_width=200, x_range=None, cut=None, show=True):
        self.format_statbox(all_stat=True)
        values = self.get_charges(cut=cut, flat=True) * self.DUT.VcalToEl
        return self.draw_disto(values, 'Pulse Height', bins.get_electrons(bin_width), x_tit='Pulse Height [e]', x_range=x_range, show=show)

    def draw_trigger_phase(self):
        self.format_statbox(all_stat=True)
        self.draw_disto(self.get_data('TriggerPhase'), 'Trigger Phase', bins.make(0, 10), x_tit='Trigger Phase', y_off=1.8, lm=.13)

    def draw_charge_vs_trigger_phase(self):
        pass

    def draw_charge_trend(self, bin_width=30, e=False, y_range=None, show=True, stats=True):
        self.format_statbox(entries=True, exe=stats)
        cut = self.get_cluster_cut()
        t = self.get_time(cut)
        charges = self.get_charges(e=e, cut=cut)
        p = self.draw_prof(t, charges, bins.make(t[0], t[-1], bin_width), x_tit='Time [hh:mm}', y_tit='Charge [{}]'.format('e' if e else 'vcal'), show=show, t_ax_off=0)
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
        cut = self.Tracks.get_n() > 0
        eff = (self.get_n(cut=cut) > 0).astype('u2') * 100
        t = self.get_time(cut)
        return self.draw_prof(t, eff, bins.make(t[0], t[-1], bin_width), x_tit='Time [hh:mm]', y_tit='Efficiency [%]', t_ax_off=0, y_range=[0, 105], show=show, stats=0)

    def fit_efficiency(self, bin_width=30):
        self.format_statbox(only_fit=True)
        h = self.draw_efficiency(bin_width=bin_width)
        format_histo(h, stats=1, name='Fit Result')
        fit = h.Fit('pol0', 'sq')
        return fit

    def draw_efficiency_vs_trigger_phase(self):
        # add profile
        pass

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
    args = aparser.parse_args()
    z = DUTAnalysis(args.run, args.dut, test_campaign=args.testcampaign, single_mode=args.single_mode, verbose=args.verbose)
    cal = z.Calibration
    c = z.Converter
    t = z.Tracks
    z.add_info(t_start, prnt=True)
