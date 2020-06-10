#!/usr/bin/env python
# --------------------------------------------------------
#       small script to read simple text files written by pXar
# created on August 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from __future__ import print_function

from analysis import *
from ROOT import TH2I, TH1I, TProfile2D, TCut, TProfile, TH1F, TF1
from ROOT import TBrowser
from cern_run import CERNRun
from desy_run import DESYRun
from pickle import load as pload
from currents import Currents
from desy_converter import DESYConverter
from converter import Converter
from fit import *
from binning import *
from tracks import TrackAnalysis


class DUTAnalysis(Analysis):

    def __init__(self, run_number, dut, test_campaign, single_mode=True, verbose=True):

        self.RunNumber = run_number
        Analysis.__init__(self, test_campaign, verbose)

        self.Run = self.init_run()(run_number, dut, self.TCDir, self.Config, single_mode)
        self.DUT = self.Run.DUT
        self.Data = self.Run.Data
        self.Bins = Bins(self.Config)

        self.Converter = self.init_converter()(self.TCDir, self.Run.Number, self.Config)

        self.NEntries = self.get_entries()
        self.StartTime = self.get_start_time()
        self.EndTime = self.get_end_time()

        # Subclasses
        self.Tracks = TrackAnalysis(self)
        self.Currents = Currents(self)

        # TODO add calibration in hdf5 file (no charge for MIMOSA anyway)
        # Calibration
        self.Fit = TF1('ErFit', '[3] * (TMath::Erf((x - [0]) / [1]) + [2])', -500, 255 * 7)
        # self.FitParameters = None
        # self.get_calibration_data()
        self.print_start(self.RunNumber)

    # ----------------------------------------
    # region INIT
    def init_run(self):
        return DESYRun if self.Location == 'DESY' else CERNRun

    def init_converter(self):
        return DESYConverter if self.Location == 'DESY' else Converter

    def get_entries(self):
        pass
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region GET
    def get_start_time(self):
        return self.Run.StartTime

    def get_end_time(self):
        return self.Run.EndTime

    def get_calibration_number(self):
        numbers = sorted(int(remove_letters(basename(name))) for name in glob(join(self.TCDir, self.DUT.Name, 'calibrations', 'phCal*.dat')) if basename(name)[5].isdigit())
        first_run = int(self.Run.RunLogs.keys()[0])
        if first_run > numbers[-1]:
            return numbers[-1]
        next_number = next(nr for nr in numbers if nr >= first_run)
        return numbers[numbers.index(next_number) - 1]

    def get_calibration_data(self):
        pickle_name = join(self.TCDir, self.DUT.Name, 'calibrations', 'fitpars{}.pickle'.format(self.get_calibration_number()))
        with open(pickle_name, 'r') as f:
            self.FitParameters = pload(f)
    # endregion GET
    # ----------------------------------------

    # ----------------------------------------
    # region DRAW
    def draw_current(self, v_range=None, f_range=None, c_range=None, show=True, draw_opt='al'):
        self.Currents.draw_indep_graphs(v_range=v_range, f_range=f_range, c_range=c_range, show=show, draw_opt=draw_opt)

    def draw_number_of_hits(self):
        h = TH1I('hnh', 'Number of Hits', 4180, 0, 4180)
        self.Tree.Draw('NHits>>hnh', '', 'goff')
        self.format_statbox(entries=True)
        x_range = [0, h.GetBinCenter(h.FindLastBinAbove(0)) + 1]
        format_histo(h, x_tit='Number of Hits', y_tit='Number of Entries', y_off=1.3, x_range=x_range)
        self.draw_histo(h, draw_opt='colz', logy=True)

    def draw_occupancy(self, start=0, n=1e9, cluster=True, res=1, cut=''):
        h = TH2I('hhm', 'Hit Map', int(52 * res), .5, 52.5, int(80 * res), .5, 80.5)
        self.Tree.Draw('{0}Y:{0}X>>hhm'.format('Cluster' if cluster else 'Pix'), TCut(cut), 'goff', int(n), start)
        self.format_statbox(entries=True, x=.78)
        format_histo(h, x_tit='Column', y_tit='Row', z_tit='Number of Entries', y_off=1.2, z_off=1.5)
        self.draw_histo(h, draw_opt='colz', rm=.18)

    def draw_charge_map(self, res=1, cut=''):
        h = TProfile2D('pam', 'Charge Map', int(52 * res), .5, 52.5, int(80 * res), .5, 80.5)
        self.Tree.Draw('ClusterVcal:ClusterY:ClusterX>>pam', TCut(cut), 'goff')
        self.format_statbox(entries=1, x=.78)
        format_histo(h, x_tit='Column', y_tit='Row', z_tit='Charge [vcal]', y_off=1.2, z_off=1.5)
        self.draw_histo(h, draw_opt='colz', rm=.18)

    def draw_charge_distribution(self, vcal=True, cut='', x_range=None):
        binning = [3756 / 4, -256, 3500] if vcal else [255, 0, 255]
        title = 'VCAL' if vcal else 'ADC'
        h = TH1I('had', '{} Distribution'.format(title), *binning)
        self.Tree.Draw('{}>>had'.format('ClusterVcal' if vcal else 'Value'), TCut(cut), 'goff')
        self.format_statbox(entries=1)
        x_range = [h.GetBinCenter(i) for i in [h.FindFirstBinAbove(0) - 1, h.FindLastBinAbove(0) + 1]] if x_range is None else x_range
        format_histo(h, x_tit=title, y_tit='Number of Entries', y_off=1.2, x_range=x_range)
        self.draw_histo(h)

    def draw_signal_distribution(self, x_range=None, cut=None, threshold=80 * 46.5, show=True):
        h = TH1F('hph', 'Pulse Height Distribution', 50000 / 200, -5000, 45000)
        self.Tree.Draw('ClusterVcal * 46.5 >> hph', self.Cut() if cut is None else cut, 'goff')
        self.format_statbox(all_stat=True)
        format_histo(h, x_tit='Pulse Height [e]', y_tit='Number of Entries', y_off=1.6, x_range=x_range, ndivx=505, fill_color=self.FillColor)
        self.draw_histo(h, lm=.12, show=show)
        if threshold and show:
            self.draw_y_axis(threshold, h.GetYaxis().GetXmin(), h.GetMaximum(), 'threshold #approx {}e  '.format(int(round_down_to(threshold, 100))), off=.3, line=True, opt='-L')
        return h

    def draw_cluster_size(self):
        h = TH1I('hcs', 'Cluster Size', 20, 0, 20)
        self.Tree.Draw('ClusterSize>>hcs', '', 'goff')
        format_histo(h, x_tit='Cluster Size', y_tit='Number of Entries', y_off=2, stats=0)
        self.draw_histo(h, lm=.14)

    def draw_trigger_phase(self):
        h = TH1I('htp', 'Trigger Phase', 10, 0, 10)
        self.Tree.Draw('Timing>>htp', '', 'goff')
        format_histo(h, x_tit='Trigger Phase', y_tit='Number of Entries', y_off=1.2, stats=0, fill_color=self.FillColor, x_range=[0, h.GetMaximum() * 1.2])
        self.draw_histo(h)

    def draw_charge_vs_trigger_phase(self):
        prof = TProfile('pctp', 'Charge @ Trigger Phase', 10, 0, 10)
        self.Tree.Draw('ClusterVcal:Timing>>pctp', '', 'goff')
        format_histo(prof, x_tit='Trigger Phase', y_tit='Charge [vcal]', y_off=1.2, stats=0)
        self.draw_histo(prof)

    def draw_charge_vs_time(self, bin_width=30, y_range=None, show=True):
        prof = TProfile('pct', 'Charge vs. Time - Run {}'.format(self.RunNumber), int((self.EndTime - self.StartTime) / bin_width), self.StartTime, self.EndTime)
        self.Tree.Draw('ClusterVcal:TimeStamp / 1000>>pct', '', 'goff')
        y_range = [0, 1.2 * prof.GetMaximum()] if y_range is None else y_range
        format_histo(prof, x_tit='Time [hh:mm]', y_tit='Charge [vcal]', y_off=1.4, stats=0, fill_color=self.FillColor, t_ax_off=0, y_range=y_range)
        self.draw_histo(prof, lm=.11, show=show)
        return prof

    def fit_charge(self, bin_width=30, y_range=None):
        self.format_statbox(only_fit=True, y=.5)
        h = self.draw_charge_vs_time(bin_width=bin_width, y_range=y_range)
        format_histo(h, stats=1, name='Fit Result')
        fit = h.Fit('pol0', 'sq')
        return fit

    def draw_efficiency(self, bin_width=30, cut='', show=True):
        prof = TProfile('pet', 'Trigger Efficiency - Run {}'.format(self.RunNumber), int((self.EndTime - self.StartTime) / bin_width), self.StartTime, self.EndTime)
        self.Tree.Draw('(NHits>0)*100.:TimeStamp / 1000>>pet', TCut(cut), 'goff')
        format_histo(prof, x_tit='Time [hh:mm]', y_tit='Efficiency [%]', y_off=1.4, stats=0, fill_color=self.FillColor, t_ax_off=0, y_range=[0, 105])
        self.draw_histo(prof, lm=.11, show=show)
        return prof

    def fit_efficiency(self, bin_width=30, cut=''):
        self.format_statbox(only_fit=True, y=.3)
        h = self.draw_efficiency(bin_width=bin_width, cut=cut)
        format_histo(h, stats=1, name='Fit Result')
        fit = h.Fit('pol0', 'sq')
        return fit

    def draw_efficiency_vs_trigger_phase(self):
        prof = TProfile('petp', 'Trigger Efficiency @ Trigger Phase', 10, 0, 10)
        self.Tree.Draw('(NHits>0)*100.:Timing>>petp', '', 'goff')
        format_histo(prof, x_tit='Trigger Phase', y_tit='Trigger Efficiency [%]', y_off=1.2, stats=0, y_range=[0, 105])
        self.draw_histo(prof)
    # endregion DRAW
    # ----------------------------------------

    def draw_calibration_fit(self, col=31, row=55):
        f = open(join(self.TCDir, self.Run.DUTName, 'calibrations', 'phCal{}.dat'.format(self.get_calibration_number())))
        f.readline()
        low_range = [int(val) for val in f.readline().split(':')[-1].split()]
        high_range = [int(val) for val in f.readline().split(':')[-1].split()]
        x = [ufloat(xval, 1) for xval in low_range + [val * 7 for val in high_range]]
        f.readline()
        self.Fit.SetParameters(*self.FitParameters[col][row])
        for line in f.readlines():
            data = line.split('Pix')
            icol, irow = [int(val) for val in data[-1].split()]
            if col == icol and row == irow:
                y = [ufloat(int(val), 1) for val in data[0].split()]
                g = self.make_tgrapherrors('gcal', 'Calibration Fit for Pix {} {}'.format(col, row), x=x, y=y)
                format_histo(g, x_tit='vcal', y_tit='adc', y_off=1.4)
                self.draw_histo(g)
                self.Fit.Draw('same')
                break

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
    p = ArgumentParser()
    p.add_argument('run', nargs='?', default=11, type=int)
    p.add_argument('dut', nargs='?', default=1, type=int)
    p.add_argument('--testcampaign', '-tc', nargs='?', default=None)
    p.add_argument('--verbose', '-v', action='store_true')
    p.add_argument('--single_mode', '-s', action='store_false')
    args = p.parse_args()
    y = DUTAnalysis(args.run, args.dut, test_campaign=args.testcampaign, single_mode=args.single_mode, verbose=args.verbose)
    z = y.Converter
    t = y.Tracks
