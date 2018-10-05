#!/usr/bin/env python
# --------------------------------------------------------
#       small script to read simple text files written by pxar
# created on August 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from draw import *
from utils import *
from numpy import zeros
from ROOT import TH2I, TH1I, TProfile2D, TF1, TCut, TProfile
from progressbar import Bar, ETA, FileTransferSpeed, Percentage, ProgressBar
from argparse import ArgumentParser
from os.path import dirname, join, realpath
from sys import stdout
from pickle import load, dump
from json import load as jload


class ReadData(Draw):

    def __init__(self, filename):
        Draw.__init__(self)

        # progress bar
        self.Widgets = ['Progress: ', Percentage(), ' ', Bar(marker='>'), ' ', ETA(), ' ', FileTransferSpeed()]
        self.ProgressBar = None

        self.Dir = dirname(dirname(realpath(__file__)))
        self.RunNumber = int(filename.strip('.root').split('_')[-1])
        self.FileName = filename
        self.File = read_root_file(filename)
        self.Tree = self.File.Get('Hits')

        self.FitParameters = zeros((52, 80, 4))
        self.Fit = TF1('ErFit', '[3] * (TMath::Erf((x - [0]) / [1]) + [2])', -500, 255 * 7)
        # self.get_calibration_data(filename)

        self.NEntries = self.Tree.GetEntries()
        self.StartTime = self.get_start_time()
        self.EndTime = self.get_end_time()

        # self.RunInfo = self.load_run_info()

    def load_run_info(self):
        f = open(join(self.Dir, 'runlog.json'))
        info = jload(f, parse_int=True)
        for key, value in info.iteritems():
            if value['Run Number']['CMS'] == str(self.RunNumber):
                return info[key]


    def draw_number_of_hits(self):
        h = TH1I('hnh', 'Number of Hits', 4180, 0, 4180)
        self.Tree.Draw('NHits>>hnh', '', 'goff')
        self.format_statbox(only_entries=True)
        x_range = [0, h.GetBinCenter(h.FindLastBinAbove(0)) + 1]
        self.format_histo(h, x_tit='Number of Hits', y_tit='Number of Entries', y_off=1.3, x_range=x_range)
        self.draw_histo(h, draw_opt='colz', logy=True)

    def draw_hitmap(self, start=0, n=1e9, cluster=True, res=1, cut=''):
        h = TH2I('hhm', 'Hit Map', int(52 * res), .5, 52.5, int(80 * res), .5, 80.5)
        self.Tree.Draw('{0}Y:{0}X>>hhm'.format('Cluster' if cluster else 'Pix'), TCut(cut), 'goff', int(n), start)
        self.format_statbox(only_entries=True, x=.78)
        self.format_histo(h, x_tit='Column', y_tit='Row', z_tit='Number of Entries', y_off=1.2, z_off=1.5)
        self.draw_histo(h, draw_opt='colz', rm=.18)

    def draw_charge_map(self, res=1, cut=''):
        h = TProfile2D('pam', 'Charge Map', int(52 * res), .5, 52.5, int(80 * res), .5, 80.5)
        self.Tree.Draw('ClusterVcal:ClusterY:ClusterX>>pam', TCut(cut), 'goff')
        self.format_statbox(only_entries=1, x=.78)
        self.format_histo(h, x_tit='Column', y_tit='Row', z_tit='Charge [vcal]', y_off=1.2, z_off=1.5)
        self.draw_histo(h, draw_opt='colz', rm=.18)

    def draw_charge_distribution(self, vcal=True, cut='', x_range=None):
        h = TH1I('had', 'ADC Distribution', 3756 / 4, -256, 3500)
        self.Tree.Draw('{}>>had'.format('ClusterVcal' if vcal else ''), TCut(cut), 'goff')
        self.format_statbox(only_entries=1)
        x_range = [h.GetBinCenter(i) for i in [h.FindFirstBinAbove(0) - 1, h.FindLastBinAbove(0) + 1]] if x_range is None else x_range
        self.format_histo(h, x_tit='VCAL', y_tit='Number of Entries', y_off=1.2, x_range=x_range)
        self.draw_histo(h)

    def draw_cluster_size(self):
        h = TH1I('hcs', 'Cluster Size', 20, 0, 20)
        self.Tree.Draw('ClusterSize>>hcs', '', 'goff')
        self.format_histo(h, x_tit='Cluster Size', y_tit='Number of Entries', y_off=1.2, stats=0)
        self.draw_histo(h)

    def draw_trigger_phase(self):
        h = TH1I('htp', 'Trigger Phase', 10, 0, 10)
        self.Tree.Draw('Timing>>htp', '', 'goff')
        self.format_histo(h, x_tit='Trigger Phase', y_tit='Number of Entries', y_off=1.2, stats=0, fill_color=self.FillColor, x_range=[0, h.GetMaximum() * 1.2])
        self.draw_histo(h)

    def draw_charge_vs_trigger_phase(self):
        p = TProfile('pctp', 'Charge @ Trigger Phase', 10, 0, 10)
        self.Tree.Draw('ClusterVcal:Timing>>pctp', '', 'goff')
        self.format_histo(p, x_tit='Trigger Phase', y_tit='Charge [vcal]', y_off=1.2, stats=0)
        self.draw_histo(p)

    def draw_charge_vs_time(self, bin_width=30, y_range=None, show=True):
        p = TProfile('pct', 'Charge vs. Time - Run {}'.format(self.RunNumber), int((self.EndTime - self.StartTime) / bin_width), self.StartTime, self.EndTime)
        self.Tree.Draw('ClusterVcal:TimeStamp / 1000>>pct', '', 'goff')
        y_range = [0, 1.2 * p.GetMaximum()] if y_range is None else y_range
        self.format_histo(p, x_tit='Time [hh:mm]', y_tit='Charge [vcal]', y_off=1.4, stats=0, fill_color=self.FillColor, t_ax_off=0, y_range=y_range)
        self.draw_histo(p, lm=.11, show=show)
        return p

    def fit_charge(self, bin_width=30, y_range=None):
        set_statbox(only_fit=True, y=.5)
        h = self.draw_charge_vs_time(bin_width=bin_width, y_range=y_range)
        self.format_histo(h, stats=1, name='Fit Result')
        fit = h.Fit('pol0', 'sq')
        return fit

    def draw_efficiency(self, bin_width=30, cut='', show=True):
        p = TProfile('pet', 'Trigger Efficiency - Run {}'.format(self.RunNumber), int((self.EndTime - self.StartTime) / bin_width), self.StartTime, self.EndTime)
        self.Tree.Draw('(NHits>0)*100.:TimeStamp / 1000>>pet', TCut(cut), 'goff')
        self.format_histo(p, x_tit='Time [hh:mm]', y_tit='Efficiency [%]', y_off=1.4, stats=0, fill_color=self.FillColor, t_ax_off=0, y_range=[0,105])
        self.draw_histo(p, lm=.11, show=show)
        return p

    def fit_efficiency(self, bin_width=30, cut=''):
        set_statbox(only_fit=True, y=.3)
        h = self.draw_efficiency(bin_width=bin_width, cut=cut)
        self.format_histo(h, stats=1, name='Fit Result')
        fit = h.Fit('pol0', 'sq')
        return fit

    def draw_efficiency_vs_trigger_phase(self):
        p = TProfile('petp', 'Trigger Efficiency @ Trigger Phase', 10, 0, 10)
        self.Tree.Draw('(NHits>0)*100.:Timing>>petp', '', 'goff')
        self.format_histo(p, x_tit='Trigger Phase', y_tit='Trigger Efficiency [%]', y_off=1.2, stats=0, y_range=[0, 105])
        self.draw_histo(p)

    def get_start_time(self):
        self.Tree.GetEntry(0)
        return self.Tree.TimeStamp / 1000

    def get_end_time(self):
        self.Tree.GetEntry(self.NEntries - 1)
        t = self.Tree.TimeStamp / 1000
        return t

    def get_calibration_data(self, filename):
        pickle_name = 'fitpars.pickle'
        if file_exists(pickle_name):
            f = open(pickle_name, 'r')
            self.FitParameters = load(f)
            f.close()
        else:
            f = open(join(dirname(dirname(filename)), 'phCalibration_C0.dat'))
            f.readline()
            low_range = [int(val) for val in f.readline().split(':')[-1].split()]
            high_range = [int(val) for val in f.readline().split(':')[-1].split()]
            x = low_range + [val * 7 for val in high_range]
            f.readline()
            self.Fit.SetParameters(309.2062, 112.8961, 1.022439, 35.89524)
            for line in f.readlines():
                data = line.split('Pix')
                y = [int(val) for val in data[0].split()]
                x1 = [ufloat(ix, 1) for (ix, iy) in zip(x, y) if iy]
                y1 = [ufloat(iy, 1) for iy in y if iy]
                g = self.make_tgrapherrors('gcal', 'gcal', x=x1, y=y1)
                g.Fit(self.Fit, 'q', '', 0, 3000)
                print '\r{}'.format(data[-1]).strip('\n'),
                stdout.flush()
                col, row = [int(val) for val in data[-1].split()]
                self.FitParameters[col][row] = [self.Fit.GetParameter(i) for i in xrange(4)]
            fp = open(pickle_name, 'w')
            dump(self.FitParameters, fp)
            fp.close()
            f.close()

    def draw_calibration_fit(self, col=14, row=14):
        f = open(join(dirname(dirname(self.FileName)), 'phCalibration_C0.dat'))
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
                self.format_histo(g, x_tit='vcal', y_tit='adc', y_off=1.4)
                self.draw_histo(g)
                self.Fit.Draw('same')
                break

    def start_pbar(self, n):
        self.ProgressBar = ProgressBar(widgets=self.Widgets, maxval=n)
        self.ProgressBar.start()


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('filename', nargs='?', default='')
    args = p.parse_args()
    z = ReadData(args.filename)
