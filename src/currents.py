from os.path import getsize
from time import sleep
from ROOT import TGraph
from numpy import genfromtxt, isnan, datetime64, invert, where, concatenate, char, sign, uint32
from pytz import timezone
import bins
from analysis import *
from utils import *

# ====================================
# CONSTANTS
axis_title_size = 0.06
label_size = .05
title_offset = 0.8
col_vol = 602  # 807
col_cur = 899  # 418
pad_margins = [.065, .09, .15, .1]
marker_size = .3


class Currents(Analysis):
    """reads in information from the keithley log file"""

    def __init__(self, analysis=None, test_campaign=None, dut=None, begin=None, end=None, averaging=None, verbose=False):
        Analysis.__init__(self, test_campaign if analysis is None else analysis.TestCampaign, verbose=verbose)

        # Settings
        self.Averaging = averaging
        self.TimeZone = timezone('Europe/Zurich')
        self.DataDir = join(self.TCDir, 'hvdata')

        # Config
        self.Ana = analysis
        self.IsCollection = hasattr(analysis, 'Runs')
        self.Collection = None
        self.RunNumber = self.load_run_number()
        self.RunPlan = self.load_run_plan()  # required for plotting
        self.RunLogs = self.Ana.Run.RunLogs
        self.Run = self.Ana.Run
        self.HVConfig = load_config(join(self.DataDir, 'config'))
        self.set_save_directory('currents')
        self.Bias = self.load_bias()

        # Times
        self.Begin, self.End = self.load_times(begin, end)

        # DUT
        self.DUTNumber = choose(self.Run.DUT.Number, dut, self.Ana)
        self.DUTName = self.get_dut_name()

        # HV Device Info
        self.Number = self.get_device_number()
        self.Channel = self.get_device_channel()
        self.Name = self.HVConfig.get('HV{}'.format(self.Number), 'name')
        self.Brand = remove_digits(self.Name.split('-')[0])
        self.Model = self.HVConfig.get('HV{}'.format(self.Number), 'model')
        self.Precision = .005 if '237' in self.Name else .05

        # data
        self.IgnoreJumps = True
        self.Data = self.load_data()

        # plotting
        self.Graphs = []

    # ----------------------------------------
    # region INIT
    def load_data(self):
        data_file = join(self.DataDir, 'data.hdf5')
        if not file_exists(data_file):
            self.convert_data()
        data = h5py.File(data_file, 'r')['{}_CH{}'.format(self.Name, self.Channel)]
        data = data[where((data['timestamps'] >= time_stamp(self.Begin, off=True)) & (data['timestamps'] <= time_stamp(self.End, off=True)))]
        if self.IgnoreJumps:  # filter out jumps
            data = data[where(abs(data['currents'][:-1]) * 100 > abs(data['currents'][1:]))[0] + 1]  # take out the events that are 100 larger than the previous
        data['currents'] *= 1e9 * sign(mean(data['currents']))  # convert to nA and flip sign if current is negative
        if self.Ana is not None:
            data['timestamps'] -= uint32(data['timestamps'][0] - self.Ana.StartTime)  # synchronise time vectors
        return data

    def load_bias(self):
        return self.Run.DUT.Bias if hasattr(self.Run, 'Bias') else None

    def load_run_number(self):
        return None if self.Ana is None else self.Ana.RunNumber if not self.IsCollection else self.Ana.RunPlan

    def load_run_plan(self):
        return self.Collection.SelectedRunplan if self.Ana is None else self.Ana.RunPlan if self.IsCollection else None

    def load_parser(self):
        parser = ConfigParser()
        file_path = join(self.DataDir, 'config.ini')
        if not file_exists(file_path):
            critical('HV info file "{f}" does not exist'.format(f=file_path))
        parser.read(file_path)
        self.info('HV Devices: {}'.format(', '.join(name for name in parser.sections() if name.startswith('HV'))))
        return parser

    def load_times(self, begin, end):
        if self.Ana is None:
            if str(begin).isdigit():  # run number or run plan is provided
                self.Collection.select_runs_in_range(begin, end if end is not None else begin) if end or end is None else self.Collection.select_runs_from_runplan(begin)
                return self.Collection.get_start_time(), self.Collection.get_end_time()
            else:  # actual time strings are provided
                return (self.TimeZone.localize(datetime.strptime('{}-{}'.format(self.TestCampaign.year, t), '%Y-%m/%d-%H:%M:%S')) for t in [begin, end])
        return [self.TimeZone.localize(datetime.fromtimestamp(self.RunLogs[key])) for key in ['start', 'end']]

    def get_dut_name(self):
        if self.Ana is not None:
            return self.Ana.DUT.Name
        elif self.Collection.has_selected_runs():
            return self.Collection.get_diamond_names(sel=True)[0]
        return next(log['dia{}'.format(self.DUTNumber)] for log in self.RunLogs.itervalues() if (log['starttime0']) > self.Begin)

    def get_device_str(self):
        if self.Ana is not None:
            run_info = self.RunLogs
        elif self.Collection.has_selected_runs():
            run_info = self.RunLogs[str(self.Collection.get_selected_runs()[0])]
        else:
            run_info = next(log for log in self.RunLogs.itervalues() if datetime.timestamp(log['start']) > self.Begin)
        return str(run_info['hvsupply{}'.format(self.DUTNumber)])

    def get_device_number(self):
        return self.get_device_str().split('-')[0]

    def get_device_channel(self):
        words = self.get_device_str().split('-')
        return words[1] if len(words) > 1 else '0'

    def find_data_path(self):
        data_dir = join(self.Run.TCDir, 'hvdata', '{}_CH{}'.format(self.Name, self.Channel))
        if not dir_exists(data_dir):
            critical('HV data path "{}" does not exist!'.format(data_dir))
        return data_dir

    def load_time(self, t, t_log):
        t = self.TimeZone.localize(datetime.fromtimestamp(t)) if t is not None else t_log
        return t_log if t.year < 2000 or t.day != t_log.day else t

    def load_ana_start_time(self):
        ana = self.Ana if not self.IsCollection else self.Ana.FirstAnalysis
        return self.load_time(ana.Run.StartTime if hasattr(ana.Run, 'StartTime') else None, ana.Run.LogStart)

    def load_ana_end_time(self):
        ana = self.Ana if not self.IsCollection else self.Ana.LastAnalysis
        return self.load_time(ana.Run.EndTime if hasattr(ana.Run, 'EndTime') else None, ana.Run.LogEnd)
    # endregion INIT
    # ----------------------------------------

    # ----------------------------------------
    # region DATA ACQUISITION
    def get_log_date(self, name):
        log_date = ''.join(basename(name).split('_')[-6:])
        return self.TimeZone.localize(datetime.strptime(log_date, '%Y%m%d%H%M%S.log'))

    def convert_data(self):
        info('converting hv text files to hdf5 ...')
        self.PBar.start(len(glob(join(self.DataDir, '*', '*.log'))))
        f = h5py.File(join(self.DataDir, 'data.hdf5'), 'w')
        for d in glob(join(self.DataDir, '*_*')):
            arrays = []
            for file_name in glob(join(d, '*.log')):
                if getsize(file_name) == 0:
                    remove_file(file_name)
                log_date = self.get_log_date(file_name)
                data = genfromtxt(file_name, usecols=arange(3), dtype=[('timestamps', object), ('voltages', 'f2'), ('currents', 'f4')])
                data = data[where(invert(isnan(data['voltages'])))[0]]  # remove text entries
                date_times = array(log_date.strftime('%Y-%m-%d ') + char.array(data['timestamps'].astype('U'))).astype(datetime64).astype(datetime)
                data['timestamps'] = (array([time_stamp(dt, log_date.utcoffset().seconds) for dt in date_times]).astype('u4'))
                data = data.astype([('timestamps', 'u4'), ('voltages', 'f2'), ('currents', 'f4')])
                arrays.append(data)
                self.PBar.update()
            if len(arrays):
                f.create_dataset(basename(d), data=concatenate(arrays))

    # endregion DATA ACQUISITION
    # ----------------------------------------

    # ----------------------------------------
    # region PLOTTING
    def draw_profile(self, bin_width=5, show=True):
        x, y = self.Data['timestamps'], self.Data['currents']
        self.format_statbox(entries=True)
        return self.draw_prof(x, y, bins.make(x[0], x[-1], bin_width), 'Leakage Current', x_tit='Time [hh:mm]', y_tit='Current [nA]', t_ax_off=0, markersize=.7, cx=1.5,
                              cy=.75, lm=.08, y_off=.8, show=show)

    def draw_distribution(self, show=True):
        m, s = mean_sigma(self.Data['currents'])
        xmin, xmax = m - 4 * max(s, .1), m + 4 * max(s, .1)
        self.format_statbox(all_stat=True)
        return self.draw_disto(self.Data['currents'], 'Current Distribution', bins.make(xmin, xmax, self.Precision * 2), show=show, x_tit='Current [nA]')

    def get_current(self):
        if self.Ana is not None and not self.Ana.DUT.Bias:
            warning('Bias of run {} is 0!'.format(self.Ana.RunNumber))
            return ufloat(0, 0)
        else:
            h = self.draw_distribution(show=False)
            if h.GetEntries() < 3:
                return None
            m, s = mean_sigma(*get_hist_vecs(h, err=False))
            fit = h.Fit('gaus', 'sq0', '', m - 2 * s, m + 2 * s)
            fm, fs = fit.Parameter(1), fit.Parameter(2)
            if .8 * m < fit.Parameter(1) < 1.2 * m and s > 0 and fs < fm and fit.ParError(1) < m:  # only use gauss fit if its not deviating too much from the the mean
                current = ufloat(fm, fs + self.Precision + .03 * fm)  # add .05 as uncertainty of the device and 5% systematic error
            else:
                current = ufloat(h.GetMean(), h.GetMeanError() + .05 + .05 * h.GetMean())
        return current

    def draw_iv(self, show=True):
        g = self.make_tgrapherrors('giv', 'I-V Curve for {}'.format(self.DUTName), x=self.Data['voltages'], y=self.Data['currents'])
        format_histo(g, x_tit='Voltage [V]', y_tit='Current [nA]', y_off=1.4)
        self.draw_histo(g, show, .12)
        return g

    def set_graphs(self, averaging=1):
        self.find_data()
        sleep(.1)
        self.make_graphs(averaging)
        self.set_margins()

    def draw_indep_graphs(self, rel_time=False, ignore_jumps=True, v_range=None, f_range=None, c_range=None, averaging=1, draw_opt='ap', show=True):
        self.IgnoreJumps = ignore_jumps
        self.set_graphs(averaging)
        c = self.make_canvas('cc', 'Keithley Currents for Run {0}'.format(self.RunNumber), x=1.5, y=.75, show=show)
        self.draw_flux_pad(f_range, rel_time, draw_opt) if with_flux else self.draw_voltage_pad(v_range, draw_opt)
        self.draw_title_pad()
        self.draw_current_pad(rel_time, c_range, draw_opt)
        self.Stuff.append(c)
        run = self.Ana.RunPlan if self.IsCollection else self.RunNumber
        save_name = 'Currents{}_{}_{}'.format(self.TCString, run, self.DUTNumber)
        self.save_plots(save_name, canvas=c, ftype='png')

    def zoom_pads(self, low, high):
        self.VoltageGraph.GetXaxis().SetRangeUser(low, high)
        self.CurrentGraph.GetXaxis().SetRangeUser(low, high)

    def draw_current_pad(self, rel_t, c_range, draw_opt):
        self.draw_tpad('p3', gridx=True, margins=pad_margins, transparent=True)
        g = self.CurrentGraph
        format_histo(g, x_tit='#font[22]{Time [hh:mm]}', lab_size=label_size, x_off=1.05, tit_size=axis_title_size, t_ax_off=self.Time[0] if rel_t else 0, y_off=.55, yax_col=col_cur,
                     y_tit='#font[22]{Current [nA]}', center_y=True, x_range=[self.Time[0], self.Time[-1]], y_range=c_range, color=col_cur, markersize=marker_size)
        self.CurrentGraph.Draw(draw_opt)

    def draw_voltage_pad(self, v_range, draw_opt='ap'):
        self.draw_tpad('p1', gridy=True, margins=pad_margins, transparent=True)
        g = self.VoltageGraph
        v_range = [-1100, 1100] if v_range is None else v_range
        format_histo(g, y_range=v_range, y_tit='#font[22]{Voltage [V]}', x_range=[self.Time[0], self.Time[-1]], tit_size=axis_title_size, tick_size=0, x_off=99, l_off_x=99, center_y=True,
                     color=col_vol, y_off=title_offset, markersize=marker_size, yax_col=col_vol, lw=3, lab_size=label_size)
        g.Draw('{}y+'.format(draw_opt))

    def draw_flux_pad(self, f_range, rel_t=False, draw_opt='ap'):
        pad = self.draw_tpad('p1', margins=pad_margins, transparent=True, logy=True)
        h = self.Ana.draw_flux(rel_time=rel_t, show=False)
        pad.cd()
        f_range = [1, 20000] if f_range is None else f_range
        format_histo(h, title=' ', y_tit='#font[22]{Flux [kHz/cm^{2}]}', fill_color=4000, fill_style=4000, lw=3, y_range=f_range, stats=0, x_off=99, l_off_x=99, tick_size=0,
                     center_y=True, tit_size=axis_title_size, y_off=.7)
        h.Draw('{}y+'.format(draw_opt) if 'TGraph' in h.Class_Name() else 'histy+')

    def draw_title_pad(self):
        self.draw_tpad('p2', transparent=True)
        bias_str = 'at {b} V'.format(b=self.Bias) if self.Bias else ''
        run_str = '{n}'.format(n=self.RunNumber) if not self.IsCollection else 'Plan {rp}'.format(rp=self.Ana.RunPlan)
        text = 'Currents of {dia} {b} - Run {r} - {n}'.format(dia=self.DUTName, b=bias_str, r=run_str, n=self.Name)
        self.draw_tlatex(pad_margins[0], 1.02 - pad_margins[-1], text, align=11, size=.06)

    def find_margins(self):
        x = [min(self.Time), max(self.Time)]
        dx = .05 * (x[1] - x[0])
        y = [min(self.Currents), max(self.Currents)]
        dy = .01 * (y[1] - y[0])
        return {'x': [x[0] - dx, x[1] + dx], 'y': [y[0] - dy, y[1] + dy]}

    def set_margins(self):
        self.Margins = self.find_margins()

    def make_graphs(self, averaging=1):
        xv = array(self.Time)
        xc = array(average_list(self.Time, averaging))
        # current
        y = array(average_list(self.Currents, averaging))
        g1 = TGraph(len(xc), xc, y)
        format_histo(g1, 'Current', '', color=col_cur, markersize=.5)
        g1.SetTitle('')
        # voltage
        y = array(self.Voltages)
        g2 = TGraph(len(xv), xv, y)
        format_histo(g2, 'Voltage', '', color=col_vol, markersize=.5)
        self.CurrentGraph = g1
        self.VoltageGraph = g2

    def draw_time_axis(self, y, opt=''):
        x = self.Margins['x']
        a1 = self.draw_x_axis(y, x[0], x[1], 'Time [hh:mm]    ', off=1.2, tit_size=.05, opt=opt, lab_size=.05, tick_size=.3, l_off=.01)
        a1.SetTimeFormat("%H:%M")
        a1.SetTimeOffset(-3600)

    # endregion

    def run_exists(self, run):
        if run in self.RunLogs:
            return True
        else:
            warning('Run {run} does not exist in {tc}!'.format(run=run, tc=self.TCString))
            return False

    def print_run_times(self, run):
        run = str(run)
        if not self.run_exists(run):
            return
        log = self.RunLogs[run]
        out = '{date}: {start}-{stop}'.format(date=log['begin date'], start=log['begin time'], stop=log['stop time'])
        print(out)

    def get_time_from_log(self, t_str, year_str):
        return self.TimeZone.localize(datetime.strptime(year_str.strftime('%Y%m%d') + t_str, '%Y%m%d%H:%M:%S'))


if __name__ == '__main__':

    aparser = ArgumentParser()
    aparser.add_argument('dut', nargs='?', default=1, type=int, help='dut number [default: 1] (choose from 1,2,...)')
    aparser.add_argument('begin', nargs='?', default=12)
    aparser.add_argument('end', nargs='?', default=None)
    aparser.add_argument('-tc', '--testcampaign', nargs='?', default=None, help='YYYYMM beam test [default in main.ini]')
    aparser.add_argument('-v', '--verbose', action='store_false')
    aparser.add_argument('-c', '--collection', action='store_true', help='begin analysis collection')
    pargs = aparser.parse_args()

    z = Currents(test_campaign=pargs.testcampaign, dut=pargs.dut, begin=pargs.begin, end=pargs.end if not pargs.collection else False, verbose=pargs.verbose)
