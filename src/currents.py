from os.path import getsize

from numpy import genfromtxt, isnan, datetime64, invert, char, uint32
from pytz import timezone

import src.bins as bins
from src.analysis import *
from src.utils import *
from glob import glob


class Currents(Analysis):
    """reads in information from the keithley log file"""

    def __init__(self, analysis=None, test_campaign=None, dut=None, begin=None, end=None, averaging=None, verbose=False):
        Analysis.__init__(self, test_campaign if analysis is None else analysis.BeamTest.Tag, verbose=verbose)

        # Settings
        self.Averaging = averaging
        self.TimeZone = timezone('Europe/Zurich')
        self.DataDir = self.BeamTest.Path.joinpath('hvdata')

        # Config
        self.Ana = analysis
        self.IsCollection = hasattr(analysis, 'Runs')
        self.Collection = None
        self.RunPlan = self.load_run_plan()  # required for plotting
        self.RunLogs = self.Ana.Run.Logs
        self.Run = self.load_run()
        self.HVConfig = Config(self.DataDir.joinpath('config.ini'))
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
        if self.Ana is not None and data.size:
            data['timestamps'] -= uint32(data['timestamps'][0] - self.Run.StartTime)  # synchronise time vectors
        return data

    def reload_data(self, ignore_jumps):
        if ignore_jumps != self.IgnoreJumps:
            self.IgnoreJumps = ignore_jumps
            self.Data = self.load_data()

    def load_bias(self):
        return self.Run.DUT.Bias if hasattr(self.Run, 'Bias') else None

    def load_run(self):
        return None if self.Ana is None else self.Ana.Run if not self.IsCollection else self.Ana.RunPlan

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
                return (self.TimeZone.localize(datetime.strptime('{}-{}'.format(self.BeamTest.year, t), '%Y-%m/%d-%H:%M:%S')) for t in [begin, end])
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

    def get_title(self):
        bias_str = 'at {b} V'.format(b=self.Bias) if self.Bias else ''
        run_str = '{n}'.format(n=self.Run.Number) if not self.IsCollection else 'Plan {rp}'.format(rp=self.Ana.RunPlan)
        return 'Currents of {dia} {b} - Run {r} - {n}'.format(dia=self.DUTName, b=bias_str, r=run_str, n=self.Name)
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
    def draw_profile(self, bw=None, **dkw):
        x, y = self.Data['timestamps'], self.Data['currents']
        return self.Draw.profile(x, y, find_bins(x, w=bw), 'Leakage Current', **prep_kw(dkw, x_tit='Time [hh:mm]', y_tit='Current [nA]', t_ax_off=0, markersize=.7, **Draw.mode(2)))

    def draw_distribution(self, show=True):
        m, s = mean_sigma(self.Data['currents'])
        xmin, xmax = m - 4 * max(s, .1), m + 4 * max(s, .1)
        return self.Draw.distribution(self.Data['currents'], 'Current Distribution', bins.make(xmin, xmax, self.Precision * 2), show=show, x_tit='Current [nA]')

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

    def draw_iv(self, **dkw):
        return self.Draw.graph(self.Data['voltages'], self.Data['currents'], f'I-V Curve for {self.DUTName}', **prep_kw(dkw, x_tit='Voltage [V]', y_tit='Current [nA]'))

    def draw(self, rel_time=False, ignore_jumps=True, v_range=None, c_range=None, averaging=1, **dkw):
        self.reload_data(ignore_jumps)
        t, c, v = (average_list(self.Data[n], averaging) for n in ['timestamps', 'currents', 'voltages'])
        gv = self.Draw.graph(t, v, self.get_title(), y_tit='Voltage [nA]', yax_col=602, color=602, y_range=choose(v_range, [-100, 0]), l_off_x=10, x_ticks=0, show=False)
        gc = self.Draw.graph(t, c, x_tit='Time [hh:mm]', y_tit='Current [nA]', yax_col=899, color=899, y_range=choose(c_range, [round_down_to(min(c)), round_up_to(max(c))]), show=False)
        for g in [gc, gv]:
            format_histo(g, lab_size=.05, x_off=1.05, tit_size=.06, t_ax_off=t[0] if rel_time else 0, y_off=.8, center_y=True, x_range=[t[0], t[-1]], markersize=.3)
        self.Draw(gv, **prep_kw(dkw, **Draw.mode(2, lm=.1, rm=.1), draw_opt='aly+'))
        self.Draw.tpad('pc', transparent=True, c=get_last_canvas())
        gc.Draw(dkw['draw_opt'] if 'draw_opt' in dkw else 'al')
        update_canvas()
    # endregion PLOTTING
    # ----------------------------------------

    def run_exists(self, run):
        if run in self.RunLogs:
            return True
        else:
            warning(f'Run {run} does not exist in {self.BeamTest}!')
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

    from argparse import ArgumentParser

    aparser = ArgumentParser()
    aparser.add_argument('dut', nargs='?', default=1, type=int, help='dut number [default: 1] (choose from 1,2,...)')
    aparser.add_argument('begin', nargs='?', default=12)
    aparser.add_argument('end', nargs='?', default=None)
    aparser.add_argument('-tc', '--testcampaign', nargs='?', default=None, help='YYYYMM beam test [default in main.ini]')
    aparser.add_argument('-v', '--verbose', action='store_false')
    aparser.add_argument('-c', '--collection', action='store_true', help='begin analysis collection')
    pargs = aparser.parse_args()

    z = Currents(test_campaign=pargs.testcampaign, dut=pargs.dut, begin=pargs.begin, end=pargs.end if not pargs.collection else False, verbose=pargs.verbose)
