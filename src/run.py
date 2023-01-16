#!/usr/bin/env python
# --------------------------------------------------------
#       run module that handles the files and information of a single run
# created on October 5th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from utility.utils import print_table, datetime, ev2str, remove_letters, Dir, array, small_banner, isint
from plotting.utils import load_json, warning, critical
import plotting.latex as tex
from src.analysis import Analysis, Path, choose
from src.dut import DUT


def load_runlog(p: Path):
    f = p.joinpath(Analysis.Config.get('data', 'runlog file'))
    if not f.exists():
        warning('runlog file does not exist! -> creating new one!')
        from src.spreadsheet import make
        make(p.stem.replace('-', ''))
    return load_json(f)


def load_nrs(p: Path):
    log = load_runlog(p)
    return [key for key, dic in log.items() if 'status' not in dic or dic['status'] == 'green']


def init_batch(name, dut_nr, beam_test, log=None):
    custom = Batch.load_custom(beam_test.Name)
    return DUTBatch(name, beam_test.Path, log) if name in custom and type(custom[name]) == dict else Batch(name, dut_nr, beam_test.Path, log)


class Batch:
    """ class containing the run infos of a single batch. """

    def __init__(self, name, dut_nr, data_dir: Path, log=None):
        self.Name = name
        self.DataDir = data_dir
        self.Log = load_runlog(self.DataDir) if log is None else log
        self.LogNames = self.load_log_names()
        self.Custom = Batch.load_custom(data_dir)
        self.Runs = self.load_runs(dut_nr)
        self.FirstRun = self.Runs[0]
        self.RunNrs = array([run.Number for run in self.Runs])
        self.Size = len(self.Runs)
        self.DUT = self.Runs[0].DUT

        self.FileName = data_dir.joinpath('data', f'batch{self}.hdf5')

        if not self.verify():
            critical('the duts of the individual runs do not agree! This is not implemented yet ...')

    def __str__(self):
        return str(self.Name)

    def __repr__(self):
        return f'Batch {self} ({self.Size} runs {self.min_run}-{self.max_run}, {ev2str(self.n_ev)} ev)'

    def __getitem__(self, item):
        return self.Runs[item]

    @classmethod
    def from_beamtest(cls, bt, name, dut_nr, log=None):
        return cls(name, dut_nr, bt.Path, log)

    def verify(self):
        if self.Name is None:
            return True
        duts = [run.DUTs for run in self.Runs]
        return all([d == duts[0] for d in duts])

    def load_log_names(self):
        return sorted(list(set([dic['batch'] for dic in self.Log.values() if dic['batch']])), key=lambda x: (int(remove_letters(x)), x))

    @staticmethod
    def load_custom(name):
        """ user defined batches in ./ensembles"""
        return load_json(Dir.joinpath('ensembles', f'{Path(name).stem}.json'))

    def run_is_good(self, nr):
        return self.Log[str(nr)]['status'] == 'green'

    def load_runs(self, dut_nr):
        if self.Name in self.LogNames or self.Name is None:
            is_good = lambda dic: (self.Name is None or dic['batch'] == self.Name) and dic['status'] == 'green'
            return [Run(key, dut_nr, self.DataDir, log=self.Log) for key, dic in self.Log.items() if is_good(dic)]
        if self.Name in self.Custom:
            return [Run(nr, dut_nr, self.DataDir, log=self.Log) for nr in self.Custom[self.Name] if self.run_is_good(nr)]
        critical('unknown batch name')

    @property
    def n_ev(self):
        return sum(run.n_ev for run in self.Runs)

    @property
    def min_run(self):
        return self[0]

    @property
    def max_run(self):
        return self[-1]

    def show(self):
        small_banner(f'{self!r}')
        print_table([run.info for run in self.Runs], header=['Nr.', 'Events', 'DUTs', 'Begin', 'End'])

    def show_all(self):
        data = {n: [] for n in self.LogNames}
        for run, log in self.Log.items():
            if log['status'] == 'green':
                log['run_nr'] = run
                data[log['batch']].append(log)
        r_str = lambda x: f'{x[0]["run_nr"]:>03}-{x[-1]["run_nr"]:>03}'
        t_str = lambda x: [f'{datetime.fromtimestamp(t)}'[-8:-3] for t in [x[0]['start'], x[-1]['end']]]
        rows = [[n, r_str(logs), ev2str(sum([log['events'] for log in logs])), ', '.join(logs[0]['duts'])] + t_str(logs) for n, logs in data.items() if len(logs)]
        print_table(rows, header=['Batch', 'Events', 'Runs', 'DUTs', 'Begin', 'End'])

    def show_custom(self, dut_nr=0, dut=None):
        rows = []
        batches = [Batch(n, dut_nr, self.DataDir, self.Log) for n, runs in self.Custom.items() if len(self.Log[str(runs[0])]['duts']) > dut_nr]
        batches = list(filter(lambda x: x.DUT.Name == dut, batches)) if dut is not None else batches
        irrs = array([i.DUT.get_irradiation(self.DataDir.name.replace('-', '')) for i in batches]).astype('d')
        for i, b_ in enumerate(batches):
            tc = self.DataDir.name.replace('-', '')
            bt = tex.multirow(datetime.strptime(tc, '%Y%m').strftime('%b%y'), len(batches)) if not i else ''
            irr = f'{irrs[i]:.1f}' if not all(irrs == irrs[0]) else tex.multirow(f'{irrs[i]:.1f}', len(batches)) if not i else ''
            rows.append([bt, irr, b_.Name, f'{tex.num_range(*b_.RunNrs[[0, -1]])} ({b_.Size})', 'data', b_.DUT.Bias, ev2str(b_.n_ev)])
        print(tex.table(None, rows=rows))

    def find_runs(self, dut, bias, min_run=0):
        return array([run.Number for run in self.Runs if dut in run.DUTs and int(run.Info['hv'][run.DUTs.index(dut)]) == bias and run >= min_run])  # noqa

    @staticmethod
    def find_dut_numbers(batch_name, dut_name, log, bt):
        new_batches = load_json(Dir.joinpath('ensembles', f'{bt.Path.stem}.json'))
        logs = [log[str(run)] for run in new_batches[batch_name]] if batch_name in new_batches else filter(lambda x: x['batch'] == batch_name and x['status'] == 'green', log.values())
        return [d['duts'].index(dut_name) for d in logs]


class DUTBatch(Batch):
    """ extension of batch class for a single DUT (for runs with mismatching duts). """

    def __init__(self, name, data_dir: Path, log=None):
        super().__init__(name, 0, data_dir, log)

    def verify(self):
        return True

    def load_runs(self, dut_nr=None):
        return [Run.from_dut_name(nr, dut_name, self.DataDir, self.Log) for dut_name, nrs in self.Custom[self.Name].items() for nr in nrs if self.run_is_good(nr)]


class Run:
    """ Run class containing all the information for a single run from the tree and the json file. """

    def __init__(self, run_number, dut_number, tc_dir: Path, log=None):

        # Main Info
        self.Number = int(run_number)
        self.TCDir = tc_dir

        # Files
        self.FileName = self.TCDir.joinpath('data', f'run{self.Number:04d}.hdf5')

        # Info
        self.Info = self.load_info(log)
        self.DUTs = self.Info['duts']
        self.NDUTs = len(self.DUTs)
        self.Positions = self.Info['dut position']  # which of the DUT slots are occupied
        self.DUT = DUT(dut_number, self.Info, has_ref=tc_dir.stem in Analysis.Config.get('REF', 'dates'))

        # Times for the run log
        self.LogStart = self.Info['start']
        self.LogEnd = self.Info['end']

    def __str__(self):
        return str(self.Number)

    def __repr__(self):
        return f'Run {self} ({ev2str(self.n_ev)} ev)'

    def __format__(self, format_spec):
        return f'{self.Number:{format_spec}}'

    def __le__(self, other):
        return self.Number <= (other.Number if isinstance(other, Run) else other)

    def __lt__(self, other):
        return self.Number < (other.Number if isinstance(other, Run) else other)

    def __ge__(self, other):
        return self.Number >= (other.Number if isinstance(other, Run) else other)

    def __gt__(self, other):
        return self.Number > (other.Number if isinstance(other, Run) else other)

    @property
    def info(self):
        return [str(self.Number), ev2str(self.n_ev), ', '.join(self.DUTs)] + [f'{datetime.fromtimestamp(t)}'[-8:-3] for t in [self.LogStart, self.LogEnd]]

    def load_info(self, log=None) -> dict:
        return (load_runlog(self.TCDir) if log is None else log)[str(self)]

    def print_info(self):
        print(f'{self!r}')
        print_table(rows=[[key, str(datetime.fromtimestamp(value) if key in ['start', 'end'] else value)] for key, value in self.Info.items()])

    @classmethod
    def from_ana(cls, run_number, dut=0, ana: Analysis = None):
        ana = choose(ana, Analysis)
        return cls(run_number, dut, ana.BeamTest.Path)

    @classmethod
    def from_dut_name(cls, run_number, dut_name, tc_dir, log):
        return cls(run_number, log[str(run_number)]['duts'].index(dut_name), tc_dir, log)

    @property
    def n_ev(self):
        return self.Info['events']


class Ensemble(object):
    """ General enseble class for runs. """

    FilePath = Dir.joinpath('ensembles', 'scans.json')

    def __init__(self, name: str):

        self.Name = name
        self.Dic = self.load_dic()
        self.DUTName = self.Dic.pop('dut')
        self.BeamTests = {bt: Analysis.load_test_campaign(bt) for bt in self.Dic}
        self.Logs = {bt: load_runlog(Analysis.load_test_campaign(bt).Path) for bt in self.Dic}
        self.Units = self.init_units()
        self.Size = len(self.Units)
        self.DUT = self.Units[0].DUT

    def __getitem__(self, item):
        return self.Units[item]

    def __str__(self):
        return self.Name

    def __repr__(self):
        return f'{self.__class__.__name__} {self.Name} with {self.Size} units'

    def load_dic(self):
        d = load_json(Ensemble.FilePath)
        if self.Name not in d:
            critical(f'could not find ensemble "{self.Name}" in {self.FilePath}')
        return d[self.Name]

    def init_units(self):
        return [self.init_run(name, bt) if isint(name) else self.init_batch(name, bt) for bt in self.Dic for name in self.Dic[bt]]

    def init_run(self, run_number, bt):
        return Run(run_number, self.Logs[bt][str(run_number)]['duts'].index(self.DUTName), self.BeamTests[bt].Path, self.Logs[bt])

    def init_batch(self, name, bt):
        dut_nr = Batch.find_dut_numbers(name, self.DUTName, self.Logs[bt], self.BeamTests[bt])[0]  # TODO: update if batches with varying dut_nrs are implemented
        return Batch(name, dut_nr, self.BeamTests[bt].Path, self.Logs[bt])

    @property
    def biases(self):
        return [u.DUT.Bias for u in self.Units]


if __name__ == '__main__':

    from argparse import ArgumentParser
    p_ = ArgumentParser()
    p_.add_argument('b', nargs='?', default='16b')
    p_.add_argument('-a', action='store_true')
    p_.add_argument('-all', action='store_true')
    args = p_.parse_args()

    a = Analysis()
    b = init_batch(None if args.all else args.b, 0, a.BeamTest)
    r = b[0]
    if args.a:
        b.show_all()
    else:
        b.show()
