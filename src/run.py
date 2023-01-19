#!/usr/bin/env python
# --------------------------------------------------------
#       run module that handles the files and information of a single run
# created on October 5th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from numpy import where, roll, split
from utility.utils import print_table, datetime, ev2str, remove_letters, Dir, array, small_banner, isint, do_pickle, file_hash
from plotting.utils import load_json, warning, critical
from src.analysis import Analysis, Path, choose, BeamTest
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


def init_batch(name, dut, beam_test: BeamTest, log=None):
    batches = load_batches(beam_test)
    return DUTBatch(name, beam_test, log) if type(batches[name]) == dict else Batch(name, dut, beam_test, log) if isint(dut) else Batch.from_dut_name(name, dut, beam_test, log)


def load_batches(bt: BeamTest, redo=False):
    """ unify batches from run logs and custom batches and save them in a tmp file. """
    log_file = bt.Path.joinpath(Analysis.Config.get('data', 'runlog file'))
    custom_file = Dir.joinpath('ensembles', f'{bt.Name}.json')
    tmp_file = Analysis.MetaDir.joinpath('ensembles', f'{bt.Tag}.pickle')
    tmp_file.parent.mkdir(exist_ok=True)


    if tmp_file.exists() and not redo:
        h1, h2, bs = do_pickle(tmp_file)
        if file_hash(log_file) == h1 and file_hash(custom_file) == h2:
            return bs

    def load():
        log = load_json(log_file)

        # batches from log
        d = array(sorted([(run, dic['batch']) for run, dic in log.items() if dic['status'] == 'green'], key=lambda x: x[1])).T
        s = where(d[1] != roll(d[1], 1))[0]  # indices where the new batches start
        batches = {batch: runs.astype('i8') for batch, runs in zip(d[1][s], split(d[0], s[1:]))}

        # custom batches
        if custom_file.exists():
            good_runs = lambda runs: array([run for run in runs if log[str(run)]['status'] == 'green'])
            batches.update({batch: {dut: good_runs(d_runs) for dut, d_runs in runs.items()} if type(runs) is dict else good_runs(runs) for batch, runs in load_json(custom_file).items()})

        return file_hash(log_file), file_hash(custom_file), dict(sorted(batches.items(), key=lambda dic_pair: int(remove_letters(dic_pair[0]))))

    return do_pickle(tmp_file, load, redo=True)[2]



class Batch:
    """ class containing the run infos of a single batch. """

    DUTName = None

    def __init__(self, name, dut_nr, beam_test: BeamTest, log=None):
        self.Name = name
        self.BeamTest = beam_test
        self.DataDir = beam_test.Path
        self.Log = load_runlog(self.DataDir) if log is None else log

        self.Runs = self.load_runs(dut_nr)
        self.FirstRun = self.Runs[0]
        self.RunNrs = array([run.Number for run in self.Runs])
        self.Size = len(self.Runs)
        self.DUT = self.Runs[0].DUT

        self.FileName = self.DataDir.joinpath('data', f'batch{self}.hdf5')

        if not self.verify():
            critical('the duts of the individual runs do not agree! This is not implemented yet ...')

    def __str__(self):
        return str(self.Name)

    def __repr__(self):
        return f'{self.__class__.__name__} {self} ({self.Size} runs {self.min_run}-{self.max_run}, {ev2str(self.n_ev)} ev)'

    def __getitem__(self, item):
        return self.Runs[item]

    @classmethod
    def from_dut_name(cls, name, dut_name, bt: BeamTest, log=None):
        log = load_runlog(bt.Path) if log is None else log
        dut_nrs = list(set([log[str(nr)]['duts'].index(dut_name) for nr in load_batches(bt)[name]]))
        if len(dut_nrs) == 1:
            return cls(name, dut_nrs[0], bt, log)
        critical(f'cannot instanciate Batch {name} from DUT {dut_name}, varying DUT numbers ...')

    def verify(self):
        if self.Name is None:
            return True
        duts = [run.DUTs for run in self.Runs]
        return all([d == duts[0] for d in duts])

    def load_runs(self, dut_nr):
        batches = load_batches(self.BeamTest)
        if self.Name in batches:
            return [Run(nr, dut_nr, self.DataDir, log=self.Log) for nr in batches[self.Name]]
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
        rows = []
        log = {int(i): dic for i, dic in self.Log.items()}
        for n, rs in load_batches(self.BeamTest).items():
            duts, rs = (rs.keys(), list(rs.values())[0]) if type(rs) is dict else (log[rs[0]]['duts'], rs)
            t_str = [f'{datetime.fromtimestamp(t)}'[-8:-3] for t in [log[rs[0]]['start'], log[rs[-1]]['end']]]
            rows.append([n, f'{rs[0]:03d}-{rs[-1]:03d}', ev2str(sum([log[i]['events'] for i in rs])), ', '.join(duts)] + t_str)
        print_table(rows, header=['Batch', 'Runs', 'Events', 'DUTs', 'Begin', 'End'])

    def find_runs(self, dut, bias, min_run=0):
        return array([run.Number for run in self.Runs if dut in run.DUTs and int(run.Info['hv'][run.DUTs.index(dut)]) == bias and run >= min_run])  # noqa

    @staticmethod
    def find_dut_numbers(batch_name, dut_name, log, bt):
        new_batches = load_json(Dir.joinpath('ensembles', f'{bt.Path.stem}.json'))
        logs = [log[str(run)] for run in new_batches[batch_name]] if batch_name in new_batches else filter(lambda x: x['batch'] == batch_name and x['status'] == 'green', log.values())
        return [d['duts'].index(dut_name) for d in logs]


class DUTBatch(Batch):
    """ extension of batch class for a single DUT (for runs with mismatching duts). """

    def __init__(self, name, beam_test: BeamTest, log=None):
        super().__init__(name, 0, beam_test, log)
        self.DUTName = self.DUT.Name

    def verify(self):
        return True

    def load_runs(self, dut_nr=None):
        batches = load_batches(self.BeamTest)
        return [Run.from_dut_name(nr, dut_name, self.DataDir, self.Log) for dut_name, nrs in batches[self.Name].items() for nr in nrs]


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
    """ General enseble class for runs and batches. """

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
        return [self.init_run(n, bt) if isint(n) else self.init_batch(n, bt) for bt in self.Dic for n in self.Dic[bt]]

    def init_run(self, run_number, bt):
        return Run(run_number, self.Logs[bt][str(run_number)]['duts'].index(self.DUTName), self.BeamTests[bt].Path, self.Logs[bt])

    def init_batch(self, name, bt):
        return init_batch(name, self.DUTName, self.BeamTests[bt], self.Logs[bt])

    @property
    def biases(self):
        return [u.DUT.Bias for u in self.Units]


if __name__ == '__main__':

    from argparse import ArgumentParser
    p_ = ArgumentParser()
    p_.add_argument('b', nargs='?', default='1a')
    p_.add_argument('-a', action='store_true')
    p_.add_argument('-all', action='store_true')
    args = p_.parse_args()

    beam_tests = [Analysis.load_test_campaign(str(bt)) for ls in Analysis.Locations.values() for bt in ls]

    a = Analysis()
    b = init_batch(None if args.all else args.b, 0, a.BeamTest)
    r = b[0]
    if args.a:
        b.show_all()
    else:
        b.show()
