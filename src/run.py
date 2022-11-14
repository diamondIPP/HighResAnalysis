#!/usr/bin/env python
# --------------------------------------------------------
#       run module that handles the files and information of a single run
# created on October 5th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from utility.utils import print_table, datetime, ev2str, remove_letters, Dir, array
from plotting.utils import load_json, warning, critical
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
    return [key for key, dic in log.items() if dic['status'] == 'green']


class Batch:
    """ class containing the run infos of a single batch. """

    def __init__(self, name, dut_nr, data_dir: Path):
        self.Name = name
        self.DataDir = data_dir
        self.Log = load_runlog(self.DataDir)
        self.LogNames = self.load_log_names()
        self.Runs = self.load_runs(dut_nr)
        self.RunNrs = [run.Number for run in self.Runs]
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

    def verify(self):
        duts = array([run.DUTs for run in self.Runs])
        return (duts == duts[0]).all()

    def load_log_names(self):
        return sorted(list(set([dic['batch'] for dic in self.Log.values()])), key=lambda x: (int(remove_letters(x)), x))

    def load_runs(self, dut_nr):
        if self.Name in self.LogNames:
            is_good = lambda dic: (self.Name is None or dic['batch'] == self.Name) and dic['status'] == 'green'
            return [Run(key, dut_nr, self.DataDir, log=self.Log) for key, dic in self.Log.items() if is_good(dic)]
        dic = load_json(Dir.joinpath('config', 'batches.json'))
        if self.Name in dic:
            return [Run(nr, dut_nr, self.DataDir, log=self.Log) for nr in dic[self.Name] if self.Log[str(nr)]['status'] == 'green']
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
        self.DUT = DUT(dut_number, self.Info)

        # Times  TODO: Fix later with real timestamps from the data
        self.StartTime = self.Info['start']
        self.EndTime = self.Info['end']

    def __str__(self):
        return str(self.Number)

    def __repr__(self):
        return f'Run {self} ({ev2str(self.n_ev)} ev)'

    def __format__(self, format_spec):
        return f'{self.Number:{format_spec}}'

    def __le__(self, other):
        return self.Number <= (other.Number if isinstance(other, Run) else other)

    def __ge__(self, other):
        return self.Number >= (other.Number if isinstance(other, Run) else other)

    def load_info(self, log=None) -> dict:
        return (load_runlog(self.TCDir) if log is None else log)[str(self)]

    def print_info(self):
        print(f'{self!r}')
        print_table(rows=[[key, str(datetime.fromtimestamp(value) if key in ['start', 'end'] else value)] for key, value in self.Info.items()])

    @classmethod
    def from_ana(cls, run_number, dut=0, ana: Analysis = None):
        ana = choose(ana, Analysis)
        return cls(run_number, dut, ana.BeamTest.Path)

    @property
    def n_ev(self):
        return self.Info['events']


if __name__ == '__main__':
    a = Analysis()
    b = Batch('23b', 0, a.BeamTest.Path)
    r = b[0]
