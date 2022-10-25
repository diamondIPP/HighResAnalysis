#!/usr/bin/env python
# --------------------------------------------------------
#       Script to automatically convert several files
# created on October 19th 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------


from src.run import Batch
from plotting.utils import choose, info, colored, GREEN, RED
from utility.utils import print_banner, PBAR, small_banner
from src.analysis import Analysis
from cern.converter import CERNConverter, Converter
from multiprocessing import cpu_count, Pool


class AutoConvert:

    def __init__(self, first_run=None, last_run=None, batch=None, beamtest=None, verbose=False, force=False):

        self.Ana = Analysis(beamtest, verbose=verbose)
        self.Batch = Batch(batch, dut_nr=0, data_dir=self.Ana.BeamTest.Path)
        self.Converter = CERNConverter if self.Ana.BeamTest.Location == 'CERN' else Converter

        self.FirstRun, self.LastRun = first_run, last_run
        self.Force = force

    def __repr__(self):
        return f'Converter for {self.Batch!r}'

    @property
    def runs(self):
        rmin, rmax = self.first_run, self.last_run
        return [r for r in self.Batch.Runs if (not r.FileName.exists() or self.Force) and rmin <= r <= rmax]

    @property
    def converters(self):
        return [self.Converter.from_run(r) for r in self.runs]

    @property
    def first_run(self):
        return choose(self.FirstRun, self.Batch.min_run)

    @property
    def last_run(self):
        return choose(self.LastRun, self.Batch.max_run)

    def copy_raw_files(self):
        conv = self.converters
        n_missing = sum([not f.exists() for c in conv for f in c.raw_files])
        if n_missing > 0:
            info(f'try to copy {n_missing} raw files from {Analysis.Config.get("data", "server")} ...')
            PBAR.start(n_missing)
            for c in conv:
                for f in c.raw_files:
                    if c.download_raw_file(f, out=False) is not None:
                        PBAR.update()

    def run(self):
        """parallel conversion"""
        self.copy_raw_files()
        info(f'Creating pool with {cpu_count()} processes')
        with Pool() as pool:
            result = pool.map_async(convert_run, self.converters)
            conv = result.get()
            small_banner('Summary:')
            for c in conv:
                speed = f'{c.Run.n_ev}, {c.Run.n_ev / c.T1.total_seconds():1.0f} Events/s' if c.finished else 'NOT CONVERTED!'
                print(colored(f'{c} --> {c.T1} ({speed})'), GREEN if c.finished else RED)


def convert_run(c: Converter):
    return c


from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-m', action='store_true', help='turn parallel processing ON')
parser.add_argument('-tc', nargs='?', default=None)
parser.add_argument('s', nargs='?', default=None, help='run number where to start, default [None], = stop if no end is provided', type=int)
parser.add_argument('e', nargs='?', default=None, help='run number where to stop, default [None]')
parser.add_argument('-b', nargs='?', default=None, help='batch number, default [None]')
parser.add_argument('-v', action='store_false', help='turn verbose OFF')
parser.add_argument('-t', action='store_true', help='turn test mode ON')
parser.add_argument('-f', action='store_true', help='force conversion')

args = parser.parse_args()

z = AutoConvert(args.s, args.e, args.b, args.tc, args.v, args.f)
if not args.t:
    runs = z.runs
    if len(runs):
        print_banner(f'Start converting runs {runs[0]} - {runs[-1]}', color=GREEN)
        z.run()
        print_banner('Finished Conversion!', color=GREEN)
    else:
        info('There is nothing to convert :-)\n', blank_lines=1)
