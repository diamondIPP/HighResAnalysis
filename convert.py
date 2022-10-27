#!/usr/bin/env python
# --------------------------------------------------------
#       Script to automatically convert several files
# created on October 19th 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------


from src.run import Batch
from plotting.utils import choose, info, colored, GREEN, RED, check_call
from utility.utils import print_banner, PBAR, small_banner
from src.analysis import Analysis
from cern.converter import CERNConverter, Converter
from src.converter import batch_converter
from multiprocessing import cpu_count, Pool


class AutoConvert:

    Force = False

    def __init__(self, first_run=None, last_run=None, batch=None, beamtest=None, verbose=False, force=False):

        self.Ana = Analysis(beamtest, verbose=verbose)
        self.Batch = Batch(batch, dut_nr=0, data_dir=self.Ana.BeamTest.Path)
        self.Converter = CERNConverter if self.Ana.BeamTest.Location == 'CERN' else Converter

        self.FirstRun, self.LastRun = first_run, last_run
        AutoConvert.Force = force

    def __repr__(self):
        return f'Converter for {self.Batch!r}'

    @property
    def converters(self):
        return [self.Converter.from_run(r) for r in self.Batch.Runs if (not r.FileName.exists() or self.Force) and self.first_run <= r <= self.last_run]

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

    @staticmethod
    def convert_run(c: Converter):
        c.run(force=AutoConvert.Force)
        return c

    def run(self):
        """parallel conversion"""
        self.copy_raw_files()
        info(f'Creating pool with {cpu_count()} processes')
        with Pool() as pool:
            result = pool.map_async(self.convert_run, self.converters)
            conv = result.get()
            small_banner('Summary:')
            for c in conv:
                speed = f'{c.Run.n_ev}, {c.Run.n_ev / c.T1.total_seconds():1.0f} Events/s' if c.finished else 'NOT CONVERTED!'
                print(colored(f'{c} --> {c.T1} ({speed})'), GREEN if c.finished else RED)


class BatchConvert(AutoConvert):

    def __init__(self, batch=None, beamtest=None, verbose=False, force=False):
        super().__init__(None, None, batch, beamtest, verbose, force)
        self.Converters = [self.Converter.from_run(r) for r in self.Batch.Runs]
        self.Converter = batch_converter(self.Converters[0].__class__).from_batch(self.Batch)

    @staticmethod
    def convert_run(c: Converter):
        c.raw2root()
        return c

    @property
    def converters(self):
        return list(filter(lambda c: not c.proteus_raw_file_path().exists() or not c.trigger_info_file().exists(), self.Converters))

    def run(self):
        super().run()  # create root files of the single runs
        self.merge_files()
        self.Converter.run(force=self.Force)  # tracking and hdf5 conversion of the single merged file

    def merge_files(self):
        """merge proteus raw files of all batch runs"""
        cmd = f'hadd -f {self.Converter.proteus_raw_file_path()} {" ".join(str(c_.proteus_raw_file_path()) for c_ in self.Converters)}'
        info(cmd)
        check_call(cmd, shell=True)
        if self.Converter.trigger_info_file() != self.Converter.proteus_raw_file_path():
            cmd = f'hadd -f {self.Converter.trigger_info_file()} {" ".join(str(c_.trigger_info_file()) for c_ in self.Converters)}'
            info(cmd)
            check_call(cmd, shell=True)


from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument('-m', action='store_true', help='turn parallel processing ON')
parser.add_argument('-tc', nargs='?', default=None)
parser.add_argument('s', nargs='?', default=None, help='run number where to start, default [None], = stop if no end is provided', type=int)
parser.add_argument('e', nargs='?', default=None, help='run number where to stop, default [None]', type=int)
parser.add_argument('-b', nargs='?', default='23b', help='batch number, default [None]')
parser.add_argument('-v', action='store_false', help='turn verbose OFF')
parser.add_argument('-t', action='store_true', help='turn test mode ON')
parser.add_argument('-f', action='store_true', help='force conversion')

args = parser.parse_args()

z = AutoConvert(args.s, args.e, args.b, args.tc, args.v, args.f) if args.s is not None else BatchConvert(args.b, args.tc, args.v, args.f)
a = z.Converter
cs = z.Converters

if not args.t:
    cs = z.converters
    if len(cs):
        print_banner(f'Start converting runs {cs[0].Run} - {cs[-1].Run}', color=GREEN)
        z.run()
        print_banner('Finished Conversion!', color=GREEN)
    else:
        info('There is nothing to convert :-)\n', blank_lines=1)
