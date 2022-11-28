#!/usr/bin/env python
# --------------------------------------------------------
#       Script to automatically convert several files
# created on October 19th 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
import h5py
import uproot

from src.run import Batch
from plotting.utils import choose, info, colored, GREEN, RED, check_call, critical, warning
from utility.utils import print_banner, PBAR, small_banner, byte2str
from src.analysis import Analysis
from cern.converter import CERNConverter, Converter
from src.converter import batch_converter
from multiprocessing import cpu_count, Pool
from numpy import where, diff, array, concatenate, cumsum, append, mean
from functools import partial


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
    def final_files_exist(self):
        return all([r.FileName.exists() and self.first_run <= r <= self.last_run for r in self.Batch.Runs])

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
        conv = self.converters
        if not conv:
            return True
        self.copy_raw_files()
        info(f'Creating pool with {cpu_count()} processes')
        with Pool() as pool:
            result = pool.map_async(self.convert_run, conv)
            conv = result.get()
            small_banner('Summary:')
            for c in conv:
                speed = f'{c.Run.n_ev}, {c.Run.n_ev / c.T1.total_seconds():1.0f} Events/s' if c.finished else 'NOT CONVERTED!'
                print(colored(f'{c} --> {str(c.T1)[:-5]} ({speed})', GREEN if c.finished else RED))
            return all([c.finished for c in conv])


class BatchConvert(AutoConvert):

    def __init__(self, batch=None, beamtest=None, verbose=False, force=False):
        super().__init__(None, None, batch, beamtest, verbose, force)
        self.Converters = [self.Converter.from_run(r) for r in self.Batch.Runs]
        self.Converter = batch_converter(self.Converters[0].__class__).from_batch(self.Batch)

    @staticmethod
    def convert_run(c: Converter):
        f, out = c.Proteus.Steps[-1]
        c.run(force=True, steps=c.first_steps + [(partial(f, progress=False), out)], rm=False)
        return c

    @property
    def converters(self):
        return list(filter(lambda c: not c.proteus_raw_file_path().exists() or not c.Proteus.OutFilePath.exists(), self.Converters))

    def single(self):
        return super().run()

    def run(self):
        if super().run():  # create root files of the single runs
            self.merge_files()
            self.Converter.run(force=self.Force, rm=False)  # tracking and hdf5 conversion of the single merged file
            self.fix_event_nrs()
            self.remove_aux_files()
        else:
            critical('Not all runs were converted ... please re-run')

    @property
    def file_size_ratios(self):
        return array([c.Proteus.OutFilePath.stat().st_size / c.Run.n_ev for c in self.Converters])

    def check_file_sizes(self):
        r = self.file_size_ratios
        bad_runs = array(self.Batch.Runs)[r < mean(r) * .9]
        return True if bad_runs.size == 0 else bad_runs

    def check_proteus_files(self, pl=0, branch='trk_v', min_run=0):
        conv = [c for c in self.Converters if c.Run >= min_run]
        PBAR.start(len(conv))
        for c in conv:
            with uproot.open(c.Proteus.OutFilePath) as f:
                key = f.keys(recursive=False)[pl]
                try:
                    array(f[f'{key}/tracks_clusters_matched'][branch])
                except Exception as err:
                    warning(f'{c} has corrupted data ({err})')
                PBAR.update()

    def merge_proteus_files(self):
        """merge proteus out files of all batch runs"""
        out = self.Converter.Proteus.OutFilePath
        if not out.exists() or self.Force:
            cmd = f'hadd -f {out} {" ".join(str(c_.Proteus.OutFilePath) for c_ in self.Converters)}'
            info(cmd)
            check_call(cmd, shell=True)
        else:
            info(f'found {out}')

    def merge_trigger_info_files(self):
        """merge dut files with trigger info for the CERN analysis"""
        out = self.Converter.trigger_info_file()
        if not out.exists() or self.Force:
            cmd = f'hadd -f {out} {" ".join(str(c_.trigger_info_file()) for c_ in self.Converters)}'
            info(cmd)
            check_call(cmd, shell=True)
        else:
            info(f'found {out}')

    def merge_files(self):
        self.merge_proteus_files()
        if self.Converter.trigger_info_file() != self.Converter.proteus_raw_file_path():
            self.merge_trigger_info_files()

    def fix_event_nrs(self):
        """fix the events numbers in the merged hdf5 file"""
        n_ev = [uproot.open(c.trigger_info_file())['Event'].num_entries for c in self.Converters]  # number of events for each run
        with h5py.File(self.Batch.FileName, 'r+') as f:
            ev = array(f['Tracks']['Events']).astype('i8')   # event number of the tracks which needs to be fixed
            i_tr = diff(concatenate([[0], where(diff(ev) < 0)[0] + 1, [ev.size]]))  # number of tracks for each run
            f['Tracks']['Events'][...] = ev + cumsum(append(0, n_ev[:-1])).repeat(i_tr).astype('u4')  # add the cumulated number of events from the previous runs

    def remove_aux_files(self):
        for conv in self.Converters:
            conv.remove_aux_files()
        self.Converter.remove_aux_files()

    def remove_raw_files(self):
        f = [i for c in self.Converters for i in c.raw_files if i.exists()]
        info(f'removing {len(f)} raw files ({byte2str(sum([i.stat().st_size for i in f]))})')
        for conv in self.Converters:
            conv.remove_raw_files(warn=False)


if __name__ == '__main__':

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
    cs = z.Converters if hasattr(z, 'Converters') else None

    if not args.t:
        cs = z.converters
        if len(cs):
            print_banner(f'Start converting runs {cs[0].Run} - {cs[-1].Run}', color=GREEN)
            z.run()
            print_banner('Finished Conversion!', color=GREEN)
        else:
            info('There is nothing to convert :-)\n', blank_lines=1)
