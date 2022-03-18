#!/usr/bin/env python
# --------------------------------------------------------
#       run module that handles the files and information of a single run
# created on October 5th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from src.analysis import Analysis
from src.converter import Converter
from src.run import Run
from utils import *


# todo revise
class CERNRun(Run):
    """ Run class containing all the information for a single run from the tree and the json file. """

    def __init__(self, run_number, dut_nr, tc_dir, config, single_mode=False):

        Run.__init__(self, run_number, dut_nr, tc_dir, config, single_mode)
        self.Tree = self.File.Get('Hits')
        
    def load_run_info(self):
        if self.SingleMode:
            info('Using single raw files')
            return None
        data = load_json(join(self.TCDir, self.Config.get('MAIN', 'run plan file')))
        data = data.values()[self.DUTNr][str(self.Number)]
        data['Batches'] = data['Batches'].split(',')
        return data

    def load_run_logs(self):
        data = load_json(join(self.TCDir, self.Config.get('MAIN', 'run log file')))
        data = {key: value for key, value in data.iteritems() if value['Run Number']['CMS'] not in ['-', '/', '\\', '']}
        if self.Info is None:
            return next(value for value in data.itervalues() if int(value['Run Number']['CMS']) == self.Number)
        return OrderedDict(sorted([(int(value['Run Number']['CMS']), value) for value in data.itervalues() if value['Batch'] in self.Info['Batches']]))

    def load_dut_name(self):
        return load_json(join(self.TCDir, self.Config.get('MAIN', 'run plan file')))['Plane{}'.format(self.DUTNr)]['Name']

    def load_raw_file_name(self):
        if self.SingleMode:
            return join(self.TCDir, 'cms-raw', 'ljutel_{}.root'.format(str(self.Number).zfill(3)))
        return join(self.TCDir, self.DUTName, 'run_{}.root'.format(str(self.Number).zfill(2)))
    
    def load_file_name(self):
        join(self.TCDir, self.DUTName, 'Clustered_{}.root'.format(str(self.Number).zfill(3 if self.SingleMode else 2)))

    def load_file(self):
        if not file_exists(self.FileName):
            self.convert_file()
        return TFile(self.FileName)

    def convert_file(self):
        if not file_exists(self.RawFileName) and not self.SingleMode:
            self.merge_root_files()
        warning('final root file "{}" does not exist. Starting Converter!'.format(self.FileName))
        converter = Converter(self.RawFileName, self.DUTNr, join(self.TCDir, self.DUTName), first_run=int(self.Logs.keys()[0]))
        converter.run()
        if not self.SingleMode:
            info('removing raw file "{}"'.format(self.RawFileName))
            remove(self.RawFileName)

    def merge_root_files(self):
        warning('raw file "{}" does not exist. Starting single file merger!'.format(self.RawFileName))
        single_files = [join(self.TCDir, 'cms-raw', 'ljutel_{}.root'.format(n)) for n in self.Logs.keys()]
        new_file = join(self.TCDir, self.DUTName, 'run_{}.root'.format(str(self.Number).zfill(2)))
        with open(devnull, 'w') as f:
            call([join(environ.get('ROOTSYS'), 'bin', 'hadd'), '-f', new_file] + single_files, stdout=f)
        info('successfully merged the single files to "{}"'.format(basename(new_file)))

    def print_run_info(self):
        for key, value in sorted(load_json(join(self.TCDir, self.Config.get('MAIN', 'run plan file'))).values()[self.DUTNr].iteritems(), key=lambda k: int(k[0])):
            print(key, value)


if __name__ == '__main__':
    a = Analysis()
    z = Run(301, 0, a.TCDir, a.Config, single_mode=True)
