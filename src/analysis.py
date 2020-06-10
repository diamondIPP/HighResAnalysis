from draw import *
from glob import glob
from os.path import realpath, basename, expanduser
from shutil import copyfile
from sys import stdout
from time import time
from json import loads

# global test campaign
g_test_campaign = None


class Analysis(Draw):
    """ The analysis class provides default behaviour objects in the analysis framework and is the parent of all other analysis objects. """

    def __init__(self, testcampaign=None, verbose=False):
        Draw.__init__(self)

        self.Verbose = verbose
        self.Dir = get_program_dir()
        self.Config = self.load_config()
        self.Locations = loads(self.Config.get('MAIN', 'locations'))

        # test campaign
        self.TestCampaign = self.load_test_campaign(testcampaign)
        self.Location = self.get_test_campaigns()[self.TestCampaign]
        self.TCString = self.generate_tc_str()

        # Directories
        self.DataDir = join(self.get_raw_data_dir(), str(self.Location.lower()))
        self.TCDir = self.generate_tc_directory()
        self.PickleDir = join(self.Dir, self.Config.get('SAVE', 'pickle directory'))
        self.PickleSubDir = ''
        self.ResultsDir = join(self.Dir, 'results')

        # progress bar
        self.PBar = PBar()

    # ----------------------------------------
    # region INIT
    def load_config(self):
        parser = ConfigParser()
        config_file_path = join(self.Dir, 'config', 'main.ini')
        if not isfile(config_file_path):
            warning('The main config file "config/main.ini" does not exist! Using the default!')
            copyfile(join(self.Dir, 'config', 'main.default'), config_file_path)
        parser.read(config_file_path)
        return parser

    def get_raw_data_dir(self):
        return expanduser(self.Config.get('MAIN', 'data directory'))

    def generate_tc_str(self):
        return self.TestCampaign.strftime('%b %Y')

    def generate_tc_directory(self):
        return join(self.DataDir, self.TestCampaign.strftime('%Y-%m'))

    def load_test_campaign(self, testcampaign):
        global g_test_campaign
        if g_test_campaign is None:
            g_test_campaign = datetime.strptime(self.Config.get('MAIN', 'default test campaign') if testcampaign is None else testcampaign, '%Y%m')
        if g_test_campaign not in self.get_test_campaigns():
            critical('The Testcampaign {} does not exist!'.format(g_test_campaign))
        return g_test_campaign

    def get_test_campaigns(self):
        return {datetime.strptime(basename(path), '%Y-%m'): loc for loc in self.Locations for path in glob(join(self.get_raw_data_dir(), loc.lower(),  '*'))}

    def print_testcampaign(self):
        self.info('TESTCAMPAIGN: {}'.format(self.generate_tc_str()))
    # endregion INIT
    # ----------------------------------------

    def info(self, msg, next_line=True, prnt=True):
        if prnt and self.Verbose:
            t1 = time()
            t = datetime.now().strftime('%H:%M:%S')
            print('INFO: {t} --> {msg}'.format(t=t, msg=msg),)
            stdout.flush()
            if next_line:
                print()
            return t1

    def add_info(self, t, msg='Done'):
        if self.Verbose:
            print('{m} ({t:2.2f} s)'.format(m=msg, t=time() - t))

    def set_pickle_sub_dir(self, name):
        self.PickleSubDir = name

    def make_pickle_path(self, name='', suf='', sub_dir=None, run=None, dut=None, camp=None):
        directory = join(self.PickleDir, self.PickleSubDir if sub_dir is None else sub_dir)
        ensure_dir(directory)
        campaign = self.TCString if camp is None else camp
        run_str = str(run) if run is not None else self.RunPlan if hasattr(self, 'RunPlan') else ''
        run_str = run_str if run is not None or run_str else str(self.RunNumber) if hasattr(self, 'RunNumber') else ''
        # noinspection PyUnresolvedReferences
        dut = str(dut if dut is not None else self.DUT.Number if hasattr(self, 'DUT') and hasattr(self.DUT, 'Number') else '')
        return join(directory, '{}.pickle'.format('_'.join([v for v in [name, campaign, run_str, dut, str(suf)] if v])))

    def make_hdf5_path(self, *args, **kwargs):
        return self.make_pickle_path(*args, **kwargs).replace('pickle', 'hdf5')

    def print_start(self, run=None, prnt=True, tc=True):
        if prnt:
            ana_name = self.__class__.__name__.replace('Analysis', '')
            run = ' FOR RUN{} {}'.format('PLAN' if 'Coll' in ana_name else '', run) if run is not None else ''
            tc = ' OF {}'.format(self.TCString) if tc else ''
            print_banner('STARTING {} ANALYSIS{}{}'.format(ana_name.upper(), run, tc), symbol='~', color=GREEN)


def get_program_dir():
    return dirname(dirname(realpath(__file__)))


if __name__ == '__main__':
    z = Analysis()
