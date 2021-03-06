from draw import *
from glob import glob
from os.path import realpath, basename, expanduser
from shutil import copyfile
from json import loads
from cut import Cuts

# global test campaign
g_test_campaign = None


class Analysis(Draw):
    """ The analysis class provides default behaviour objects in the analysis framework and is the parent of all other analysis objects. """

    def __init__(self, testcampaign=None, verbose=False):

        self.Verbose = verbose
        self.Dir = get_program_dir()
        self.Config = self.load_config()
        self.Locations = loads(self.Config.get('MAIN', 'locations'))
        Draw.__init__(self, verbose, self.Config)

        # test campaign
        self.TestCampaign = self.load_test_campaign(testcampaign)
        self.Location = self.get_test_campaigns()[self.TestCampaign]
        self.TCString = self.generate_tc_str()

        # Directories
        self.DataDir = join(self.get_raw_data_dir(), str(self.Location.lower()))
        self.TCDir = self.generate_tc_directory()
        self.MetaDir = join(self.Dir, self.Config.get('SAVE', 'meta directory'))
        self.MetaSubDir = ''
        self.ResultsDir = join(self.Dir, 'results')

        self.Cuts = Cuts()

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

    def info(self, msg, overlay=False, prnt=None):
        return info(msg, overlay, choose(prnt, self.Verbose))

    def add_info(self, t, msg='Done', prnt=None):
        add_to_info(t, msg, choose(prnt, self.Verbose))

    def set_pickle_sub_dir(self, name):
        self.MetaSubDir = name

    def make_pickle_path(self, name='', suf='', sub_dir=None, run=None, dut=None, camp=None):
        directory = join(self.MetaDir, self.MetaSubDir if sub_dir is None else sub_dir)
        ensure_dir(directory)
        campaign = choose(camp, self.TestCampaign.strftime('%Y%m'))
        # noinspection PyUnresolvedReferences
        dut = str(dut if dut is not None else self.DUT.Number if hasattr(self, 'DUT') and hasattr(self.DUT, 'Number') else '')
        return join(directory, '{}.pickle'.format('_'.join([v for v in [name, campaign, self.make_run_str(run), dut, str(suf)] if v])))

    def make_run_str(self, run_str=None):
        run = self.RunPlan if hasattr(self, 'RunPlan') else self.Run.Number if hasattr(self, 'Run') and hasattr(self.Run, 'Number') else ''
        return choose(run_str, str(run))

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
