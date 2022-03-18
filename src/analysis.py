from os.path import expanduser
from os import getcwd
from shutil import copyfile

from plotting.draw import *
from src.utils import PBar, Dir, print_banner, byte2str
from src.cut import Cuts


class BeamTest:

    def __init__(self, p: Path):

        self.Path = p
        self.Location = p.parts[-2].upper()
        self.T = datetime.strptime(p.stem, '%Y-%m')
        self.Year = self.T.year
        self.Tag = self.T.strftime('%Y%m')

    def __str__(self):
        return self.T.strftime('%b %Y')

    def __repr__(self):
        return f'Beam Test in {self} at {self.Location}'


def load_config():
    config_file_path = join(Dir, 'config', 'main.ini')
    if not isfile(config_file_path):
        warning('The main config file "config/main.ini" does not exist! Using the default!')
        copyfile(join(Dir, 'config', 'default.ini'), config_file_path)
    return Config(config_file_path)


class Analysis:
    """ The analysis class provides default behaviour objects in the analysis framework and is the parent of all other analysis objects. """

    Config = load_config()
    Locations = Config.get_list('MAIN', 'locations')

    DataDir = expanduser(Config.get('MAIN', 'data directory'))
    ResultsDir = join(Dir, 'results')
    MetaDir = join(Dir, Config.get('SAVE', 'meta directory'))

    def __init__(self, beamtest=None, meta_sub_dir='', verbose=False):

        self.Verbose = verbose

        self.BeamTest = self.load_test_campaign(beamtest)
        self.MetaSubDir = meta_sub_dir

        self.Cuts = Cuts()

        self.PBar = PBar()
        self.Draw = Draw(Analysis.Config.FilePath, self.Verbose)

    def __str__(self):
        return f'{self.__class__.__name__.replace("Analysis", "").upper()} ANALYSIS'

    def __repr__(self):
        return f'{self} of the {self.BeamTest!r}'

    # ----------------------------------------
    # region INIT
    @staticmethod
    def load_test_campaign(beamtest=None):
        bt = choose(beamtest, Analysis.Config.get('MAIN', 'default test campaign'))
        ps = [p for loc in Analysis.Locations for p in Path(Analysis.DataDir, loc.lower()).glob('*')]
        p = next((p for p in ps if bt == p.stem.replace('-', '')), None)
        return critical(f'The beamtest "{bt}" does not exist!') if p is None else BeamTest(p)

    @staticmethod
    def find_testcampaign():
        p = Path(getcwd())
        return BeamTest(p).Tag if p.parts[-2].upper() in Analysis.Locations else Analysis.Config.get('MAIN', 'default test campaign')

    def print_testcampaign(self):
        self.info(f'{self.BeamTest!r}')
    # endregion INIT
    # ----------------------------------------

    def info(self, msg, overlay=False, prnt=None):
        return info(msg, overlay, choose(prnt, self.Verbose))

    def add_info(self, t, msg='Done', prnt=None):
        add_to_info(t, msg, choose(prnt, self.Verbose))

    def set_meta_sub_dir(self, name):
        self.MetaSubDir = name

    def make_pickle_path(self, name='', suf='', sub_dir=None, run=None, dut=None, camp=None):
        directory = join(self.MetaDir, self.MetaSubDir if sub_dir is None else sub_dir)
        ensure_dir(directory)
        campaign = choose(camp, self.BeamTest.T.strftime('%Y%m'))
        dut = str(dut if dut is not None else self.DUT.Number if hasattr(self, 'DUT') and hasattr(self.DUT, 'Number') else '')
        run = choose(run, self.run_str)
        return join(directory, f'{"_".join([v for v in [name, campaign, run, dut, str(suf)] if v])}.pickle')

    @property
    def run_str(self):
        return self.RunPlan if hasattr(self, 'RunPlan') else str(self.Run) if hasattr(self, 'Run') else ''

    def get_meta_files(self):
        return [*Path(self.MetaDir).rglob(f'*_{self.BeamTest.Tag}_{self.run_str}*')] if self.run_str else []

    def meta_file_size(self):
        info(f'total size of metadata: {byte2str(sum(p.stat().st_size for p in self.get_meta_files()))}')

    def make_hdf5_path(self, *args, **kwargs):
        return self.make_pickle_path(*args, **kwargs).replace('pickle', 'hdf5')

    def print_start(self, run=None, prnt=True, tc=True):
        if prnt:
            ana_name = self.__class__.__name__.replace('Analysis', '')
            run = ' FOR RUN{} {}'.format('PLAN' if 'Coll' in ana_name else '', run) if run is not None else ''
            tc = f' OF {self.BeamTest}' if tc else ''
            print_banner('STARTING {} ANALYSIS{}{}'.format(ana_name.upper(), run, tc), symbol='~', color=GREEN)


if __name__ == '__main__':
    z = Analysis()
