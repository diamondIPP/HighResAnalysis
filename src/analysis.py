from os import getcwd
from shutil import copyfile

from plotting.save import *
from utility.utils import Dir, print_banner, byte2str


class BeamTest:

    def __init__(self, p: Path):

        self.Path = p
        self.Location = p.parts[-2].upper()
        self.T = datetime.strptime(p.stem, '%Y-%m')
        self.Year = self.T.year
        self.Tag = self.T.strftime('%Y%m')
        self.Name = self.T.strftime('%Y-%m')

    def __str__(self):
        return self.T.strftime('%b %Y')

    def __repr__(self):
        return f'Beam Test in {self} at {self.Location}'


def load_config():
    config_file_path = Dir.joinpath('config', 'main.ini')
    if not isfile(config_file_path):
        warning('The main config file "config/main.ini" does not exist! Using the default!')
        copyfile(Dir.joinpath('config', 'default.ini'), config_file_path)
    return Config(config_file_path)


class Analysis:
    """ The analysis class provides default behaviour objects in the analysis framework and is the parent of all other analysis objects. """

    Config = load_config()
    Locations = Config.get_value('data', 'beam tests', type)

    DataDir = Path(Config.get('data', 'dir')).expanduser()
    ResultsDir = Dir.joinpath('results')
    MetaDir = Dir.joinpath(Config.get('SAVE', 'meta directory'))

    def __init__(self, beamtest=None, meta_sub_dir='', verbose=False):

        self.Verbose = verbose

        self.BeamTest = self.load_test_campaign(beamtest)
        self.MetaSubDir = meta_sub_dir

        self.Draw = SaveDraw(self, results_dir=self.BeamTest.Tag)

    def __str__(self):
        return f'{self.__class__.__name__.replace("Analysis", "").upper()} ANALYSIS'.strip(' ')

    def __repr__(self):
        return f'{self} ({self.BeamTest})'

    # ----------------------------------------
    # region INIT
    @property
    def server_save_dir(self):
        return

    @staticmethod
    def init_locations():
        for loc, beam_tests in Analysis.Locations.items():
            p = Analysis.DataDir.joinpath(loc.lower())
            p.mkdir(exist_ok=True)
            for bt in beam_tests:
                bt = str(bt)
                p.joinpath(f'{bt[:4]}-{bt[4:]}').mkdir(exist_ok=True)

    @staticmethod
    def load_test_campaign(beamtest=None):
        bt = choose(beamtest, Analysis.Config.get('data', 'default test campaign')).replace('-', '')
        Analysis.init_locations()
        ps = [p for loc in Analysis.Locations for p in Path(Analysis.DataDir, loc.lower()).glob('*')]
        p = next((p for p in ps if bt == p.stem.replace('-', '')), None)
        return BeamTest(p)

    @staticmethod
    def find_testcampaign():
        p = Path(getcwd())
        return BeamTest(p).Tag if p.parts[-2].upper() in Analysis.Locations else Analysis.Config.get('data', 'default test campaign')

    def print_testcampaign(self):
        self.info(f'{self.BeamTest!r}')
    # endregion INIT
    # ----------------------------------------

    def info(self, msg, blank_lines=0, endl=True, prnt=None):
        return info(msg, blank_lines, endl, choose(prnt, self.Verbose))

    def add_info(self, t, msg='Done', prnt=None):
        add_to_info(t, msg, prnt=choose(prnt, self.Verbose))

    def make_pickle_path(self, name='', suf='', sub_dir=None, run=None, dut=None, camp=None):
        directory = self.MetaDir.joinpath(self.MetaSubDir if sub_dir is None else sub_dir)
        ensure_dir(directory)
        campaign = choose(camp, self.BeamTest.T.strftime('%Y%m'))
        dut = str(dut if dut is not None else self.DUT.Number if hasattr(self, 'DUT') and hasattr(self.DUT, 'Number') else '')
        run = choose(run, self.run_str)
        return directory.joinpath(f'{"_".join([v for v in [name, campaign, run, dut, str(suf)] if v])}.pickle')

    @property
    def unit_str(self):
        return f'run {self.Run}' if hasattr(self, 'Run') else ''

    @property
    def run_str(self):
        return str(self.Run) if hasattr(self, 'Run') else ''

    def get_meta_files(self):
        return [*Path(self.MetaDir).rglob(f'*_{self.BeamTest.Tag}_{self.run_str}*')] if self.run_str else []

    @property
    def meta_file_size(self):
        return sum(p.stat().st_size for p in self.get_meta_files())

    def print_meta_file_size(self):
        info(f'total size of metadata: {byte2str(self.meta_file_size)}')

    def print_meta_file_sizes(self):
        for p in self.get_meta_files():
            info(f'{p}: {byte2str(p.stat().st_size)}')
        self.print_meta_file_size()

    def remove_metadata(self):
        s = self.meta_file_size
        for p in self.get_meta_files():
            remove_file(p)
        info(f'removed {byte2str(s)} of meta files')

    def make_hdf5_path(self, *args, **kwargs):
        return self.make_pickle_path(*args, **kwargs).with_suffix('.hdf5')

    def print_start(self):
        print_banner(f'STARTING {self!r}', symbol='~', color=GREEN)

    def create_run_config(self):
        if self.BeamTest.Location == 'CERN':
            from src.spreadsheet import make_cern_run_log
            make_cern_run_log(self.BeamTest.Path.stem)
        elif self.BeamTest.Location == 'DESY':
            from src.spreadsheet import make_desy_run_log
            make_desy_run_log()


if __name__ == '__main__':
    z = Analysis()
