#!/usr/bin/env python
# --------------------------------------------------------
#       alignment and tracking of telescope data with proteus
# created on April 14th 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from pathlib import Path
from os import chdir
import toml
from numpy import array
from plotting.utils import info, warning, choose, remove_file
from utility.utils import print_banner, GREEN, print_elapsed_time, wraps
from subprocess import check_call, CalledProcessError
from shutil import copytree


def init_toml(name):
    def inner(func):
        @wraps(func)
        def wrapper(self, arg=None):
            default = self.ConfigDir.joinpath(f'{name}.toml')
            if arg is None:
                return default
            d = toml.load(default)
            func(self, arg, d)
            tmp = self.ConfigDir.joinpath(f'tmp-{name[:3]}.toml')
            with open(tmp, 'w') as f:
                toml.dump(d, f)
            return tmp
        return wrapper
    return inner


class Proteus:
    """ Alignment and tracking of telescope data.
        STEP 1: noisescan
        STEP 2: alignment
        STEP 3: tracking
        STEP 4: track-matching  """

    Si_Detectors = ['D8']

    def __init__(self, soft_dir, data_dir, cfg_dir, raw_file, max_events=None, skip_events=None, dut_pos=None, duts=None):

        # DIRECTORIES
        self.SoftDir = Path(soft_dir)
        self.DataDir = Path(data_dir)
        self.ConfigDir = Path(cfg_dir)

        # CONFIG
        self.NTelPlanes = sum([d['name'].startswith('M') for d in toml.load(self.ConfigDir.joinpath('device.toml'))['sensors']])
        self.MaxDUTs = len(toml.load(self.ConfigDir.joinpath('geometry.toml'))['sensors']) - self.NTelPlanes  # default geo has all sensors
        self.Geo = self.init_geo(dut_pos)
        self.Device = self.init_device(duts)
        self.Ana = self.init_ana(duts)

        self.RawFilePath = Path(raw_file)
        self.RunNumber = int(''.join(filter(lambda x: x.isdigit(), self.RawFilePath.stem)))
        self.Out = self.DataDir.joinpath(f'tracked-{self.RunNumber:04d}')           # name for proteus
        self.OutFilePath = self.Out.with_name(f'{self.Out.name}-trees.root')        # final file
        self.HistFilePath = self.Out.with_name(f'{self.Out.name}-hists.root')       # file with histograms
        self.TrackName = self.DataDir.joinpath(f'clustered-{self.RunNumber:04d}')   # tracking file for alignment

        self.N = max_events
        self.S = skip_events

        self.AlignSteps = self.align_steps()
        self.Steps = [(self.noise_scan, self.toml_name('tel', 'mask', 'mask')), (self.align, self.toml_name()), (self.recon, self.OutFilePath)]

    def __repr__(self):
        return f'Proteus interface for run {self.RunNumber} ({self.RawFilePath.name})'

    def __del__(self):
        remove_file(self.ConfigDir.joinpath('tmp-geo.toml'), prnt=False)
        remove_file(self.ConfigDir.joinpath('tmp-dev.toml'), prnt=False)

    # ----------------------------------------
    # region INIT
    def __create_default_cfg(self):
        copytree(self.ConfigDir.with_name('default'), self.ConfigDir)

    @init_toml('geometry')
    def init_geo(self, dut_pos, data=None):
        for i in range(self.MaxDUTs):  # remove non-existing DUTs
            if i not in dut_pos:
                data['sensors'].pop(i - self.MaxDUTs)
        for i, dic in enumerate(data['sensors']):  # fix ids
            dic['id'] = i

    @init_toml('device')
    def init_device(self, duts, data=None):
        s = [dic for dic in data['sensors'] if dic['name'].startswith('M')]  # only select TEL planes
        s += [{'type': 'CMSPixel-Si' if dut in Proteus.Si_Detectors else 'CMSPixel-Dia', 'name': f'C{i}'} for i, dut in enumerate(duts)]  # add all given DUTs
        data['sensors'] = s

    @init_toml('analysis')
    def init_ana(self, duts, data=None):
        data['recon']['extrapolation_ids'] = list(range(self.NTelPlanes + len(duts)))

    # endregion INIT
    # ----------------------------------------

    def align_steps(self):
        if not self.ConfigDir.is_dir():
            warning(f'no proteus config found for {self.ConfigDir.stem}, creating default!')
            self.__create_default_cfg()
        return list(toml.load(str(self.ConfigDir.joinpath('align.toml')))['align'])

    def align_files(self):
        return [self.ConfigDir.joinpath(name) for name in [self.toml_name('all', 'mask', 'mask'), self.toml_name()]]

    def alignment(self, step=-1):
        return toml.load(str(self.ConfigDir.joinpath(self.toml_name(self.AlignSteps[step]))))

    def z_positions(self, raw=False):
        d = toml.load(str(self.ConfigDir.joinpath('geometry.toml' if raw else self.toml_name())))
        return array([s['offset'][-1] for s in d['sensors']])

    # ----------------------------------------
    # region MISC
    def toml_name(self, name=None, d='alignment', typ='geo'):
        return self.ConfigDir.joinpath(d, f'{self.AlignSteps[-1] if name is None else name}-{typ}.toml')

    def make_empty_masks(self, cfg):
        for section in cfg:
            m = self.toml_name(section, d='mask', typ='mask')
            if not m.exists():
                m.write_text('[[sensors]]\nid = 0\nmasked_pixels = []\n')

    def remove_root_files(self):
        for f in self.ConfigDir.rglob('*.root'):
            remove_file(f)

    def remove_alignment(self):
        for f in self.ConfigDir.joinpath('alignment').glob('*.toml'):
            remove_file(f)

    def remove_mask(self):
        for f in self.ConfigDir.joinpath('mask').glob('*.toml'):
            remove_file(f)
    # endregion MISC
    # ----------------------------------------

    # ----------------------------------------
    # region RUN
    def run(self, prog, out: Path, cfg=None, geo=None, dev=None, section=None, f=None, n=None, s=None):
        old_dir = Path.cwd()
        chdir(self.ConfigDir)  # proteus needs to be in the directory where all the toml files are (for the default args)...
        cfg = '' if cfg is None else f' -c {str(cfg).replace(".toml", "")}.toml'
        section = '' if section is None else f' -u {section}'
        geo = f' -g {choose(geo, self.Geo)}'
        dev = f' -d {choose(dev, self.Device)}'
        n = f' -n {choose(n, self.N)}' if choose(n, self.N) is not None else ''
        s = f' -s {choose(s, self.S)}' if choose(s, self.S) is not None else ''
        cmd = f'{self.SoftDir.joinpath("bin", prog)} {choose(f, self.RawFilePath)} {out}{cfg}{dev}{geo}{section}{n}{s}'
        info(cmd)
        try:
            check_call(cmd, shell=True)
        except CalledProcessError:
            warning(f'{prog} failed!')
        chdir(old_dir)

    def noise_scan(self):
        """ step 1: find noisy pixels. """
        d = Path('mask')
        cfg = toml.load(str(self.ConfigDir.joinpath('noisescan.toml')))['noisescan']
        self.ConfigDir.joinpath(d).mkdir(exist_ok=True)
        self.make_empty_masks(cfg)
        for section in cfg:
            print_banner(f'Starting noise scan for {section}', color=GREEN)
            self.run('pt-noisescan', out=d.joinpath(section), cfg='noisescan', section=section)

    def align(self, step=None, force=False):
        """ step 2: align the telescope in several steps. """
        t = info('Starting alignment ...')
        d = Path('alignment')
        self.ConfigDir.joinpath(d).mkdir(exist_ok=True)
        for i in range(len(self.AlignSteps)) if step is None else [step]:
            step = self.AlignSteps[i]
            if not self.toml_name(step).exists() or force:
                self.run('pt-align', out=d.joinpath(step), geo=self.toml_name(self.AlignSteps[i - 1]) if i else None, section=step, cfg='align')
            else:
                warning('geo file already exists!')
        print_elapsed_time(t)

    def recon(self, cfg=None, geo=-1):
        """ step 3: based on the alignment generate the tracks with proteus. """
        self.Out.parent.mkdir(exist_ok=True)
        self.run('pt-recon', out=self.Out, cfg=choose(cfg, self.Ana), geo=None if geo is None else self.toml_name(self.AlignSteps[geo]))

    def track(self):
        """ tracking and clustering for the event alignment. """
        self.Out.parent.mkdir(exist_ok=True)
        self.run('pt-track', out=self.TrackName)
    # endregion RUN
    # ----------------------------------------


if __name__ == '__main__':
    from src.analysis import Analysis, Dir

    a = Analysis()
    sdir = Path(a.Config.get('SOFTWARE', 'dir')).expanduser().joinpath(a.Config.get('SOFTWARE', 'proteus'))
    f_ = a.BeamTest.Path.joinpath('data', f'run{11:06d}.root')
    z = Proteus(sdir, a.BeamTest.Path.joinpath('proteus'), Dir.joinpath('proteus'), f_, a.Config.getint('align', 'max events'), a.Config.getint('align', 'skip events'))
