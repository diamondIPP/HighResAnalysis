# --------------------------------------------------------
#       UTILITY FUNCTIONS
# created on June 19th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from os.path import isfile, exists
from os import makedirs, _exit
from ConfigParser import ConfigParser
from datetime import datetime
from ROOT import TFile
from json import load
from collections import OrderedDict


type_dict = {'int32': 'I',
             'uint16': 's',
             'float64': 'D',
             'int64': 'L'}


GREEN = '\033[92m'
ENDC = '\033[0m'
YELLOW = '\033[93m'
RED = '\033[91m'


def get_t_str():
    return datetime.now().strftime('%H:%M:%S')


def info(msg, overlay=False, prnt=True):
    if prnt:
        print '{ov}{head} {t} --> {msg}'.format(t=get_t_str(), msg=msg, head='{}INFO:{}'.format(GREEN, ENDC), ov='\033[1A\r' if overlay else '')


def warning(msg):
    print '{head} {t} --> {msg}'.format(t=get_t_str(), msg=msg, head='{}WARNING:{}'.format(YELLOW, ENDC))


def critical(msg):
    print '{head} {t} --> {msg}\n'.format(t=get_t_str(), msg=msg, head='{}CRITICAL:{}'.format(RED, ENDC))
    _exit(1)


def file_exists(filename):
    return isfile(filename)


def round_down_to(num, val):
    return int(num) / val * val


def ensure_dir(path):
    if not exists(path):
        info('Creating directory: {d}'.format(d=path))
        makedirs(path)


def is_num(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


def print_banner(msg, symbol='=', new_lines=True):
    print '{n}{delim}\n{msg}\n{delim}{n}'.format(delim=(len(str(msg)) + 10) * symbol, msg=msg, n='\n' if new_lines else '')


def do_nothing():
    pass


def load_config(name, ext='ini'):
    parser = ConfigParser()
    parser.read('{}.{}'.format(name, ext))
    return parser


def has_root():
    try:
        import ROOT
        return True
    except ImportError:
        return False


def read_root_file(filename):
    if file_exists(filename):
        return TFile(filename)
    critical('The file: "{}" does not exist...'.format(filename))


def load_json(filename, ordered=None):
    if not file_exists(filename):
        return {}
    with open(filename) as f:
        return load(f, object_hook=None if ordered is None else OrderedDict)


def do(fs, pars, exe=-1):
    fs, pars = ([fs], [pars]) if type(fs) is not list else (fs, pars)
    exe = pars if exe == -1 else [exe]
    for f, p, e in zip(fs, pars, exe):
        f(p) if e is not None else do_nothing()
