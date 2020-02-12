# --------------------------------------------------------
#       UTILITY FUNCTIONS
# created on June 19th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from __future__ import print_function

import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True  # disable ROOT overwriting the help settings...

from os.path import isfile, exists, isdir
from os import makedirs, _exit
from ConfigParser import ConfigParser
from datetime import datetime
from ROOT import TFile, gROOT
from json import load
from collections import OrderedDict
from uncertainties import ufloat
from uncertainties.core import Variable, AffineScalarFunc
from numpy import average, sqrt, array, arange, mean
from progressbar import Bar, ETA, FileTransferSpeed, Percentage, ProgressBar


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
        print('{ov}{head} {t} --> {msg}'.format(t=get_t_str(), msg=msg, head='{}INFO:{}'.format(GREEN, ENDC), ov='\033[1A\r' if overlay else ''))


def warning(msg):
    print('{head} {t} --> {msg}'.format(t=get_t_str(), msg=msg, head='{}WARNING:{}'.format(YELLOW, ENDC)))


def critical(msg):
    print('{head} {t} --> {msg}\n'.format(t=get_t_str(), msg=msg, head='{}CRITICAL:{}'.format(RED, ENDC)))
    _exit(1)


def file_exists(filename):
    return isfile(filename)


def dir_exists(path):
    return isdir(path)


def time_stamp(dt, off=None):
    t = float(dt.strftime('%s'))
    return t if off is None else t - (off if off > 1 else dt.utcoffset().seconds)


def average_list(lst, n):
    return [mean(lst[i:i+n]) for i in arange(0, len(lst), n)] if n > 1 else lst


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


def colored(string, color):
    return '{}{}{}'.format(color, string, ENDC)


def print_banner(msg, symbol='~', new_lines=1, color=None):
    msg = '{} |'.format(msg)
    print(colored('{n}{delim}\n{msg}\n{delim}{n}'.format(delim=len(str(msg)) * symbol, msg=msg, n='\n' * new_lines), color))


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


def get_object(name):
    return gROOT.FindObject(name)


def remove_letters(string):
    return filter(lambda x: x.isdigit(), string)


def remove_digits(string):
    return filter(lambda x: not x.isdigit(), string)


def get_last_canvas():
    try:
        return gROOT.GetListOfCanvases()[-1]
    except IndexError:
        warning('There is no canvas is in the list...')


def set_z_range(zmin, zmax):
    c = get_last_canvas()
    h = c.GetListOfPrimitives()[1]
    h.GetZaxis().SetRangeUser(zmin, zmax)


def set_axes_range(xmin, xmax, ymin, ymax):
    set_x_range(xmin, xmax)
    set_y_range(ymin, ymax)


def set_x_range(xmin, xmax):
    c = get_last_canvas()
    h = c.GetListOfPrimitives()[1]
    h.GetXaxis().SetRangeUser(xmin, xmax)


def set_y_range(ymin, ymax):
    c = get_last_canvas()
    h = c.GetListOfPrimitives()[1]
    h.GetYaxis().SetRangeUser(ymin, ymax)


def normalise_histo(histo, x_range=None, from_min=False):
    h = histo
    x_axis = h.GetXaxis()
    x_axis.SetRangeUser(*x_range) if x_range is not None else do_nothing()
    min_bin = h.GetMinimumBin() if from_min else 0
    integral = h.Integral(min_bin, h.GetNbinsX() - 1)
    return scale_histo(h, integral)


def scale_histo(histo, value=None, to_max=False, x_range=None):
    h = histo
    maximum = h.GetBinContent(h.GetMaximumBin())
    if x_range is not None:
        h.GetXaxis().SetRangeUser(*x_range) if x_range is not None else do_nothing()
        maximum = h.GetBinContent(h.GetMaximumBin())
        h.GetXaxis().UnZoom()
    value = maximum if to_max else value
    if value:
        h.Scale(1. / value)
    return h


def mean_sigma(values, weights=None):
    """ Return the weighted average and standard deviation. values, weights -- Numpy ndarrays with the same shape. """
    if len(values) == 1:
        value = make_ufloat(values[0])
        return value.n, value.s
    weights = [1] * len(values) if weights is None else weights
    if type(values[0]) in [Variable, AffineScalarFunc]:
        weights = [1 / v.s for v in values]
        values = array([v.n for v in values], 'd')
    if all(weight == 0 for weight in weights):
        return [0, 0]
    avrg = average(values, weights=weights)
    variance = average((values - avrg) ** 2, weights=weights)  # Fast and numerically precise
    return avrg, sqrt(variance)


def make_ufloat(tup):
    if type(tup) in [Variable, AffineScalarFunc]:
        return tup
    return ufloat(tup[0], tup[1]) if type(tup) in [tuple, list] else ufloat(tup, 0)


def get_root_vec(tree, n, ind=0, dtype=None):
    vec = tree.GetVal(ind)
    vec.SetSize(n)
    return array(vec, dtype=dtype)


def get_root_vecs(tree, n, n_ind, dtype=None):
    return [get_root_vec(tree, n, i, dtype) for i in xrange(n_ind)]


class PBar:
    def __init__(self):
        self.PBar = None
        self.Widgets = ['Progress: ', Percentage(), ' ', Bar(marker='>'), ' ', ETA(), ' ', FileTransferSpeed()]

    def start(self, n):
        self.PBar = ProgressBar(widgets=self.Widgets, maxval=n).start()

    def update(self, i):
        self.PBar.update(i + 1)
        if i == self.PBar.maxval - 1:
            self.finish()

    def finish(self):
        self.PBar.finish()
