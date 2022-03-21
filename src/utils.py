# --------------------------------------------------------
#       UTILITY FUNCTIONS
# created on June 19th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True  # disable ROOT overwriting the help settings...

from os.path import isfile, exists, isdir, dirname, realpath, join, basename
from os import makedirs, environ, remove, devnull
from subprocess import call
from ROOT import TFile
from json import load, loads
from collections import OrderedDict
from uncertainties import ufloat
from uncertainties.core import Variable, AffineScalarFunc
from numpy import average, sqrt, array, arange, mean, exp, concatenate, count_nonzero, zeros, sin, cos, dot, log2, log10, array_split
from progressbar import Bar, ETA, FileTransferSpeed, Percentage, ProgressBar, SimpleProgress, Widget
import h5py
import pickle
from copy import deepcopy
from inspect import signature
from functools import wraps
from plotting.utils import info, warning, critical, add_to_info
from datetime import timedelta, datetime
from multiprocessing import Pool, cpu_count


Dir = dirname(dirname(realpath(__file__)))

type_dict = {'int32': 'I',
             'uint16': 's',
             'float64': 'D',
             'int64': 'L'}


GREEN = '\033[92m'
WHITE = '\033[98m'
ENDC = '\033[0m'
YELLOW = '\033[93m'
CYAN = '\033[96m'
RED = '\033[91m'
UP1 = '\033[1A'
ERASE = '\033[K'


def move_up(n):
    print('\033[{}A'.format(n))


def file_exists(filename):
    return isfile(filename)


def dir_exists(path):
    return isdir(path)


def time_stamp(dt, off=None):
    t = float(dt.strftime('%s'))
    return t if off is None else t - (off if off > 1 else dt.utcoffset().seconds)


def average_list(lst, n):
    return [mean(lst[i:i+n]) for i in arange(0, len(lst), n)] if n > 1 else lst


def round_down_to(num, val=1):
    return int(num) // val * val


def round_up_to(num, val=1):
    return int(num) // val * val + val


def get_base_dir():
    return dirname(dirname(realpath(__file__)))


def ensure_dir(path):
    if not exists(path):
        info('Creating directory: {d}'.format(d=path))
        makedirs(path)


def remove_file(file_path, string=None):
    if file_exists(file_path):
        warning('removing {}'.format(choose(string, file_path)))
        remove(file_path)


def isint(x):
    try:
        return float(x) == int(x)
    except (ValueError, TypeError):
        return False


def is_iter(v):
    try:
        iter(v)
        return True
    except TypeError:
        return False


def is_num(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


def colored(string, color):
    return '{}{}{}'.format(color, string, ENDC)


def print_banner(msg, symbol='~', new_lines=1, color=WHITE):
    msg = '{} |'.format(msg)
    print(colored('{n}{delim}\n{msg}\n{delim}{n}'.format(delim=len(str(msg)) * symbol, msg=msg, n='\n' * new_lines), color))


def prime_factors(n):
    factors = []
    while n % 2 == 0:  # even dividers
        factors.append(2)
        n /= 2
    for i in range(3, int(sqrt(n)) + 1, 2):  # odd dividers
        while n % i == 0:
            factors.append(i)
            n /= i
    if n > 2:
        factors.append(int(n))
    return factors


def do_nothing():
    pass


def has_root():
    try:
        import ROOT
        return True
    except ImportError:
        return False


def get_rot_matrix(alpha):
    return array([[cos(alpha), -sin(alpha)], [sin(alpha), cos(alpha)]])


def rotate(alpha, v):
    return dot(get_rot_matrix(alpha), v)


def read_root_file(filename):
    if file_exists(filename):
        return TFile(filename)
    critical('The file: "{}" does not exist...'.format(filename))


def load_json(filename, ordered=None):
    if not file_exists(filename):
        warning('json file does not exist: {}'.format(filename))
        return {}
    with open(filename) as f:
        return load(f, object_hook=None if ordered is None else OrderedDict)


def do(fs, pars, exe=-1):
    fs, pars = ([fs], [pars]) if type(fs) is not list else (fs, pars)
    exe = pars if exe == -1 else [exe]
    for f, p, e in zip(fs, pars, exe):
        f(p) if e is not None else do_nothing()


def choose(v, default, decider='None', *args, **kwargs):
    use_default = decider is None if decider != 'None' else v is None
    if callable(default) and use_default:
        default = default(*args, **kwargs)
    return default if use_default else v


def remove_letters(string):
    return ''.join(filter(lambda x: x.isdigit(), string))


def remove_digits(string):
    return ''.join(filter(lambda x: not x.isdigit(), string))


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
    values = values.astype('d')
    avrg = average(values, weights=weights)
    variance = average((values - avrg) ** 2, weights=weights)  # Fast and numerically precise
    return avrg, sqrt(variance)


def interpolate_two_points(x1, y1, x2, y2, name=''):
    # f = p1*x + p0
    p1 = (y1 - y2) / (x1 - x2)
    p0 = y1 - x1 * p1
    w = abs(x2 - x1)
    fit_range = array(sorted([x1, x2])) + [-w / 3., w / 3.]
    f = TF1('fpol1{}'.format(name), 'pol1', *fit_range)
    f.SetParameters(p0, p1)
    return f


def interpolate_x(x1, x2, y1, y2, y):
    p1 = get_p1(x1, x2, y1, y2)
    p0 = get_p0(x1, y1, p1)
    return (y - p0) / p1 if p1 else 0


def interpolate(x1, x2, y1, y2, x):
    x1, x2, y1, y2 = [float(i) for i in [x1, x2, y1, y2]]
    p1 = get_p1(float(x1), x2, y1, y2)
    p0 = get_p0(x1, y1, p1)
    return p1 * x + p0


def get_p1(x1, x2, y1, y2):
    return (y1 - y2) / (x1 - x2) if x1 != x2 else 0


def get_p0(x1, y1, p1):
    return y1 - x1 * p1


def calc_eff(k=0, n=0, values=None):
    values = array(values) if values is not None else None
    if n == 0 and not values.size:
        return zeros(3)
    k = float(k if values is None else count_nonzero(values))
    n = float(n if values is None else values.size)
    m = (k + 1) / (n + 2)
    mode = k / n
    s = sqrt(((k + 1) / (n + 2) * (k + 2) / (n + 3) - ((k + 1) ** 2) / ((n + 2) ** 2)))
    return array([mode, max(s + (mode - m), 0), max(s - (mode - m), 0)]) * 100


def make_ufloat(tup):
    if type(tup) in [Variable, AffineScalarFunc]:
        return tup
    return ufloat(tup[0], tup[1]) if type(tup) in [tuple, list] else ufloat(tup, 0)


def byte2str(v):
    n = int(log2(v) // 10) if v else 0
    return '{:1.1f} {}'.format(v / 2 ** (10 * n), ['B', 'kB', 'MB', 'GB'][n])


def ev2str(v):
    n = int(log10(v) // 3)
    return f'{v / 10 ** (3 * n):.{2 if n > 1 else 0}f}{["", "k", "M"][n]}'


def get_root_vec(tree, n=0, ind=0, dtype=None, var=None, cut=''):
    if var is not None:
        n = tree.Draw(var, cut, 'goff')
    vec = tree.GetVal(ind)
    vec.SetSize(n)
    return array(list(vec), dtype=dtype)


def get_root_vecs(tree, n, n_ind, dtype=None):
    dtypes = [None] * n_ind if dtype is None else dtype
    return [get_root_vec(tree, n, i, dtypes[i]) for i in range(n_ind)]


def make_list(value, dtype=None):
    v = array([choose(value, [])]).flatten()
    return v.tolist() if dtype == list else v.astype(dtype) if dtype is not None else v


def gauss(x, scale, mean_, sigma, off=0):
    return scale * exp(-.5 * ((x - mean_) / sigma) ** 2) + off


def do_hdf5(path, func, redo=False, *args, **kwargs):
    if file_exists(path) and redo:
        remove_file(path)
    if file_exists(path) and not redo:
        return h5py.File(path, 'r')['data']
    else:
        data = func(*args, **kwargs)
        f = h5py.File(path, 'w')
        info('creating {}'.format(join(basename(dirname(path)), basename(path))))
        f.create_dataset('data', data=data)
        return f['data']


def do_pickle(path, func, value=None, redo=False, *args, **kwargs):
    if value is not None:
        with open(path, 'wb') as f:
            pickle.dump(value, f)
        return value
    try:
        if file_exists(path) and not redo:
            with open(path, 'rb') as f:
                return pickle.load(f)
    except ImportError:
        pass
    ret_val = func(*args, **kwargs)
    with open(path, 'wb') as f:
        pickle.dump(ret_val, f)
    return ret_val


def print_table(rows, header=None, footer=None, prnt=True):
    head, foot = [choose([v], zeros((0, len(rows[0]))), v) for v in [header, footer]]
    t = concatenate([head, rows, foot]).astype('str')
    # t = array(rows, dtype=str) if header is None else concatenate((array([header], dtype=str), array(rows, dtype=str)))
    col_width = [len(max(t[:, i], key=len)) for i in range(t.shape[1])]
    total_width = sum(col_width) + len(col_width) * 3 + 1
    hline = '{}'.format('~' * total_width)
    if prnt:
        for i, row in enumerate(t):
            if i in [0] + choose([1], [], header) + choose([t.shape[0] - 1], [], footer):
                print(hline)
            print('| {r} |'.format(r=' | '.join(word.ljust(n) for word, n in zip(row, col_width))))
        print('{}\n'.format(hline))
    return rows


def merge_root_files(files, new_file_name):
    with open(devnull, 'w') as f:
        call([join(environ.get('ROOTSYS'), 'bin', 'hadd'), '-f', new_file_name] + files, stdout=f)
    info('successfully merged the single files to "{}"'.format(basename(new_file_name)))


def update_pbar(func):
    @wraps(func)
    def my_func(*args, **kwargs):
        value = func(*args, **kwargs)
        if args[0].PBar is not None and args[0].PBar.PBar is not None and not args[0].PBar.is_finished():
            args[0].PBar.update()
        return value
    return my_func


class PBar(object):
    def __init__(self, start=None, counter=False, t=None):
        self.PBar = None
        self.Widgets = self.init_widgets(counter, t)
        self.Step = 0
        self.N = 0
        self.start(start)

    def __reduce__(self):
        return self.__class__, (None, False, None), (self.Widgets, self.Step, self.N)

    def __setstate__(self, state):
        self.Widgets, self.Step, self.N = state
        if self.N:
            self.PBar = ProgressBar(widgets=self.Widgets, maxval=self.N).start()
            self.update(self.Step) if self.Step > 0 else do_nothing()

    @staticmethod
    def init_widgets(counter, t):
        return ['Progress: ', SimpleProgress('/') if counter else Percentage(), ' ', Bar(marker='>'), ' ', ETA(), ' ', FileTransferSpeed() if t is None else EventSpeed(t)]

    def start(self, n, counter=None, t=None):
        if n is not None:
            self.Step = 0
            self.PBar = ProgressBar(widgets=self.Widgets if t is None and counter is None else self.init_widgets(counter, t), maxval=n).start()
            self.N = n

    def update(self, i=None):
        i = self.Step if i is None else i
        if i >= self.PBar.maxval:
            return
        self.PBar.update(i + 1)
        self.Step += 1
        if i == self.PBar.maxval - 1:
            self.finish()

    def set_last(self):
        if self.PBar:
            self.PBar.currval = self.N
            self.PBar.finished = True

    def finish(self):
        self.PBar.finish()

    def is_finished(self):
        return self.PBar.currval == self.N

    def eta(self, i, h, m, s=0):
        self.PBar.start_time = time_stamp(datetime.now() - timedelta(hours=h, minutes=m, seconds=s))
        self.update(i - 1)


class EventSpeed(Widget):
    """Widget for showing the event speed (useful for slow updates)."""

    def __init__(self, t='s'):
        self.unit = t
        self.factor = {'s': 1, 'min': 60, 'h': 60 * 60}[t]

    def update(self, pbar):
        value = 0
        if pbar.seconds_elapsed > 2e-6 and pbar.currval > 2e-6:
            value = pbar.currval / pbar.seconds_elapsed * self.factor
        return f'{value:4.1f} E/{self.unit}'


def prep_kw(dic, **default):
    d = deepcopy(dic)
    for kw, value in default.items():
        if kw not in d:
            d[kw] = value
    return d


def make_suffix(ana, *values):
    suf_vals = [ana.get_short_name(suf) if type(suf) is str and suf.startswith('TimeIntegralValues') else suf for suf in values]
    return '_'.join(str(int(val) if isint(val) else val.GetName() if hasattr(val, 'GetName') else val) for val in suf_vals if val is not None)


def prep_suffix(f, args, kwargs, suf_args, field=None):
    def_pars = signature(f).parameters
    names, values = list(def_pars.keys()), [par.default for par in def_pars.values()]
    i_arg = (arange(len([n for n in names if n not in ['self', '_redo']])) if suf_args == 'all' else make_list(loads(str(suf_args)))) + 1
    suf_vals = [args[i] if len(args) > i else kwargs[names[i]] if names[i] in kwargs else values[i] for i in i_arg]
    suf_vals += [getattr(args[0], str(field))] if field is not None and hasattr(args[0], field) else []
    return make_suffix(args[0], *suf_vals)


def load_pickle(file_name):
    with open(file_name, 'rb') as f:
        return pickle.load(f)


def save_pickle(*pargs, print_dur=False, low_rate=False, high_rate=False, suf_args='[]', field=None, verbose=False, **pkwargs):
    def inner(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            if '_no_save' in kwargs:
                return func(*args, **kwargs)
            run = args[0].Run.get_high_rate_run(high=not low_rate) if low_rate or high_rate else None
            pickle_path = args[0].make_pickle_path(*pargs, **prep_kw(pkwargs, run=run, suf=prep_suffix(func, args, kwargs, suf_args, field)))
            info(f'Pickle path: {pickle_path}', prnt=verbose)
            redo = (kwargs['_redo'] if '_redo' in kwargs else False) or (kwargs['show'] if 'show' in kwargs else False)
            if file_exists(pickle_path) and not redo:
                return load_pickle(pickle_path)
            prnt = print_dur and (kwargs['prnt'] if 'prnt' in kwargs else True)
            t = (args[0].info if hasattr(args[0], 'info') else info)(f'{args[0].__class__.__name__}: {func.__name__.replace("_", " ")} ...', endl=False, prnt=prnt)
            value = func(*args, **kwargs)
            with open(pickle_path, 'wb') as f:
                pickle.dump(value, f)
            (args[0].add_to_info if hasattr(args[0], 'add_to_info') else add_to_info)(t, prnt=prnt)
            return value
        return wrapper
    return inner


def save_hdf5(*pargs, arr=False, dtype=None, suf_args='[]', **pkwargs):
    def inner(f):
        @wraps(f)
        def wrapper(*args, **kwargs):
            file_path = args[0].make_hdf5_path(*pargs, **prep_kw(pkwargs, suf=prep_suffix(f, args, kwargs, suf_args)))
            redo = kwargs['_redo'] if '_redo' in kwargs else False
            if file_exists(file_path) and not redo:
                d = h5py.File(file_path, 'r')['data']
                return array(d) if arr else d
            remove_file(file_path)
            data = f(*args, **kwargs)
            hf = h5py.File(file_path, 'w')
            hf.create_dataset('data', data=data.astype(choose(dtype, data.dtype)))
            return array(hf['data']) if arr else hf['data']
        return wrapper
    return inner


def parallel(fp, what=''):
    def inner(f):
        @wraps(f)
        def my_f(ana, *args, **kwargs):
            with Pool() as pool:
                ana.info(f'generate {what}')
                r = f(ana, *args, **kwargs)
                d, fargs = (r[0], r[1:]) if len(r) > 1 else (r, [])
                pbar = PBar(d.shape[0] // cpu_count())
                f_ = getattr(ana.__class__, fp)
                result = pool.starmap(_parallel, [(f_, d, i, pbar, *fargs) for i in array_split(arange(d.shape[0]), cpu_count())])
                return concatenate(result)
        return my_f
    return inner


def _parallel(f, d, i, pbar, *args):
    ret = []
    use_pbar = i[0] == 0
    for v in d[i]:
        if use_pbar:
            pbar.update()
        ret.append(f(v, *args))
    return ret
