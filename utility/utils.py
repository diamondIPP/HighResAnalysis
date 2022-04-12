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
from numpy import sqrt, array, arange, mean, exp, concatenate, zeros, log2, log10, array_split, ndarray, full, frombuffer
from progressbar import Bar, ETA, FileTransferSpeed, Percentage, ProgressBar, SimpleProgress, Widget
import h5py
import pickle
from copy import deepcopy
from inspect import signature
from functools import wraps
from plotting.utils import info, warning, critical, add_to_info, get_kw
from datetime import timedelta, datetime
from multiprocessing import Pool, cpu_count
from hashlib import md5
from pathlib import Path


Dir = Path(__file__).resolve().parent.parent

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
    return path


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


def choose(v, default, decider='None', *args, **kwargs):
    use_default = decider is None if decider != 'None' else v is None
    if callable(default) and use_default:
        default = default(*args, **kwargs)
    return default if use_default else v


def remove_letters(string):
    return ''.join(filter(lambda x: x.isdigit(), string))


def remove_digits(string):
    return ''.join(filter(lambda x: not x.isdigit(), string))


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


def make_ufloat(n, s=0):
    return array([ufloat(*v) for v in array([n, s]).T]) if is_iter(n) else n if is_ufloat(n) else ufloat(n, s)


def is_ufloat(value):
    return type(value) in [Variable, AffineScalarFunc]


def byte2str(v):
    n = int(log2(v) // 10) if v else 0
    return '{:1.1f} {}'.format(v / 2 ** (10 * n), ['B', 'kB', 'MB', 'GB'][n])


def ev2str(v):
    n = int(log10(v) // 3)
    return f'{v / 10 ** (3 * n):.{2 if n > 1 else 0}f}{["", "k", "M"][n]}'


def get_buf(buf, n, dtype=None):
    return frombuffer(buf, dtype=buf.typecode, count=n).astype(dtype)


def get_tree_vec(tree, var, cut='', dtype=None, nentries=None, firstentry=0):
    strings = make_list(var)
    n = tree.Draw(':'.join(strings), cut, 'goff', choose(nentries, tree.kMaxEntries), firstentry)
    dtypes = dtype if type(dtype) in [list, ndarray] else full(len(strings), dtype)
    vals = [get_buf(tree.GetVal(i), n, dtypes[i]) for i in range(len(strings))]
    return vals[0] if len(vals) == 1 else vals


def make_list(value, dtype=None):
    v = value if is_iter(value) and not type(value) is str else array([choose(value, [])]).flatten()
    return v.tolist() if dtype == list else v.astype(dtype) if dtype is not None else v


def uarr2n(arr):
    return array([i.n for i in arr]) if len(arr) and is_ufloat(arr[0]) else arr


def uarr2s(arr):
    return array([i.s for i in arr]) if len(arr) and is_ufloat(arr[0]) else arr


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
        if args[0].PBar is not None and args[0].PBar.PBar is not None and not args[0].PBar.is_finished() and not get_kw('_no_update', kwargs):
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


def get_field(obj, field: str):
    if '.' in field:
        return get_field(getattr(obj, field.split('.')[0]), '.'.join(field.split('.')[1:]))
    return getattr(obj, field) if hasattr(obj, field) else None


def make_suffix(*values):
    vals = [md5(val).hexdigest() if type(val) is ndarray else f'{val:.0f}' if isint(val) else val for val in values if val is not None]
    return '_'.join(str(val) for val in vals)


def prep_suffix(f, ana, args, kwargs, suf_args, field=None):
    def_pars = signature(f).parameters
    names, values = list(def_pars.keys())[1:], [par.default for par in def_pars.values()][1:]  # first par is class instance
    i_arg = arange(len([n for n in names if n not in ['self', '_redo']])) if suf_args == 'all' else make_list(loads(str(suf_args)))
    suf_vals = [args[i] if len(args) > i else kwargs[names[i]] if names[i] in kwargs else values[i] for i in i_arg]
    suf_vals += [] if field is None else [get_field(ana, field)]
    return make_suffix(*suf_vals)


def load_pickle(file_name):
    with open(file_name, 'rb') as f:
        return pickle.load(f)


def save_pickle(*pargs, print_dur=False, low_rate=False, high_rate=False, suf_args='[]', field=None, verbose=False, **pkwargs):
    def inner(func):
        @wraps(func)
        def wrapper(ana, *args, **kwargs):
            if not get_kw('_save', kwargs, default=True):
                return func(ana, *args, **kwargs)
            run = ana.Run.get_high_rate_run(high=not low_rate) if low_rate or high_rate else None
            pickle_path = ana.make_pickle_path(*pargs, **prep_kw(pkwargs, run=run, suf=prep_suffix(func, ana, args, kwargs, suf_args, field)))
            info(f'Pickle path: {pickle_path}', prnt=verbose)
            redo = (kwargs['_redo'] if '_redo' in kwargs else False) or (kwargs['show'] if 'show' in kwargs else False)
            if file_exists(pickle_path) and not redo:
                return load_pickle(pickle_path)
            prnt = print_dur and (kwargs['prnt'] if 'prnt' in kwargs else True)
            t = (ana.info if hasattr(ana, 'info') else info)(f'{ana.__class__.__name__}: {func.__name__.replace("_", " ")} ...', endl=False, prnt=prnt)
            value = func(ana, *args, **kwargs)
            with open(pickle_path, 'wb') as f:
                pickle.dump(value, f)
            (ana.add_to_info if hasattr(ana, 'add_to_info') else add_to_info)(t, prnt=prnt)
            return value
        return wrapper
    return inner


def save_hdf5(*pargs, arr=False, dtype=None, suf_args='[]', field=None, verbose=False, **pkwargs):
    def inner(f):
        @wraps(f)
        def wrapper(ana, *args, **kwargs):
            file_path = ana.make_hdf5_path(*pargs, **prep_kw(pkwargs, suf=prep_suffix(f, ana, args, kwargs, suf_args, field)))
            info(f'HDF5 path: {file_path}', prnt=verbose)
            redo = kwargs['_redo'] if '_redo' in kwargs else False
            if file_exists(file_path) and not redo:
                d = h5py.File(file_path, 'r')['data']
                return array(d) if arr else d
            remove_file(file_path)
            data = f(ana, *args, **kwargs)
            hf = h5py.File(file_path, 'w')
            hf.create_dataset('data', data=data.astype(choose(dtype, data.dtype)))
            return array(hf['data']) if arr else hf['data']
        return wrapper
    return inner


def parallel(fp, what='something'):
    def inner(f):
        @wraps(f)
        def my_f(ana, *args, **kwargs):
            with Pool() as pool:
                ana.info(f'generate {what} for {ana}')
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
