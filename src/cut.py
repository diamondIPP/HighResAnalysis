#!/usr/bin/env python
# --------------------------------------------------------
#       handles the cuts for the high rate analysis
# created on July 10th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from typing import Any
from numpy import array, all, invert, ones, log10, count_nonzero, cumsum, ceil, ndarray

from utility.utils import print_table, make_list, choose, Dir, is_iter, critical
from plotting.utils import Config, warning


class Cuts:
    """ Class that holds several cuts with functionality to combine them """

    def __init__(self):

        self.Dir = Dir.joinpath('cuts')
        self.Config = self.init_config()
        self.Cuts = {}
        self.make()

    def __call__(self, cut=None, **k):
        cut = cut.Values if isinstance(cut, Cut) else cut
        values = array(cut)
        return values if values.size > 1 else self.generate() if cut is None else ...

    def __add__(self, other=None):
        return self.generate() if other is None else all([self.generate(), other], axis=0)

    def __getitem__(self, v):
        return self.get(v)()

    def __iter__(self):
        return iter(self.Cuts.values())

    def __str__(self):
        return self.__class__.__name__

    def __repr__(self):
        k, n = (count_nonzero(self()), self().size) if self.n else (0, 0)
        return f'{self} instance with {self.n} cuts' + (f', selecting {k}/{n} events ({100 * k / n:.1f} %)' if n else '')

    @property
    def config_file(self):
        return self.Dir.joinpath('cut.ini')

    def init_config(self):
        if not self.config_file.exists():
            critical(f'analysis config file "{self.config_file}" does not exist!')
        return Config(self.config_file)

    def make(self, redo=False):
        pass

    def get_config(self, option, dtype=str, default=None):
        v = self.Config.get_value(option, dtype=dtype, default=default)
        return array(v) if is_iter(v) else v

    def generate(self):
        cuts = [cut.Values for cut in self.Cuts.values() if cut.Level < 80]
        return all(cuts, axis=0).flatten() if len(cuts) else ...

    def register(self, name, values=None, level=None, description=None):
        if isinstance(name, Cut):
            cut = name
            cut.set_level(choose(level, cut.Level))
            self.Cuts[cut.Name] = cut
        else:
            if not all(values):
                self.Cuts[name] = Cut(name, values, level, description)

    def get(self, name):
        return self.Cuts[name]

    def set(self, name, values):
        if name not in self.Cuts:
            return warning(f'{name} does not exist in {self} {list(self.Cuts.keys())}')
        cut = self.Cuts[name]
        self.Cuts[name] = Cut(name, values, cut.Level, cut.Description)

    @property
    def n(self):
        return len([cut for cut in self.Cuts.values() if cut.Level < Cut.MaxLevel])

    @property
    def size(self):
        return self().size

    @property
    def none(self):
        return ones(self.size, '?')

    @property
    def names(self):
        return [c.Name for c in self]

    def get_consecutive(self):
        return cumsum(sorted(filter(lambda x: x.Level < Cut.MaxLevel, self.Cuts.values())))

    def get_nofid(self, cut=None, fid=False):
        return self(cut) if fid else self.exclude('fid', cut)

    def include(self, *names):
        return all([self[n] for n in names], axis=0)

    def exclude(self, exclude, cut=None):
        if cut is not None:
            return self(cut)
        exclude = make_list(exclude)
        cuts = [cut.Values for cut in self.Cuts.values() if cut.Level < 80 and cut.Name not in exclude]
        return all(cuts, axis=0).flatten() if len(cuts) else ...

    def add(self, cut: ndarray):
        return self() & cut

    def remove(self, name):
        return self.Cuts.pop(name) if name in self.Cuts else None

    def show(self, raw=False):
        rows = [[cut.Name, '{:5d}'.format(cut.Level), cut.Size, cut.nstr, cut.exstr, cut.pstr, cut.Value if raw else cut.Description] for cut in self.get_consecutive()]
        c = Cut('all', self.generate(), 0, 'final cut')
        print_table([row for row in rows if row[2]], ['Cut Name', 'Level', 'Size', 'N', 'Exl', 'P', 'Description'], [c.Name, '', c.Size, c.nstr, c.exstr, c.pstr, c.Description])


class Cut:
    """ Base class to describe a single cut """

    MaxLevel = 80

    def __init__(self, name='', values: Any = 1, level=99, description=None, n0=None):

        self.Name = name
        self.Values = ones(values, '?') if type(values) is int else values.Values if isinstance(values, Cut) else array(values)
        self.Level = level
        self.Description = description
        self.Size = self.Values.size
        self.P = self.calc_p()
        self.N = round(self.Size * (1 - self.P))
        self.N0 = choose(n0, self.Size)  # n excluded before adding

    def __call__(self):
        return self.Values

    def __add__(self, other=None):
        if type(other) is bool or other is ...:
            return self
        values = other.Values if isinstance(other, Cut) else array(other)
        if values.size != self.Size:
            warning(f'could not add cuts! Array has incorrect size ({values.size}), {self.Size} required')
            return self
        n, d, lev = (other.Name, other.Description, other.Level) if isinstance(other, Cut) else (self.Name, self.Description, self.Level)
        return Cut(n, self.Values & values, lev, d, self.N)

    def __gt__(self, other):
        return self.Level > other.Level

    def __str__(self):
        return self.Name

    def __repr__(self):
        return f'{self.Level}, {self.Name} cut, {self.pstr}: {self.Description}'

    def calc_p(self):
        return invert(self.Values).nonzero()[0].size / self.Size

    @property
    def pstr(self):
        return f'{self.P * 100:.1f}%'

    @property
    def nstr(self):
        return f'{self.N:{int(ceil(log10(self.Size)))}d}'

    @property
    def exstr(self):
        return f'{self.N0 - self.N:{int(ceil(log10(self.Size)))}d}'

    def set_p(self, p):
        self.P = p

    def set_level(self, level):
        if level is not None:
            self.Level = level

    def set_values(self, values):
        self.Values = values
        self.P = self.calc_p()

    @staticmethod
    def make(cut):
        if isinstance(cut, Cut):
            return cut.Values
        return cut if cut else ...

    @staticmethod
    def add(*cuts):
        c = Cut('added', cuts[0])
        for cut in cuts[1:]:
            c += cut
        return c()
