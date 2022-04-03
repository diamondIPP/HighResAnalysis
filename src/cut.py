#!/usr/bin/env python
# --------------------------------------------------------
#       handles the cuts for the high rate analysis
# created on July 10th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from typing import Any
from numpy import array, all, in1d, invert, ones, count_nonzero, max

from src.utils import print_table, warning, join, make_list, choose, Dir, is_iter
from plotting.utils import Config
from plotting.draw import make_box_args, Draw, prep_kw, TCutG


class Cuts:
    """ Class that holds several cuts with functionality to combine them """

    def __init__(self):

        self.Config = self.init_config()
        self.Cuts = {}
        self.make()

    def __call__(self, cut=None, pl=None):
        cut = cut.Values if isinstance(cut, Cut) else cut
        values = array(cut)
        return values if values.size > 1 else self.generate() if cut is None else ...

    def __add__(self, other=None):
        return self.generate() if other is None else all([self.generate(), other], axis=0)

    def __getitem__(self, v):
        return self.get(v)()

    def __repr__(self):
        k, n = (count_nonzero(self()), self().size) if self.n else (0, 0)
        return f'{self.__class__.__name__} instance with {self.n} cuts' + (f', selecting {k}/{n} events ({100 * k / n:.1f} %)' if n else '')

    def init_config(self):
        return Config(join(Dir, 'cuts', 'cut.ini'))

    def make(self, redo=False):
        pass

    def get_config(self, option, dtype=str):
        v = self.Config.get_value(option, dtype=dtype)
        return array(v) if is_iter(v) else v

    def get_fid_config(self, surface=False):
        p = self.get_config('surface fiducial' if surface else 'fiducial')
        p = make_box_args(*p[[0, 2, 1, 3]]) if p.size == 4 else p  # unpack short box notation
        p[p == max(p, axis=1).reshape((-1, 1))] += 1  # extend one pixel to the top and right
        return p - .5  # pixel centre is at the integer

    def get_fid(self, surface=False, **dkw):
        x, y = self.get_fid_config(surface)
        return Draw.polygon(x, y, **prep_kw(dkw, show=False, line_color=2, width=2, name=f'fid{surface:d}'))

    def draw_fid(self, surface=False, **dkw):
        self.get_fid(surface, **prep_kw(dkw, show=True))

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

    @property
    def n(self):
        return len([cut for cut in self.Cuts.values() if cut.Level < Cut.MaxLevel])

    def get_consecutive(self):
        cuts = sorted([cut for cut in self.Cuts.values() if cut.Level < Cut.MaxLevel])
        consecutive = ones(cuts[0].Size, dtype=bool)
        p0 = 0
        for cut in cuts:
            consecutive = Cut.add(cut, consecutive)
            cut.set_values(consecutive)
            cut.set_p(cut.P - p0)
            p0 += cut.P
        return cuts

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

    def show(self, raw=False):
        rows = [[cut.Name, '{:5d}'.format(cut.Level), cut.Size, cut.get_p_str(), cut.Value if raw else cut.Description] for cut in self.get_consecutive()]
        c = Cut('all', self.generate(), 0, 'final cut')
        print_table([row for row in rows if row[2]], ['Cut Name', 'Level', 'Size', 'P', 'Description'], [c.Name, '', c.Size, c.get_p_str(), c.Description])

    def make_pixel_mask(self, x, y):
        data = x.astype('i') * 10000 + y  # make unique number out of the tuple... Is there a way to compare tuples?
        mx, my = array(self.get_config('mask')).T
        mask = mx * 10000 + my
        return in1d(data, mask, invert=True)

    @staticmethod
    def point_in_polygon(p, poly: TCutG):
        return poly.IsInside(*p)


class Cut:
    """ Base class to describe a single cut """

    MaxLevel = 80

    def __init__(self, name, values: Any, level=99, description=None):

        self.Name = name
        self.Values = values.Values if isinstance(values, Cut) else array(values)
        self.Level = level
        self.Description = description
        self.Size = self.Values.size
        self.P = self.calc_p()
        self.N = round(self.Size * (1 - self.P))

    def __call__(self):
        return self.Values

    def __add__(self, other=None):
        values = array(other)
        if isinstance(other, Cut):
            values = other.Values
        elif array(other).size == 1:
            return self
        if values.size != self.Size:
            warning('cut array has incorrect size ({}), {} required'.format(values.size, self.Size))
            return self
        return Cut('add', all([self.Values, values], axis=0))

    def __gt__(self, other):
        return self.Level > other.Level

    def __str__(self):
        return '{}, {} cut, {}: {}'.format(self.Level, self.Name, self.get_p_str(), self.Description)

    def __repr__(self):
        return self.__str__()

    def calc_p(self):
        return invert(self.Values).nonzero()[0].size / self.Size

    def get_p_str(self):
        return '{:.1f}%'.format(self.P * 100)

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
