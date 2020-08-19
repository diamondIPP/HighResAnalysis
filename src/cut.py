#!/usr/bin/env python
# --------------------------------------------------------
#       handles the cuts for the high rate analysis
# created on July 10th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from numpy import array, all, in1d, invert
from utils import print_table, warning, get_base_dir, load_config, join, critical, make_list
from json import loads


class Cuts:
    """ Class that holds several cuts with functionality to combine them """

    def __init__(self):

        self.Config = load_config(join(get_base_dir(), 'config', 'cut'))
        self.Cuts = {}

    def __call__(self, cut=None):
        cut = cut.Values if isinstance(cut, Cut) else cut
        values = array(cut)
        return values if values.size > 1 else self.generate() if cut is None else ...

    def __add__(self, other=None):
        return self.generate() if other is None else all([self.generate(), other], axis=0)

    def get_config(self, option, lst=False, dtype=None):
        if option not in self.Config.options('CUT'):
            critical('option "{}" not found in cut config, please set it!')
        value = self.Config.get('CUT', option)
        return loads(value) if lst else value if dtype is None else dtype(value)

    def generate(self):
        cuts = [cut.Values for cut in self.Cuts.values() if cut.Level < 80]
        return all(cuts, axis=0).flatten() if len(cuts) else ...

    def register(self, name, values=None, level=None, description=None):
        if isinstance(name, Cut):
            cut = name
            cut.set_level(level)
            self.Cuts[cut.Name] = cut
        else:
            self.Cuts[name] = Cut(name, values, level, description)

    def get(self, name):
        return self.Cuts[name]

    def get_special(self, exclude):
        exclude = make_list(exclude)
        cuts = [cut.Values for cut in self.Cuts.values() if cut.Level < 80 and cut.Name not in exclude]
        return all(cuts, axis=0).flatten() if len(cuts) else ...

    def show(self, raw=False):
        rows = [[cut.Name, '{:5d}'.format(cut.Level), cut.get_p_str(), cut.Value if raw else cut.Description] for cut in self.Cuts.values()]
        print_table([row for row in rows if row[2]], ['Cut Name', 'Level', 'P', 'Description'])

    @staticmethod
    def make_mask(x, y, masked_pixels):
        data = x.astype('i') * 10000 + y  # make unique number out of the tuple... Is there a way to compare tuples?
        mx, my = masked_pixels.T.astype('i')
        mask = mx * 10000 + my
        return in1d(data, mask, invert=True)


class Cut:
    """ Base class to describe a single cut """

    def __init__(self, name, values, level=99, description=None):

        self.Name = name
        self.Values = array(values)
        self.Level = level
        self.Description = description
        self.Size = self.Values.size
        self.P = invert(self.Values).nonzero()[0].size / self.Size

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

    def __str__(self):
        return '{}, {} cut, {}: {}'.format(self.Level, self.Name, self.get_p_str(), self.Description)

    def __repr__(self):
        return self.__str__()

    def get_p_str(self):
        return '{:.1f}%'.format(self.P * 100)

    def set_level(self, level):
        if level is not None:
            self.Level = level

    @staticmethod
    def make(cut):
        if isinstance(cut, Cut):
            return cut.Values
        return cut if cut else ...
