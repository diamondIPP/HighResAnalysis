#!/usr/bin/env python
# --------------------------------------------------------
#       handles the cuts for the high rate analysis
# created on July 10th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from numpy import array, all, in1d
from utils import print_table, warning


class Cuts:
    """ Class that holds several cuts with functionality to combine them """

    def __init__(self):

        self.Cuts = {}

    def __call__(self, cut=None):
        return self.generate() if cut is None else None if not cut else cut

    def __add__(self, other=None):
        return self.generate() if other is None else all([self.generate(), other], axis=0)

    def generate(self):
        cuts = [cut.Values for cut in self.Cuts.values() if cut.Level < 80]
        return all(cuts, axis=0).flatten()

    def register(self, name, values, level, description=None):
        self.Cuts[name] = Cut(name, values, level, description)

    def get(self, name):
        return self.Cuts[name]

    def show(self, raw=False):
        rows = [[cut.Name, '{:5d}'.format(cut.Level), cut.Value if raw else cut.Description] for cut in self.Cuts.values()]
        print_table([row for row in rows if row[2]], ['Cut Name', 'Level', 'Description'])

    @staticmethod
    def make_mask(x, y, masked_pixels):
        data = x.astype('i') * 10000 + y  # make unique number out of the tuple... Is there a way to compare tuples?
        mx, my = masked_pixels.T.astype('i')
        mask = mx * 10000 + my
        return in1d(data, mask, invert=True)


class Cut:
    """ Base class to describe a single cut """

    def __init__(self, name, values, level, description=None):

        self.Name = name
        self.Values = array(values)
        self.Level = level
        self.Description = description
        self.Size = self.Values.size

    def __call__(self):
        return self.Values

    def __add__(self, other=None):
        s = len(other)
        if s != self.Size:
            warning('cut array has incorrect size ({}), {} required'.format(s, self.Size))
            return self.Values
        return self.Values if other is None else all([self.Values, other], axis=0)

    def __str__(self):
        return '{}, {} cut {}/{} entries: {}'.format(self.Level, self.Name, self.Values.nonzero()[0].size, self.Size, self.Description)

    def __repr__(self):
        return self.__str__()

    def get_data(self, data):
        return data if self.Values is None else data[self.Values]
