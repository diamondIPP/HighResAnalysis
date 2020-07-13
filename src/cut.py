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

    def generate(self):
        cuts = [cut.Values for cut in self.Cuts.values()]
        return all(cuts, axis=0).flatten()

    def register(self, name, values, level, description=None):
        self.Cuts[name] = Cut(name, values, level, description)

    def get(self, name):
        return self.Cuts[name].Values

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

    def __call__(self, cut_array):
        return cut_array if self.Values is None else cut_array[self.Values]

    def get_data(self, data):
        return data if self.Values is None else data[self.Values]
