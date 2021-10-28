#!/usr/bin/env python
# --------------------------------------------------------
#       handles the cuts for the high rate analysis
# created on July 10th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from numpy import array, all, in1d, invert, ones
from utils import print_table, warning, get_base_dir, load_config, join, critical, make_list, choose, datetime
from json import loads
from typing import Any


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

    def set_config(self, test_campaign: datetime, dut_name):
        config = load_config(join(get_base_dir(), 'config', 'cut{}'.format(test_campaign.strftime('%Y%m'))))
        if dut_name not in config.sections():
            critical('detector "{}" not found in cut config, please add it!'.format(dut_name))
        self.Config = config._sections[dut_name]

    def get_config(self, option, lst=False, dtype=None):
        if option not in self.Config:
            critical('option "{}" not found in cut config, please set it!'.format(option))
        value = self.Config[option]
        return array(loads(value)) if lst else value if dtype is None else dtype(value)

    def  get_fid_config(self, surface=False):
        return self.get_config('{}fiducial'.format('surface ' if surface else ''), lst=True)

    def generate(self):
        cuts = [cut.Values for cut in self.Cuts.values() if cut.Level < 80]
        return all(cuts, axis=0).flatten() if len(cuts) else ...

    def register(self, name, values=None, level=None, description=None):
        if isinstance(name, Cut):
            cut = name
            cut.set_level(choose(level, cut.Level))
            self.Cuts[cut.Name] = cut
        else:
            self.Cuts[name] = Cut(name, values, level, description)

    def get(self, name):
        return self.Cuts[name]

    def get_consecutive(self):
        cuts = sorted([cut for cut in self.Cuts.values() if cut.Level < 80])
        consecutive = ones(cuts[0].Size, dtype=bool)
        p0 = 0
        for cut in cuts:
            consecutive = Cut.add(cut, consecutive)
            cut.set_values(consecutive)
            cut.set_p(cut.P - p0)
            p0 += cut.P
        return cuts

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

    @staticmethod
    def make_mask(x, y, masked_pixels):
        data = x.astype('i') * 10000 + y  # make unique number out of the tuple... Is there a way to compare tuples?
        mx, my = masked_pixels.T.astype('i')
        mask = mx * 10000 + my
        return in1d(data, mask, invert=True)


class Cut:
    """ Base class to describe a single cut """

    def __init__(self, name, values: Any, level=99, description=None):

        self.Name = name
        self.Values = values.Values if isinstance(values, Cut) else array(values)
        self.Level = level
        self.Description = description
        self.Size = self.Values.size
        self.P = self.calc_p()

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
