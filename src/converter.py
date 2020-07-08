#!/usr/bin/env python
# --------------------------------------------------------
#       adds clustering and charge to trees created with pXar
# created on August 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from utils import *
from os.path import expanduser
from argparse import ArgumentParser
from numpy import average, sum, append


class Converter:

    def __init__(self, data_dir, run_number, config):

        self.RunNumber = run_number
        self.DataDir = data_dir
        self.Config = config

        self.RawFilePath = self.load_raw_file_name()
        self.SoftDir = expanduser(self.Config.get('SOFTWARE', 'dir'))

        self.PBar = PBar()

    def load_raw_file_name(self):
        pass

    # TODO think about general methods for calibration...


class Hit:

    def __init__(self, x, y):

        self.X = x
        self.Y = y

        self.Charge = None

    def set_charge(self, value):
        self.Charge = value

    def __str__(self):
        return 'Hit: {0} {1}, Charge: {2:1.2f}vcal'.format(self.X, self.Y, self.Charge)


class Cluster:

    def __init__(self, seed_hit):
        self.Hits = array([seed_hit])

    def add_hit(self, hit):
        self.Hits = append(self.Hits, [hit], axis=0)

    def hit_is_adjacent(self, hit):
        """ returns: if any of the existing cluster hits has a distance of 1 to the given hit"""
        return any(sqrt(sum(abs(hit[:2] - self.Hits[:, :2]), axis=1)) <= 1.1)

    def get_charge(self):
        return sum(self.Hits[:, 2])

    def get_size(self):
        return self.Hits.shape[0]

    def get_x(self):
        return average(self.Hits[:, 0], weights=self.Hits[:, 2])

    def get_y(self):
        return average(self.Hits[:, 1], weights=self.Hits[:, 2])

    def __repr__(self):
        return 'Cluster with {} hits at [{:.1f}, {:.1f}] with charge {:.1f} vcals'.format(self.get_size(), self.get_x(), self.get_y(), self.get_charge())


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('filename', nargs='?', default='')
    # noinspection PyTypeChecker
    p.add_argument('plane', nargs='?', default=0, type=int)
    args = p.parse_args()
    z = Converter('/scratch2/cern/2018-10/cms-raw/ljutel_110.root', 0, '/scratch2/cern/2018-10/II6-B6')
    # z.run()
