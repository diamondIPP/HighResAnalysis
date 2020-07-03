#!/usr/bin/env python
# --------------------------------------------------------
#       adds clustering and charge to trees created with pXar
# created on August 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from utils import *
from os.path import expanduser
from argparse import ArgumentParser
from numpy import average


class Converter:

    def __init__(self, data_dir, run_number, config, calibration):

        self.RunNumber = run_number
        self.DataDir = data_dir
        self.Config = config
        self.Calibration = calibration

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

    def __init__(self):

        self.Hits = []

    def __str__(self):
        return 'Cluster of size {}, Charge: {}'.format(self.size(), self.charge())

    def add_hit(self, hit):
        self.Hits.append(hit)

    def size(self):
        return len(self.Hits)

    def charge(self):
        return sum(hit.Charge for hit in self.Hits)

    def seed_hit(self):
        return max(self.Hits, key=lambda hit: hit.Charge)

    def x(self):
        return average([hit.X for hit in self.Hits], weights=[hit.Charge for hit in self.Hits])

    def y(self):
        return average([hit.Y for hit in self.Hits], weights=[hit.Charge for hit in self.Hits])


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('filename', nargs='?', default='')
    # noinspection PyTypeChecker
    p.add_argument('plane', nargs='?', default=0, type=int)
    args = p.parse_args()
    z = Converter('/scratch2/cern/2018-10/cms-raw/ljutel_110.root', 0, '/scratch2/cern/2018-10/II6-B6', None)
    # z.run()
