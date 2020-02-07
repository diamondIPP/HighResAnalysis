#!/usr/bin/env python
# --------------------------------------------------------
#       adds clustering and charge to trees created with pXar
# created on August 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import TFile, vector, TF1
from utils import *
from os.path import basename, join, dirname, expanduser
from argparse import ArgumentParser
from collections import OrderedDict
from numpy import array, average, zeros
from progressbar import Bar, ETA, FileTransferSpeed, Percentage, ProgressBar
from pickle import load, dump
from draw import Draw, ufloat
from glob import glob


class Converter:

    def __init__(self, filepath, config):

        self.RawFilePath = filepath
        self.TCDir = dirname(dirname(filepath))
        self.Config = config
        self.SoftDir = expanduser(self.Config.get('SOFTWARE', 'dir'))
        self.PBar = PBar()

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
    z = Converter('/scratch2/cern/2018-10/cms-raw/ljutel_110.root', 0, '/scratch2/cern/2018-10/II6-B6')
    # z.run()
