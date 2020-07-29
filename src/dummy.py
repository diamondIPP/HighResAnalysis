#!/usr/bin/env python
# --------------------------------------------------------
#       module to create dummy file for analysis
# created on July 29th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from utils import *
from numpy.random import poisson, normal
from numpy import round
from dut import Plane
from converter import Converter


class Dummy:

    def __init__(self, tc_dir, n_tel, n_dut, config):

        self.FileName = join(tc_dir, 'data', 'dummy.hdf5')
        self.File = None
        self.NTelPlanes = n_tel
        self.NDUTPlanes = n_dut
        self.NEvents = 0
        self.Config = config

    def init_file(self):
        self.File = h5py.File(self.FileName, 'w')

    def load_file(self):
        return h5py.File(self.FileName, 'r')

    def set_n_events(self, n):
        self.NEvents = int(n)

    def make(self, n_events=1e6, beam_x=1, beam_y=1):
        self.set_n_events(n_events)
        remove_file(self.FileName)
        self.init_file()
        self.add_tracks()
        self.add_time()
        self.add_telescope(beam_x, beam_y)

    def add_tracks(self):
        g = self.File.create_group('Tracks')
        g.create_dataset('NTracks', data=poisson(lam=1.5, size=self.NEvents).astype('u1'))

    def add_time(self):
        g = self.File.create_group('Event')
        g.create_dataset('Time', data=arange(0, 1e5, .00165, dtype='f2')[:self.NEvents])

    def add_telescope(self, beam_x, beam_y):
        p = Plane(0, self.Config)
        for i in range(self.NTelPlanes):
            g = self.File.create_group('Plane{}'.format(i))
            g_hits = g.create_group('Hits')
            n_hits = poisson(lam=6, size=self.NEvents).astype('u1')
            g_hits.create_dataset('NHits', data=n_hits)
            x = round(normal(p.NCols / 2, beam_x / p.PX, sum(n_hits)))
            x[x < 0] = 0
            x[x >= p.NCols - 1] = p.NCols - 1
            y = round(normal(p.NRows / 2, beam_y / p.PY, sum(n_hits)))
            y[y < 0] = 0
            y[y >= p.NRows - 1] = p.NRows - 1
            g_hits.create_dataset('X', data=x, dtype='u2')
            g_hits.create_dataset('Y', data=y, dtype='u2')
            g_clu = g.create_group('Clusters')
