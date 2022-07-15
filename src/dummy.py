#!/usr/bin/env python
# --------------------------------------------------------
#       module to create dummy file for analysis
# created on July 29th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT.gRandom import Landau
from numpy import round, invert, delete
from numpy.random import poisson, normal, randint

from src.dut import Plane
from src.analysis import Analysis
from utility.utils import *


class Dummy:

    def __init__(self, tc_dir, n_tel, n_dut):

        self.FileName = join(tc_dir, 'data', 'dummy.hdf5')
        self.File = None
        self.NTelPlanes = n_tel
        self.NDUTPlanes = n_dut
        self.TelPlane = Plane(0, typ='TELESCOPE')
        self.DUTPlane = Plane(self.NTelPlanes, typ='DUT')
        self.NEvents = 0
        self.PBar = PBar()

    def init_file(self):
        self.File = h5py.File(self.FileName, 'w')

    def load_file(self):
        return h5py.File(self.FileName, 'r')

    def set_n_events(self, n):
        self.NEvents = int(n)

    def make(self, n_events=1e6, beam_x=.8, beam_y=.8):
        self.set_n_events(n_events)
        remove_file(self.FileName)
        self.init_file()
        self.add_tracks(beam_x, beam_y)
        self.add_time()
        self.add_telescope(beam_x, beam_y)
        self.add_dut()

    def add_tracks(self, beam_x, beam_y):
        g = self.File.create_group('Tracks')
        n_tracks = poisson(lam=1.5, size=self.NEvents)
        n = sum(n_tracks)
        # create a random origin in DUT Plane
        c_x, c_y = self.dut_2_tel(x=normal(self.DUTPlane.NCols / 2, beam_x / self.DUTPlane.PX), y=normal(self.DUTPlane.NRows / 2, beam_y / self.DUTPlane.PY))
        x = round(normal(c_x, beam_x / self.TelPlane.PX, n))
        y = round(normal(c_y, beam_y / self.TelPlane.PY, n))
        cut = (x >= 0) & (x < self.TelPlane.NCols) & (y >= 0) & (y < self.TelPlane.NRows)
        s = x[cut].size
        while sum(n_tracks) > s:
            n_tracks = delete(n_tracks, -1)
        g.create_dataset('NTracks', data=n_tracks, dtype='u1')
        g.create_dataset('X', data=x[cut], dtype='f2')
        g.create_dataset('Y', data=y[cut], dtype='f2')

    def add_time(self):
        g = self.File.create_group('Event')
        g.create_dataset('Time', data=arange(0, 1e5, .00165, dtype='f2')[:self.NEvents])

    def add_telescope(self, beam_x, beam_y):
        self.PBar.start(self.NTelPlanes)
        for i in range(self.NTelPlanes):
            g = self.File.create_group('Plane{}'.format(i))
            g_hits = g.create_group('Hits')
            n_hits = poisson(lam=6, size=self.NEvents).astype('u1')
            g_hits.create_dataset('NHits', data=n_hits)
            x = round(normal(self.TelPlane.NCols / 2, beam_x / self.TelPlane.PX, sum(n_hits)))
            x[x < 0] = 0
            x[x >= self.TelPlane.NCols - 1] = self.TelPlane.NCols - 1
            y = round(normal(self.TelPlane.NRows / 2, beam_y / self.TelPlane.PY, sum(n_hits)))
            y[y < 0] = 0
            y[y >= self.TelPlane.NRows - 1] = self.TelPlane.NRows - 1
            g_hits.create_dataset('X', data=x, dtype='u2')
            g_hits.create_dataset('Y', data=y, dtype='u2')
            g_clu = g.create_group('Clusters')
            g_clu.create_dataset('X', data=x, dtype='u2')
            g_clu.create_dataset('Y', data=y, dtype='u2')
            g_clu.create_dataset('NClusters', data=n_hits)
            self.PBar.update()

    def add_dut(self):
        lam = array([.077, .15])
        c = [[450, 40], [200, 25]]
        for i in range(self.NTelPlanes, self.NDUTPlanes + self.NTelPlanes):
            dut_plane = Plane(i, Analysis.Config('DUT'))
            g = self.File.create_group('Plane{}'.format(i))
            g_clu = g.create_group('Clusters')
            size = poisson(lam[[i - self.NTelPlanes]], self.File['Tracks']['X'].size)
            a_cut = size > 0
            size_g0 = size[a_cut]
            x_tra, y_tra = array(self.File['Tracks']['X']), array(self.File['Tracks']['Y'])
            g.create_dataset('TriggerPhase', data=randint(0, 10, x_tra.size), dtype='f2')
            x, y = self.tel_2_dut(x_tra, y_tra)
            x, y = x[a_cut], y[a_cut]
            cut = invert((x >= 0) & (x < dut_plane.NCols) & (y >= 0) & (y < dut_plane.NRows))
            size_g0[cut] = 50
            size[size > 0] = size_g0
            x[cut] = -20
            y[cut] = -20
            g_clu.create_dataset('Size', data=size, dtype='u1')
            m, s = c[i - self.NTelPlanes]
            g_clu.create_dataset('Charge', data=array([Landau(m, s) for _ in range(x.size)]), dtype='f2')
            g_clu.create_dataset('X', data=x, dtype='f2')
            g_clu.create_dataset('Y', data=y, dtype='f2')
            u, v = self.l2g(x, y, is_dut=True)
            g_clu.create_dataset('U', data=u, dtype='f2')
            g_clu.create_dataset('V', data=v, dtype='f2')
            g_tra = g.create_group('Tracks')
            x_tra_clu, y_tra_clu = x_tra + normal(0, .02 / dut_plane.PX, x_tra.size), y_tra + normal(0, .02 / dut_plane.PY, x_tra.size)
            u_tra_clu, v_tra_clu = self.l2g(x_tra_clu, y_tra_clu)
            g_tra.create_dataset('X', data=x_tra_clu, dtype='f2')
            g_tra.create_dataset('Y', data=y_tra_clu, dtype='f2')
            g_tra.create_dataset('U', data=u_tra_clu, dtype='f2')
            g_tra.create_dataset('V', data=v_tra_clu, dtype='f2')

    def l2g(self, x, y, is_dut=False):
        x, y = self.dut_2_tel(x, y) if is_dut else (x, y)
        return (x + .5 - self.TelPlane.NCols / 2) * self.TelPlane.PX, (y + .5 - self.TelPlane.NRows / 2) * self.TelPlane.PY

    def dut_2_tel(self, x, y):
        """ assume that dut is in the centre of tel"""
        return x + self.TelPlane.NCols / 2 - self.DUTPlane.NCols / 2, y + self.TelPlane.NRows / 2 - self.DUTPlane.NRows / 2

    def tel_2_dut(self, x, y):
        """ assume that dut is in the centre of tel"""
        return x - self.TelPlane.NCols / 2 + self.DUTPlane.NCols / 2, y - self.TelPlane.NRows / 2 + self.DUTPlane.NRows / 2
