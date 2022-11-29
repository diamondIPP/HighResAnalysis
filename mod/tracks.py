#!/usr/bin/env python
# --------------------------------------------------------
#       Track analysis of the high resolution tests at DESY and CERN
# created on June 10th 2020 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from numpy import rad2deg, where

import src.bins as bins
from mod.track_cuts import TrackCut
from plotting.draw import prep_kw
from utility.utils import uarr2n


def track_analysis(cls):
    class TrackAnalysis(cls):
        M = ['x', 'y']

        def __init__(self, parent):  # noqa

            self.Parent = parent
            self.__dict__.update(parent.__dict__)
            self.MetaSubDir = 'tracks'
            self.Tel.T = True

            self.Cut = TrackCut(self)

        def ax_tits(self, local=True):
            return super().ax_tits(local) if local else {'x_tit': 'Track X [mm]', 'y_tit': 'Track Y [mm]'}

        # ----------------------------------------
        # region DATA
        def get_n(self, cut=None):
            """returns: number of tracks in trk-space"""
            return self.get_data('N', cut=cut, main_grp='Tracks')

        def get_size(self, cut=None):
            """returns: number of planes per track"""
            return self.get_data('Size', cut=cut, main_grp='Tracks')

        def get_off_x(self, cut=None):  # unclear what this exactly is...
            return self.get_data('X', cut=cut, main_grp='Tracks')

        def get_off_y(self, cut=None):
            return self.get_data('Y', cut=cut, main_grp='Tracks')

        def get_xy(self, local=True, cut=None, pl=None, centre=False):
            return self.get_txy(local, cut, pl, centre)
        # endregion DATA
        # ----------------------------------------

        # ----------------------------------------
        # region MISC
        def get_offsets(self, redo=False):
            return uarr2n(self.Residuals.means(local=True, _redo=redo)).reshape(-1, 1)

        def get_z(self, raw=False):
            """returns: z-pos of the telescope planes in [cm]"""
            return self.Proteus.z_positions(raw)[:self.Tel.NPlanes] / 10
        # endregion MISC
        # ----------------------------------------

        # ----------------------------------------
        # region DRAW
        def draw_x(self, trk=0, **dkw):
            """draw track in x-z plane."""
            trk = where(self.Tel.Cut.make_all_cluster())[0][trk]
            x, y = self.get_z(), self.Tel.get_us(self.Cut.make_trk(trk)).reshape(-1)
            self.Draw.graph(x, y, **prep_kw(dkw, x_tit='Z Position [cm]', y_tit='X Position [mm]'))

        def draw_y(self, trk=0, **dkw):
            """draw track in y-z plane."""
            trk = where(self.Tel.Cut.make_all_cluster())[0][trk]
            x, y = self.get_z(), self.Tel.get_vs(self.Cut.make_trk(trk)).reshape(-1)
            self.Draw.graph(x, y, **prep_kw(dkw, x_tit='Z Position [cm]', y_tit='Y Position [mm]'))

        def draw_n(self, cut=None, **dkw):
            return self.Draw.distribution(self.get_n(cut), **prep_kw(dkw, title='NTracks', x0=-.5, w=1, x_tit='Number of Tracks'))

        def draw_map(self, bw=.3, local=True, cut=None, dut_plane=True, **dkw):
            binning = bins.get_xy(local, self.Plane if dut_plane else self.Tel.Plane, bw, aspect_ratio=True)
            return self.Draw.histo_2d(*self.get_xy(local, cut), binning, **prep_kw(dkw, title='Track Map', **self.ax_tits(local)))

        def draw_chi2(self, cut=None, **dkw):
            return self.Draw.distribution(self.get_chi2(cut), **prep_kw(dkw, x0=0, title='Chi2', x_tit='#chi^{2}'))

        def draw_slopes(self, cut=None, **dkw):
            g = [getattr(self, f'draw_slope_{m}')(cut, lf=1, rf=1, show=False) for m in self.M]
            return self.Draw.stack(g, 'TrackSlopes', self.M, **prep_kw(dkw))

        def draw_slope_x(self, cut=None, **dkw):
            return self.Draw.distribution(rad2deg(self.get_slope_x(cut)), **prep_kw(dkw, title='Track Slope X', x_tit='Track Slope [deg]'))

        def draw_slope_y(self, cut=None, **dkw):
            return self.Draw.distribution(rad2deg(self.get_slope_y(cut)), **prep_kw(dkw, title='Track Slope Y', x_tit='Track Slope [deg]'))
        # endregion DRAW
        # ----------------------------------------

    return TrackAnalysis
