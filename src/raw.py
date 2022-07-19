#!/usr/bin/env python
# --------------------------------------------------------
#       raw file conversion & analysis
# created on April 22nd 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from subprocess import check_call
import h5py
from numpy import all, array

from plotting.draw import critical, Draw, info, remove_file, Path, choose, prep_kw
from plotting.fit import Gauss
from src.converter import Converter
from utility.affine_transformations import transform
from src.dut import Plane
from src.analysis import Analysis


class Raw:

    def __init__(self, c: Converter, load_file=False, step=-1):

        self.Parent = c
        self.Run = c.Run
        self.NT = c.NTelPlanes
        self.ND = c.NDUTPlanes
        self.P = range(self.NT)
        self.Planes = [Plane(i, typ='TELESCOPE' if i < self.NT else 'DUT') for i in range(self.NT + self.ND)]

        self.SoftDir = self.load_soft_dir()
        self.DataDir = c.DataDir
        self.SaveDir = c.SaveDir

        self.RawFilePath = self.load_raw_file_path()
        self.OutFilePath = self.load_out_file_path()

        if load_file:
            self.F = c.F if c.F is not None else h5py.File(c.OutFilePath, 'r')

        self.Steps = [(self.convert, self.OutFilePath)]
        self.AtStep = step
        self.Draw = Draw(Analysis.Config.FilePath)

    def __repr__(self):
        return f'{self.__class__.__name__} file analysis run {self.Run} ({self.RawFilePath.name})'

    # ----------------------------------------
    # region CONVERT
    def load_raw_file_path(self):
        n = list(self.DataDir.joinpath('raw').glob(f'run{self.Run:06d}*.raw'))
        return n[0] if len(n) else None

    def load_out_file_path(self):
        return self.SaveDir.joinpath(f'run{self.Run:06d}.root')

    def load_soft_dir(self):
        return self.Parent.SoftDir.joinpath(Analysis.Config.get('SOFTWARE', 'eudaq2'))

    def generate_fit_files(self):
        c = None
        for dut in range(self.Parent.NDUTPlanes):
            c = self.Parent.load_calibration(dut)
            if not c.FitFileName.exists():
                c.save_fit_pars()
        return c.CalPath

    @property
    def soft(self):
        return self.SoftDir.joinpath('bin', 'euCliConverter')

    @property
    def options(self):
        return f'-c {self.generate_fit_files()}'

    def convert(self):
        """ convert binary raw file to root file with eudaq"""
        if self.RawFilePath is None:
            critical(f'raw file does not exist for run: {self.Run}')
        self.OutFilePath.parent.mkdir(exist_ok=True)
        cmd = f'{self.soft} -i {self.RawFilePath} -o {self.OutFilePath} {self.options}'
        info(f'Convert {self.RawFilePath.name} to {self.OutFilePath.name} using {self.soft.name}\n')
        info(f'{cmd}\n')
        check_call(cmd, shell=True)
        for f in Path().glob('AutoDict_vector*'):
            remove_file(f)
    # endregion CONVERT
    # ----------------------------------------

    # ----------------------------------------
    # region ANALYSIS
    def get(self, p, g, k):
        return array(self.F[f'Plane{p}'][g][k])

    def z(self):
        return self.Parent.Proteus.get_z_positions(raw=True)[:self.NT]

    def l2g(self, x, y, p, step=None):
        p = self.Planes[p]
        a = self.Parent.Proteus.get_alignment(choose(step, self.AtStep))['sensors'][p.Number]
        ox, oy = array(a['offset'][:2])
        rx, ry = array(a['unit_u']), a['unit_v']
        return transform(x, y, sx=p.PX, sy=p.PY, ox=ox, oy=oy, rx=rx, ry=ry, order='trs')

    def xy(self, p, cut=...):
        return [self.get(p, 'Clusters', k)[cut] for k in ['X', 'Y']]

    def txy(self, p, cut=...):
        return [self.get(p, 'Tracks', k)[cut] for k in ['X', 'Y']]

    def uv(self, p, cut=..., step=None):
        return self.l2g(*self.xy(p, cut), p, step)

    def uvs(self, step=None):
        return array([self.uv(p, c, step) for p, c in zip(self.P, self.cuts())])

    def tuv(self, p, cut=..., step=None):
        return self.l2g(*self.txy(p, cut), p, step)

    def tuvs(self, step=None):
        return array([self.tuv(p, self.cut(), step) for p in self.P])

    def du(self, p, step=None):
        c = self.all_planes()
        return self.uv(p, self.cuts(c)[p], step)[0] - self.tuv(p, self.cut(c), step)[0]

    def all_planes(self):
        return array(self.F['Tracks']['Size']) == self.NT

    def cuts(self, cut=True):
        n = [self.get(p, 'Clusters', 'Size') > 0 for p in self.P]
        c = all(n, axis=0) & cut
        return [c[i] for i in n]

    def res(self):
        return array([Gauss(self.draw_du(p, show=False)).fit(draw=False).get_pars(err=False)[1:] for p in self.P])

    def cut(self, cut=True):
        return all([self.get(p, 'Clusters', 'Size') > 0 for p in self.P], axis=0) & cut

    def draw_clu_u(self, i=0, step=None):
        self.Draw.graph(self.z(), self.uvs(step)[:, 0, :].T[i], x_tit='Z [mm]', y_tit='Cluster U [mm]')

    def draw_track_u(self, i=0, step=None):
        self.Draw.graph(self.z(), self.tuvs(step)[:, 0, :].T[i], x_tit='Z [mm]', y_tit='Track U [mm]')

    def draw_du(self, p=0, step=None, **dkw):
        return self.Draw.distribution(self.du(p, step) * 1e3, **prep_kw(dkw, x_tit='dU [#mum]'))

    def draw_clu_y(self, i=0, step=None):
        self.Draw.graph(self.z(),  self.uvs(step)[:, 1, :].T[i], x_tit='Z [mm]', y_tit='V [mm]')
    # ----------------------------------------
    # endregion ANALYSIS


if __name__ == '__main__':

    c_ = Converter.from_ana(41, 1)
    z = Raw(c_, load_file=True, step=2)
