#!/usr/bin/env python
# --------------------------------------------------------
#       classes for actions on several runs or batches
# created on November 30th 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from src.run import Ensemble, Dir
from src.batch_analysis import DUTAnalysis, BatchAnalysis, Batch
from numpy import array
from plotting.save import prep_kw, SaveDraw, Path, Draw, datetime, choose, rm_key
from utility.utils import bias2rootstr


class Scan(Ensemble):
    """Base class defining actions on several runs or batches"""

    XArgs = {'x_tit': 'Time', 't_ax_off': 0}

    def __init__(self, name, verbose=False, test=False):

        super().__init__(name)

        self.Anas = self.init_analyses(verbose, test)
        self.Draw = SaveDraw(self, results_dir=self.Name)

    @property
    def server_save_dir(self):
        return Path('duts', str(self.DUT), self.Name)

    @property
    def suffix(self):
        return f'{self.DUT.Name}-{self.__class__.__name__}-{self.Anas[0].BeamTest.Location}'.lower().replace('ii6-', '')

    def save_plots(self):
        old_dir = self.Draw.ResultsDir
        SaveDraw.SaveOnServer = False
        self.Draw.ResultsDir = Dir.joinpath('tmp')
        self.draw_current(fn=f'cur-{self.suffix}')
        self.draw_efficiency(fn=f'e-{self.suffix}')
        self.draw_pulse_height(fn=f'ph-{self.suffix}')
        self.Draw.ResultsDir = old_dir
        SaveDraw.SaveOnServer = True

    @property
    def n_ev(self):
        return sum([u.n_ev for u in self.Units])

    def values(self, f, *args, **kwargs):
        return array([f(ana, *args, **kwargs) for ana in self.Anas])

    def t(self):
        return self.values(DUTAnalysis.mean_time)

    def x(self):
        return self.t()

    def x2str(self):
        return [datetime.fromtimestamp(i).strftime('%H:%M') for i in self.x()]

    def legend(self, h, titles=None, **kwargs):
        return Draw.legend(h, choose(self.x2str, titles), **kwargs)

    def init_analyses(self, verbose, test):
        return [BatchAnalysis.from_batch(u, verbose, test) if type(u) is Batch else DUTAnalysis.from_run(u, verbose, test) for u in self.Units]

    def draw_graph(self, y, t=False, **dkw):
        x = self.t() if t else self.x()
        fname = f'{dkw.pop("file_name")}{"Time" if t else ""}' if 'file_name' in dkw else None
        return self.Draw.graph(x, y, **prep_kw(dkw, **Scan.XArgs if t else self.XArgs, file_name=fname))

    def draw_efficiency(self, t=False, **dkw):
        return self.draw_graph(self.values(DUTAnalysis.eff), t, **prep_kw(dkw, y_tit='Efficiency [%]', file_name='Eff'))

    def draw_current(self, t=False, **dkw):
        return self.draw_graph(self.values(DUTAnalysis.current), t, **prep_kw(dkw, y_tit='Current [nA]', file_name='Curr'))

    def draw_pulse_height(self, t=False, **dkw):
        return self.draw_graph(self.values(DUTAnalysis.ph), t, **prep_kw(dkw, y_tit='Pulse Height [vcal]', file_name='PH'))

    def draw_cluster_size(self, t=False, **dkw):
        return self.draw_graph(self.values(DUTAnalysis.cs), t, **prep_kw(dkw, y_tit='Cluster Size', file_name='CS'))

    def draw_r_ph_cols(self, t=False, r=7, **dkw):
        return self.draw_graph(self.values(DUTAnalysis.r_ph_cols, r=r), t, **prep_kw(dkw, y_tit='Pulse Height Ratio', file_name='RPHCols'))


    def draw_ph_dists(self, **dkw):
        h = [ana.draw_signal_distribution(save=False, **rm_key(dkw, 'save')) for ana in self.Anas]
        return self.Draw.stack(h, 'PHDists', self.x2str(), **prep_kw(dkw, scale=True, file_name='PhDists'))


class VScan(Scan):

    XArgs = {'x_tit': 'Bias [V]'}

    def __init__(self, name, verbose=False, test=False):
        super().__init__(name, verbose, test)

    def x(self):
        return array(self.biases)

    def x2str(self):
        return bias2rootstr(*self.x())


class TScan(Scan):

    XArgs = {'x_tit': 'Trim [vcal]'}

    def __init__(self, name, verbose=False, test=False):
        super().__init__(name, verbose, test)

    def x(self):
        return array([ana.Calibration.Trim for ana in self.Anas])

    def x2str(self):
        return [str(i) for i in self.x()]

    def draw_ph_dists(self, **dkw):
        return super().draw_ph_dists(**prep_kw(dkw, leg_head='Trim [vcal]'))


if __name__ == '__main__':

    z = VScan('v-b2')
