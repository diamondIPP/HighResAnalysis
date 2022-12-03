#!/usr/bin/env python
# --------------------------------------------------------
#       classes for actions on several runs or batches
# created on November 30th 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from src.run import Ensemble
from src.batch_analysis import DUTAnalysis, BatchAnalysis, Batch
from numpy import array
from plotting.save import prep_kw, SaveDraw, Path


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

    def values(self, f, *args, **kwargs):
        return array([f(ana, *args, **kwargs) for ana in self.Anas])

    def t(self):
        return self.values(DUTAnalysis.mean_time)

    def x(self):
        return self.t()

    def init_analyses(self, verbose, test):
        return [BatchAnalysis.from_batch(u, verbose, test) if type(u) is Batch else DUTAnalysis.from_run(u, verbose, test) for u in self.Units]

    def draw_graph(self, y, t=False, **dkw):
        x = self.t() if t else self.x()
        fname = f'{dkw.pop("file_name")}{"Time" if t else ""}' if 'file_name' in dkw else None
        return self.Draw.graph(x, y, **prep_kw(dkw, **Scan.XArgs if t else self.XArgs, file_name=fname))

    def draw_efficiency(self, t=False, **dkw):
        self.draw_graph(self.values(DUTAnalysis.eff), t, **prep_kw(dkw, y_tit='Efficiency [%]', file_name='Eff'))

    def draw_current(self, t=False, **dkw):
        self.draw_graph(self.values(DUTAnalysis.current), t, **prep_kw(dkw, y_tit='Current [nA]', file_name='Curr'))

    def draw_pulse_height(self, t=False, **dkw):
        self.draw_graph(self.values(DUTAnalysis.ph), t, **prep_kw(dkw, y_tit='Pulse Height [vcal]', file_name='PH'))


class VScan(Scan):

    XArgs = {'x_tit': 'Bias [V]'}

    def __init__(self, name, verbose=False, test=False):
        super().__init__(name, verbose, test)

    def x(self):
        return array(self.biases)


class TScan(Scan):

    XArgs = {'x_tit': 'Trim [vcal]'}

    def __init__(self, name, verbose=False, test=False):
        super().__init__(name, verbose, test)

    def x(self):
        return array([ana.Calibration.Trim for ana in self.Anas])


if __name__ == '__main__':

    z = VScan('v-b2')
