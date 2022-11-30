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

    XArgs = {'x_tit': 'Time'}

    def __init__(self, name, verbose=False, test=False):

        super().__init__(name)

        self.Anas = self.init_analyses(verbose, test)
        self.Draw = SaveDraw(self, results_dir=self.Name)

    @property
    def server_save_dir(self):
        return Path('duts', str(self.DUT), self.Name)

    def values(self, f, *args, **kwargs):
        return array([f(ana, *args, **kwargs) for ana in self.Anas])

    def x(self):
        return self.values(DUTAnalysis.mean_time)

    def init_analyses(self, verbose, test):
        return [BatchAnalysis.from_batch(u, verbose, test) if type(u) is Batch else DUTAnalysis.from_run(u, verbose, test) for u in self.Units]

    def draw_efficiency(self, **dkw):
        x, y = self.x(), self.values(DUTAnalysis.eff)
        self.Draw.graph(x, y, **prep_kw(dkw, **self.XArgs, y_tit='Efficiency [%]', file_name='Eff'))


class VScan(Scan):

    XArgs = {'x_tit': 'Bias [V]'}

    def __init__(self, name, verbose=False, test=False):
        super().__init__(name, verbose, test)

    def x(self):
        return array(self.biases)


if __name__ == '__main__':

    z = VScan('v-b2')
