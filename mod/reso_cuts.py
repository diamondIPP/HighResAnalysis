#!/usr/bin/env python
# --------------------------------------------------------
#       cuts for the reference detector
# created on March 30th 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from mod.dut_cuts import DUTCut, array
from copy import deepcopy


class ResCut(DUTCut):
    def __init__(self, parent):  # noqa

        self.__dict__.update(parent.__dict__)
        self.Cuts = deepcopy(parent.Cuts)
        self.register('cs', self.make_cs(), 65, 'cluster size = 1')

    def make_cross(self, col, row):
        x, y = self.Ana.get_xy(cut=False)
        return self.add(array((x == col) | (y == row)))

    def make_col(self, col):
        return self.add(array(self.Ana.get_x(cut=False) == col))

    def make_row(self, row):
        return self.add(array(self.Ana.get_y(cut=False) == row))
