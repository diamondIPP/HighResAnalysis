#!/usr/bin/env python
# --------------------------------------------------------
#       handles the cuts for the high rate analysis
# created on August 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from draw import Draw
from ROOT import TCut, TCutG
from utils import *
from os.path import join
from json import loads
from numpy import array


class Cut(Draw):
    """
    The Elementary class provides default behaviour objects in the analysis framework and is the Mother of all myPadAnalysis objects.
    It provides, among other things, a verbose printing method or a save plot method containing a global save directory handling.
    """

    def __init__(self, analysis):
        Draw.__init__(self)

        self.Ana = analysis
        self.Dir = analysis.Dir
        self.CutStrings = self.define_cutstrings()
        self.Config = load_config(join(self.Dir, 'config', 'cut'))

        self.NCuts = 0

        self.generate()

    def __call__(self):
        return self.generate_all_cut()

    def generate(self):
        self.CutStrings['fiducial'] += self.generate_fiducial()

    def generate_all_cut(self):
        cut = TCut('AllCuts', '')
        self.NCuts = 0
        for key, value in self.CutStrings.iteritems():
            if not key.startswith('old') and not key.startswith('AllCut'):
                cut += value
                self.NCuts += 1
        return cut

    @staticmethod
    def define_cutstrings():
        """ Defines the ordered dictionary that contains all the final cuts and the order they are going to be applied."""
        dic = OrderedDict()
        dic['raw'] = TCut('raw', '')
        # dic['cluster'] = TCut('cluster', 'NCluster == 1')
        dic['fiducial'] = TCut('fiducial', '')
        dic['AllCuts'] = TCut('AllCuts', '')
        return dic

    def generate_fiducial(self):
        name = 'fid{}'.format(self.Ana.RunNumber)
        xmin, xmax = loads(self.Config.get('CUT', 'fiducial'))[:2]
        ymin, ymax = loads(self.Config.get('CUT', 'fiducial'))[2:]
        xmin, ymin = xmin - .5, ymin - .5
        xmax, ymax = xmax + .5, ymax + .5
        x = array([xmin, xmin, xmax, xmax, xmin], 'd')
        y = array([ymin, ymax, ymax, ymin, ymin], 'd')
        cut = TCutG(name, 5, x, y)
        cut.SetVarX('ClusterX')
        cut.SetVarY('ClusterY')
        self.Objects.append(cut)
        cut.SetLineColor(2)
        cut.SetLineWidth(3)
        return name

    def draw_fid_cut(self):
        cut = get_object('fid{}'.format(self.Ana.RunNumber))
        cut.Draw()
