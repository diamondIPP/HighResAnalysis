#!/usr/bin/env python
# --------------------------------------------------------
#       adds clustering and charge to trees created with pXar
# created on August 30th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from ROOT import TFile, vector, TF1
from utils import *
from os.path import basename, join, dirname
from argparse import ArgumentParser
from collections import OrderedDict
from numpy import array, zeros
from progressbar import Bar, ETA, FileTransferSpeed, Percentage, ProgressBar
from pickle import load, dump
from draw import Draw, ufloat
from glob import glob
from converter import Converter, Hit, Cluster


class CERNConverter(Converter):

    def __init__(self, filename, plane, config, save_dir=None, first_run=None):

        Converter.__init__(self, filename, config)
        self.SaveDir = save_dir
        self.CalibrationDir = join(self.SaveDir, 'calibrations')
        self.RunNumber = basename(filename).strip('.root').split('_')[-1]
        self.FirstRun = first_run if first_run is not None else self.RunNumber

        self.OldFile = read_root_file(filename)
        self.OldTree = self.OldFile.Get('Plane{}'.format(6 + plane)).Get('Hits')
        self.EventTree = self.OldFile.Get('Event')
        self.NewFile = self.create_new_file(filename, save_dir)
        self.NewTree = self.OldTree.CloneTree(0)

        # New Branches
        self.ScalarBranches = OrderedDict([('NCluster', array([0], 'u2')),
                                           ('TimeStamp', array([0], 'f8'))])
        self.VectorBranches = OrderedDict([('VCal', vector('float')()),
                                           ('ClusterSize', vector('unsigned short')()),
                                           ('ClusterX', vector('float')()),
                                           ('ClusterY', vector('float')()),
                                           ('ClusterVcal', vector('float')())])

        # Charge Fits
        self.FitParameters = zeros((52, 80, 4))
        self.Fit = TF1('ErFit', '[3] * (TMath::Erf((x - [0]) / [1]) + [2])', -500, 255 * 7)
        self.get_fits()

        # Vars
        self.Hits = []
        self.ClusteredHits = []
        self.Clusters = []

        self.ProgressBar = None

    def start_pbar(self, n):
        self.ProgressBar = ProgressBar(widgets=['Progress: ', Percentage(), ' ', Bar(marker='>'), ' ', ETA(), ' ', FileTransferSpeed()], maxval=n)
        self.ProgressBar.start()

    def create_new_file(self, filename, save_dir):
        if save_dir is None:
            return TFile(join(dirname(filename), 'Clustered_{}.root'.format(self.RunNumber)), 'RECREATE')
        return TFile(join(save_dir, 'Clustered_{}.root'.format(self.RunNumber)), 'RECREATE')

    def set_branches(self):
        for key, value in self.ScalarBranches.iteritems():
            self.NewTree.Branch(key, value, '{}/{}'.format(key, type_dict[value[0].dtype.name]))
        for key, vec in self.VectorBranches.iteritems():
            self.NewTree.Branch(key, vec)

    def clear_vectors(self):
        for key in self.VectorBranches.iterkeys():
            self.VectorBranches[key].clear()
        self.Clusters = []
        self.ClusteredHits = []
        self.Hits = []

    def get_calibration_number(self):
        numbers = sorted(int(remove_letters(basename(name))) for name in glob(join(self.SaveDir, 'calibrations', 'phCal*.dat')) if basename(name)[5].isdigit())
        if self.FirstRun > numbers[-1]:
            return numbers[-1]
        next_number = next(nr for nr in numbers if nr >= int(self.FirstRun))
        return numbers[numbers.index(next_number) - 1]

    def get_fits(self):
        pickle_name = join(self.CalibrationDir, 'fitpars{}.pickle'.format(self.get_calibration_number()))
        if file_exists(pickle_name):
            with open(pickle_name, 'r') as f:
                self.FitParameters = load(f)
        elif file_exists(join(self.SaveDir, 'phCalibrationFits.dat')):
            self.get_fits_from_file()
        else:
            self.get_fits_from_calibration()
        with open(pickle_name, 'w') as f:
            dump(self.FitParameters, f)

    def get_fits_from_file(self):
        # TODO update
        with open(join(self.SaveDir, 'phCalibrationFits.dat')) as f:
            lines = f.readlines()[3:]
            for line in lines:
                words = line.split('Pix')
                col, row = [int(val) for val in words[-1].split()]
                data = [float(val) for val in words[0].split()]
                self.FitParameters[col][row] = data

    def get_fits_from_calibration(self):
        f = open(join(self.CalibrationDir, 'phCal{}.dat'.format(self.get_calibration_number())))
        f.readline()
        low_range = [int(val) for val in f.readline().split(':')[-1].split()]
        high_range = [int(val) for val in f.readline().split(':')[-1].split()]
        x = low_range + [val * 7 for val in high_range]
        f.readline()
        self.Fit.SetParameters(309.2062, 112.8961, 1.022439, 35.89524)
        d = Draw()
        self.start_pbar(52 * 80)
        for i, line in enumerate(f.readlines()):
            data = line.split('Pix')
            y = [int(val) for val in data[0].split()]
            x1 = [ufloat(ix, 1) for (ix, iy) in zip(x, y) if iy]
            y1 = [ufloat(iy, 1) for iy in y if iy]
            col, row = [int(val) for val in data[-1].split()]
            if not x1:
                self.FitParameters[col][row] = [0, 0, 0, 0]
                self.ProgressBar.update(i + 1)
                continue
            g = d.make_tgrapherrors('gcal', 'gcal', x=x1, y=y1)
            g.Fit(self.Fit, 'q', '', 0, 3000)
            self.FitParameters[col][row] = [self.Fit.GetParameter(j) for j in xrange(4)]
            self.ProgressBar.update(i + 1)
        self.ProgressBar.finish()

    def get_charge(self, col, row, adc):
        self.Fit.SetParameters(*self.FitParameters[col][row])
        return self.Fit.GetX(adc)

    def clusterise(self):
        if len(self.Hits) > 50:
            return
        for hit in self.Hits:
            if hit in self.ClusteredHits:
                continue
            cluster = Cluster()
            self.Clusters.append(cluster)
            cluster.add_hit(hit)
            self.ClusteredHits.append(hit)
            self.add_touching_hits(cluster, hit)

    def add_touching_hits(self, cluster, hit):
        for ihit in self.Hits:
            if ihit in self.ClusteredHits:
                continue
            if abs(ihit.X - hit.X) <= 1 and abs(ihit.Y - hit.Y) <= 1:
                cluster.add_hit(ihit)
                self.ClusteredHits.append(ihit)
                self.add_touching_hits(cluster, ihit)

    def run(self):
        self.set_branches()
        n = self.OldTree.GetEntries()
        self.start_pbar(n)
        for i, event in enumerate(self.OldTree):
            self.EventTree.GetEntry(i)
            self.clear_vectors()
            x, y, adc = event.PixX, event.PixY, event.Value
            for ix, iy, iadc in zip(x, y, adc):
                hit = Hit(ix, iy)
                hit.set_charge(self.get_charge(ix, iy, iadc))
                self.VectorBranches['VCal'].push_back(hit.Charge)
                self.Hits.append(hit)
            self.clusterise()
            self.ScalarBranches['NCluster'][0] = len(self.Clusters)
            self.ScalarBranches['TimeStamp'][0] = self.EventTree.TimeStamp
            for cluster in self.Clusters:
                self.VectorBranches['ClusterSize'].push_back(cluster.size())
                self.VectorBranches['ClusterX'].push_back(cluster.x())
                self.VectorBranches['ClusterY'].push_back(cluster.y())
                self.VectorBranches['ClusterVcal'].push_back(cluster.charge())
            self.NewTree.Fill()
            if not i % 100:
                self.ProgressBar.update(i + 1)
        self.ProgressBar.finish()
        self.NewFile.cd()
        self.NewFile.Write()
        self.NewFile.Close()


if __name__ == '__main__':
    p = ArgumentParser()
    p.add_argument('filename', nargs='?', default='')
    # noinspection PyTypeChecker
    p.add_argument('plane', nargs='?', default=0, type=int)
    args = p.parse_args()
    z = CERNConverter('/scratch2/cern/2018-10/cms-raw/ljutel_110.root', 0, '/scratch2/cern/2018-10/II6-B6')
    # z.run()
