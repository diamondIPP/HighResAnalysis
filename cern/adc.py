#!/usr/bin/env python
""" --------------------------------------------------------
       adc converter for CERN DUT data
 created on April 22nd 2022 by M. Reichmann (remichae@phys.ethz.ch)
 --------------------------------------------------------"""
import awkward as aw
import uproot
from numpy import round

from plotting.utils import info, array
from cern.converter import CERNConverter
from utility.utils import PBAR


class Adc2Vcal:

    def __init__(self, c: CERNConverter):

        self.Parent = c
        self.Run = c.Run
        self.NTelPlanes = c.NTelPlanes

        self.DataDir = c.DataDir.joinpath('dut')
        self.RawFilePath = self.DataDir.joinpath(f'ljutel_{c.Run:03d}.root')
        self.OutFilePath = self.RawFilePath.with_name(f'dut-run{c.Run:04d}.root')

        self.FlatBranches = ['NHits', 'Timing', 'TriggerCount']

        self.Steps = [(self.convert, self.OutFilePath)]

    def check_calibrations(self, n):
        return all([self.Parent.load_calibration(i).RawFileName.exists() for i in range(n)])

    def convert(self):
        """read the DUT raw file and convert the adc values to vcal"""
        if not self.RawFilePath.exists():
            CERNConverter.download_raw_file(self.RawFilePath)
        duts = array(self.Run.Logs['duts'])
        self.check_calibrations(duts.size)
        cut = self.Parent.EventAlignment.cut  # additional events in the DUT data ...
        with uproot.open(self.RawFilePath) as f:
            with uproot.recreate(self.OutFilePath) as g:
                g['Event'] = f['Event'].arrays()[cut]
                info(f'converting adc to vcal for DUTs {", ".join(duts)} ...')
                PBAR.start(duts.size)
                for i in range(duts.size):  # the data file only has no FEI4 plane
                    dir_name = f'Plane{self.NTelPlanes + i}'
                    x, y, adc = f[f'{dir_name}/Hits'].arrays(['PixX', 'PixY', 'Value'], library='np').values()
                    lut = self.Parent.load_calibration(i).get_lookup_table()
                    vcal = [round([lut[lx[i], ly[i], int(lz[i])] for i in range(lx.size)]).astype('i') for lx, ly, lz in zip(x, y, adc)]
                    data = f[f'{dir_name}/Hits'].arrays(filter_name=lambda w: not any([b in w for b in self.FlatBranches]))
                    data['Value'] = aw.values_astype(vcal, 'int32')
                    d = g.mkdir(dir_name)
                    flat = {n: arr.flatten()[cut] for n, arr in f[f'{dir_name}/Hits'].arrays(self.FlatBranches, library='np').items()}
                    flat['Timing'] = flat['Timing'].astype('i')
                    d['Hits'] = {**{k: v.flatten() for k, v in flat.items()}, 'TriggerPhase': flat['Timing'].astype('u2'), '': aw.zip({n: data[n][cut] for n in data.fields})}
                    PBAR.update()
        info(f'successfully wrote {self.OutFilePath} ({flat["Timing"].size} ev)')
