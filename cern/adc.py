#!/usr/bin/env python
""" --------------------------------------------------------
       adc converter for CERN DUT data
 created on April 22nd 2022 by M. Reichmann (remichae@phys.ethz.ch)
 --------------------------------------------------------"""
from src.converter import Converter
from plotting.utils import info, array
from utility.utils import PBAR
import uproot
import awkward as aw


class Adc2Vcal:

    def __init__(self, c: Converter):

        self.Parent = c
        self.Run = c.Run
        self.NTelPlanes = c.NTelPlanes

        self.DataDir = c.DataDir.joinpath('dut')
        self.RawFilePath = self.DataDir.joinpath(f'ljutel_{c.Run:03d}.root')
        self.OutFilePath = self.RawFilePath.with_name(f'dut-run{c.Run:04d}.root')

        self.Steps = [(self.convert, self.OutFilePath)]

    def convert(self):
        """read the DUT raw file and convert the adc values to vcal"""
        with uproot.open(self.RawFilePath) as f:
            with uproot.recreate(self.OutFilePath) as g:
                g['Event'] = f['Event'].arrays()
                dut_nrs = array([i for i, dut in enumerate(self.Run.Logs['duts']) if 'FEI4' not in dut])
                info(f'converting adc to vcal for DUTs {", ".join(self.Run.Logs["duts"][i] for i in dut_nrs)} ...')
                PBAR.start(dut_nrs.size)
                for i in range(dut_nrs.size):  # the data file only has no FEI4 plane
                    dir_name = f'Plane{self.NTelPlanes + i}'
                    x, y, adc = f[f'{dir_name}/Hits'].arrays(['PixX', 'PixY', 'Value'], library='np').values()
                    lut = self.Parent.load_calibration(dut_nrs[i]).get_lookup_table()
                    vcal = [[lut[lx[i], ly[i], int(lz[i])] for i in range(lx.size)] for lx, ly, lz in zip(x, y, adc)]
                    data = f[f'{dir_name}/Hits'].arrays(filter_name=lambda w: 'NHits' not in w)
                    data['Value'] = vcal
                    d = g.mkdir(dir_name)
                    d['Hits'] = {'NHits': array(f[f'{dir_name}/Hits']['NHits']), '': aw.zip({n: data[n] for n in data.fields})}
                    PBAR.update()
        info(f'successfully wrote {self.OutFilePath} ({x.size} ev)')
