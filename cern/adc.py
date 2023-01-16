#!/usr/bin/env python
""" --------------------------------------------------------
       adc converter for CERN DUT data
 created on April 22nd 2022 by M. Reichmann (remichae@phys.ethz.ch)
 --------------------------------------------------------"""
import awkward as aw
import uproot
from numpy import round

from plotting.utils import info
from cern.converter import CERNConverter
from src.dut import DUT
from utility.utils import PBAR


class Adc2Vcal:

    def __init__(self, c: CERNConverter, dut_name=None):
        """Reads the DUT ROOT files saved with pXar and converts adc -> vcal and removes faulty events from the telescope.
           Specific DUTs can be selected with the <dut_name> arg [default = None -> all DUTs]. """

        self.Parent = c
        self.Run = c.Run
        self.DUTs = [DUT.from_name(n, self.Run.Info) for n in (self.Run.DUTs if dut_name is None else [dut_name])]
        self.NDUTs = len(self.DUTs)
        self.NTelPlanes = c.NTelPlanes

        self.DataDir = c.DataDir.joinpath('dut')
        self.RawFilePath = self.DataDir.joinpath(f'ljutel_{c.Run:03d}.root')
        self.OutFilePath = self.RawFilePath.with_name(f'dut-run{c.Run:04d}.root')

        self.FlatBranches = ['NHits', 'Timing', 'TriggerCount']

        self.Steps = [(self.convert, self.OutFilePath)]

    def check_calibrations(self):
        return all([self.Parent.load_calibration(dut.Number).RawFileName.exists() for dut in self.DUTs])

    def convert(self):
        """read the DUT raw file and convert the adc values to vcal"""
        CERNConverter.download_raw_file(self.RawFilePath)
        self.check_calibrations()
        cut = self.Parent.EventAlignment.cut  # additional events in the DUT data ...
        with uproot.open(self.RawFilePath) as f:
            with uproot.recreate(self.OutFilePath) as g:
                g['Event'] = f['Event'].arrays()[cut]
                info(f'converting adc to vcal for DUTs {", ".join(dut.Name for dut in self.DUTs)} ...')
                PBAR.start(self.NDUTs)
                for i, dut in enumerate(self.DUTs):
                    dir_name = f'Plane{self.NTelPlanes + dut.Number}'
                    x, y, adc = f[f'{dir_name}/Hits'].arrays(['PixX', 'PixY', 'Value'], library='np').values()
                    lut = self.Parent.load_calibration(dut.Number).get_lookup_table()
                    vcal = [round([lut[lx[i], ly[i], int(lz[i])] for i in range(lx.size)]).astype('i') for lx, ly, lz in zip(x, y, adc)]
                    data = f[f'{dir_name}/Hits'].arrays(filter_name=lambda w: not any([b in w for b in self.FlatBranches]))
                    data['Value'] = aw.values_astype(vcal, 'int32')
                    d = g.mkdir(f'Plane{self.NTelPlanes + i + 1}')  # first plane after TEL is REF
                    flat = {n: arr.flatten()[cut] for n, arr in f[f'{dir_name}/Hits'].arrays(self.FlatBranches, library='np').items()}
                    flat['Timing'] = flat['Timing'].astype('i')
                    d['Hits'] = {**{k: v.flatten() for k, v in flat.items()}, 'TriggerPhase': flat['Timing'].astype('u2'), '': aw.zip({n: data[n][cut] for n in data.fields})}
                    PBAR.update()
        info(f'successfully wrote {self.OutFilePath} ({flat["Timing"].size} ev)')
