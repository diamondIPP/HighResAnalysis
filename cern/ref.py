""" --------------------------------------------------------
       onverter for FE-IV reference plane
 created on October 12th 2022 by M. Reichmann (remichae@phys.ethz.ch)
 --------------------------------------------------------"""
import awkward as aw
import uproot

from numpy import count_nonzero
from cern.converter import CERNConverter
from plotting.utils import info


class RefConverter:

    def __init__(self, c: CERNConverter):

        self.Parent = c

        self.RawFilePath = c.EventAlignment.RefPath
        self.OutFilePath = self.RawFilePath.with_name(f'ref-run{c.Run:04d}.root')

        self.FlatBranches = ['NHits']

        self.Steps = [(self.convert, self.OutFilePath)]

    def convert(self):
        """read the REF root file and remove misaligned events"""
        if not self.RawFilePath.exists():
            CERNConverter.download_raw_file(self.RawFilePath)
        cut = self.Parent.EventAlignment.cut  # additional events in the DUT data ...
        with uproot.open(self.RawFilePath) as f:
            with uproot.recreate(self.OutFilePath) as g:
                info(f'removing additional events for the REF plane')
                data = f[f'Plane0/Hits'].arrays(filter_name=lambda w: not any([b in w for b in self.FlatBranches]))
                d = g.mkdir(f'Plane{self.Parent.NTelPlanes}')
                flat = {n: arr.flatten()[cut] for n, arr in f[f'Plane0/Hits'].arrays(self.FlatBranches, library='np').items()}
                d['Hits'] = {**{k: v.flatten() for k, v in flat.items()}, '': aw.zip({n: data[n][cut] for n in data.fields})}
        info(f'successfully wrote {self.OutFilePath} ({count_nonzero(cut)} ev)')
