#!/usr/bin/env python
# --------------------------------------------------------
#       event alignment
# created on July 26th 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from src.proteus import Proteus


class EventAlignment:

    def __init__(self, p: Proteus):
        self.Proteus = p
        self.FileName = self.load_file_name()

    def __repr__(self):
        return f'{self.__class__.__name__} of Run {self.Proteus.RunNumber}'

    def load_file_name(self):
        return self.Proteus.TrackName.with_name(f'{self.Proteus.TrackName.name}-data.root')  # automatic ending from Proteus
