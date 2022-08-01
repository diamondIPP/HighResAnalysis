#!/usr/bin/env python
# --------------------------------------------------------
#       event alignment
# created on July 26th 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------
from plotting.utils import info, choose, critical
from cern.raw import Raw
import uproot
from numpy import array, roll, where, diff, abs, delete, ones


class EventAlignment:

    def __init__(self, raw: Raw):
        self.Raw = raw
        self.RefPath = raw.OutFilePath.with_name(f'anchor{raw.Run.Logs["telescope run"]}.root')

        self.X = array([])      # ts REF
        self.Y = array([])      # ts TEL
        self.Cut = array([])    # invalid events

        self.OffEvents = []

    def __repr__(self):
        return f'{self.__class__.__name__} of {self.Raw.Run!r}'

    def run(self):
        self.load_data()
        if not self.is_good and len(self.OffEvents) == 0:
            self.find_events()
            if not self.validated:
                critical(f'{self!r} failed!')

    @property
    def cut(self):
        self.run()
        c = ones(self.X.size + 1, '?')
        c[self.OffEvents] = False
        return c

    @property
    def is_good(self):
        return self.off_events().size == 0

    def load_data(self, reload=False):
        if not self.X.size or reload:
            with uproot.open(self.RefPath) as f:
                a = array(f['Event']['TimeStamp'], 'i8')
            with uproot.open(self.Raw.OutFilePath) as f:
                b = array(f['Event']['TimeStamp'], 'i8')
                c = array(f['Event']['Invalid'], '?')
            a, b, c = diff(a), diff(b // 2), ~(c[:-1] | roll(c[:-1], -1))  # always two timestamps are corrupt
            a[a < 0] += 2 ** 31  # fix resets of 32bit int
            b[(b < 0) & c] += 2 ** 31
            self.X, self.Y, self.Cut = a, b, c

    def tdiff(self, start, off, x=None):
        """returns: time difference between two time stamps for a given offset [off] """
        return abs(roll(choose(x, self.X), -off)[start:self.Y.size] - self.Y[start:])

    def off_events(self, start=0, off=0, x=None):
        """:returns: events with deviating time stamps. """
        d = self.tdiff(start, off, x)
        return where((d > self.Y[start:] * .001) & self.Cut[start:])[0]

    def find_events(self, start=0, off=0):
        off_events = self.off_events(start, off)
        if off_events.size:
            e = off_events[0] + 1 + start
            self.OffEvents.append(e + off)
            return self.find_events(e, off + 1)
        info(f'found {len(self.OffEvents)} event offsets')

    @property
    def validated(self):
        return self.off_events(x=delete(self.X, self.OffEvents)).size == len(self.OffEvents)
