#!/usr/bin/env python
# --------------------------------------------------------
#       small script to start high resolution analysis
# created on March 18th 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from argparse import ArgumentParser
from numpy import *

import src.bins as bins  # noqa
from plotting.draw import *
from src.dut_analysis import DUTAnalysis, Analysis
from src.utils import *  # noqa

t_start = time()
aparser = ArgumentParser()
aparser.add_argument('run', nargs='?', default=11, type=int)
aparser.add_argument('dut', nargs='?', default=1, type=int)
aparser.add_argument('--testcampaign', '-tc', nargs='?', default=Analysis.find_testcampaign())
aparser.add_argument('--verbose', '-v', action='store_false')
aparser.add_argument('--single_mode', '-s', action='store_false')
aparser.add_argument('--test', '-t', action='store_true')
args = aparser.parse_args()

z = DUTAnalysis(args.run, args.dut, test_campaign=args.testcampaign, single_mode=args.single_mode, verbose=args.verbose, test=args.test)
z.add_info(t_start, prnt=True)

# aliases
try:
    dut = z.DUT
    r = z.Run
    c = z.Converter
    cut = z.Cut
    cal = z.Calibration
    res = z.Residuals
    tel = z.Tel
    ref = z.REF
    t = z.Tracks
    cur = z.Currents
    e = z.Efficiency
except AttributeError:
    pass
