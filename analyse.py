#!/usr/bin/env python
# --------------------------------------------------------
#       small script to start high resolution analysis
# created on March 18th 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from argparse import ArgumentParser
from numpy import *

import src.bins as bins  # noqa
from plotting.draw import *  # noqa
from src.dut_analysis import DUTAnalysis, Analysis
from utility.utils import *  # noqa

t_start = time()
aparser = ArgumentParser()
aparser.add_argument('run', nargs='?', default=232, type=int)
aparser.add_argument('dut', nargs='?', default=0, type=int)
aparser.add_argument('--testcampaign', '-tc', nargs='?', default=Analysis.find_testcampaign())
aparser.add_argument('--verbose', '-v', action='store_false')
aparser.add_argument('--single_mode', '-s', action='store_false')
aparser.add_argument('--test', '-t', action='store_true')
aparser.add_argument('--remove_meta', '-rm', action='store_true')
aparser.add_argument('--re_convert', '-rc', action='store_true')
args = aparser.parse_args()

if args.remove_meta:
    z = DUTAnalysis(args.run, args.dut, test_campaign=args.testcampaign, single_mode=args.single_mode, verbose=False, test=True)
    z.remove_metadata()

if args.re_convert:
    z = DUTAnalysis(args.run, args.dut, test_campaign=args.testcampaign, single_mode=args.single_mode, verbose=False, test=True)
    z.remove_file()
    z.Converter.remove_raw_files()

z = DUTAnalysis(args.run, args.dut, test_campaign=args.testcampaign, single_mode=args.single_mode, verbose=args.verbose, test=args.test)

if not args.test and not z.has_alignment():
    z.Residuals.align(_save=True)
    z.REF.Residuals.align(_save=True)
    z.remove_metadata()
    z = DUTAnalysis(args.run, args.dut, test_campaign=args.testcampaign, single_mode=args.single_mode, verbose=args.verbose, test=args.test)

z.add_info(t_start, 'Init time:', prnt=True)

# aliases
try:
    d = z.Draw
    dut = z.DUT
    b = z.BeamTest
    pl = dut.Plane
    r = z.Run
    c = z.Converter
    p = z.Proteus
    cut = z.Cut
    cal = z.Calibration
    res = z.Residuals
    tel = z.Tel
    ref = z.REF
    t = z.Tracks
    cur = z.Currents
    e = z.Efficiency
    re = ref.Efficiency
    rsl = z.Resolution

except AttributeError:
    pass
