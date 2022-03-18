#!/usr/bin/env python
# --------------------------------------------------------
#       small script to start high resolution analysis
# created on March 18th 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from src.dut_analysis import DUTAnalysis, Analysis
from src.utils import *  # noqa
from plotting.draw import *
from argparse import ArgumentParser


t_start = time()
aparser = ArgumentParser()
aparser.add_argument('run', nargs='?', default=11, type=int)
aparser.add_argument('dut', nargs='?', default=1, type=int)
aparser.add_argument('--testcampaign', '-tc', nargs='?', default=Analysis.find_testcampaign())
aparser.add_argument('--verbose', '-v', action='store_true')
aparser.add_argument('--single_mode', '-s', action='store_false')
aparser.add_argument('--test', '-t', action='store_true')
args = aparser.parse_args()

z = DUTAnalysis(args.run, args.dut, test_campaign=args.testcampaign, single_mode=args.single_mode, verbose=args.verbose, test=args.test)
z.add_info(t_start, prnt=True)

# aliases
cal = z.Calibration
cut = z.Cuts
c = z.Currents
r = z.Run
dut = z.DUT
