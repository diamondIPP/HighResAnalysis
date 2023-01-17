#!/usr/bin/env python
# --------------------------------------------------------
#       small script to start high resolution analysis
# created on March 18th 2022 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from argparse import ArgumentParser
from numpy import *

import convert
from plotting.draw import *  # noqa
import src.bins as bins  # noqa
from src.dut_analysis import DUTAnalysis, Analysis
from src.batch_analysis import BatchAnalysis
from src.run import load_nrs
from src.scan import Ensemble, Scan, VScan, TScan
from utility.utils import *  # noqa
from functools import partial

t_start = time()
aparser = ArgumentParser()
aparser.add_argument('run', nargs='?', default=Analysis.Config.get_value('data', 'default run'))
aparser.add_argument('dut', nargs='?', default=Analysis.Config.get_value('data', 'default dut', default=0), type=int)
aparser.add_argument('--batch', '-b', nargs='?', default=None, help='batch name')
aparser.add_argument('--testcampaign', '-tc', nargs='?', default=Analysis.find_testcampaign())
aparser.add_argument('--verbose', '-v', action='store_false')
aparser.add_argument('--test', '-t', action='store_true')
aparser.add_argument('--remove_meta', '-rm', action='store_true')
aparser.add_argument('--convert', '-c', action='store_true', help='removes current analysis files and reconverts from the raw files')
aparser.add_argument('--runplan', '-rp', nargs='?', default=None, help='create new runplan.json for beam test <YYYYMM>')
args = aparser.parse_args()


if args.runplan is not None:
    from src.spreadsheet import make
    make(args.runplan)
    exit(2)


ensembles = load_json(Ensemble.FilePath)

if args.run in ensembles:
    s = VScan if 'v-' in args.run else TScan if 't-' in args.run else Scan
    z = s(args.run, args.verbose, args.test)

else:

    ana = Analysis(args.testcampaign)
    runs = load_nrs(ana.BeamTest.Path)
    is_batch = not (args.run in runs and args.batch is None)
    dut_ana = partial(BatchAnalysis, choose(args.batch, args.run)) if is_batch else partial(DUTAnalysis, args.run)
    dut_ana = partial(dut_ana, args.dut, args.testcampaign)

    if is_batch:
        bc = convert.BatchConvert(dut_ana.args[0], dut_ana.args[-1], verbose=False, force=False)
        if args.convert:
            remove_file(bc.Batch.FileName)
            bc.remove_aux_files()
        if not bc.Batch.FileName.exists() and not args.test:
            bc.run()

    if args.remove_meta:
        z = dut_ana(verbose=False, test=True)
        z.remove_metadata()

    if args.convert and not is_batch:
        z = dut_ana(verbose=False, test=True)
        z.remove_file()
        z.Converter.remove_aux_files()

    z = dut_ana(verbose=args.verbose, test=args.test)

    # if not args.test and z.REF is not None and not z.has_alignment():
    #     z.Residuals.align(_save=True)
    #     z.REF.Residuals.align(_save=True)
    #     z.remove_metadata()
    #     z = DUTAnalysis(args.run, args.dut, test_campaign=args.testcampaign, single_mode=args.single_mode, verbose=args.verbose, test=args.test)

    z.add_info(t_start, 'Init time:', prnt=True)

# aliases
try:
    d = z.Draw
    dut = z.DUT
    pl = dut.Plane
    b = z.BeamTest
    r = z.Run
    c = z.Converter
    raw = c.Raw
    p = z.Proteus
    if 'CERN' in str(c):
        al = c.EventAlignment
        adc = c.Adc2Vcal
    cal = z.Calibration
    cut = z.Cut
    res = z.Residuals
    tel = z.Tel
    ref = z.REF
    t = z.Tracks
    e = z.Efficiency
    re = ref.Efficiency
    rsl = z.Resolution
    cur = z.Currents

except AttributeError:
    pass
