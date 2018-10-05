#!/usr/bin/env python
# --------------------------------------------------------
#       small script to select runs by batch
# created on October 4th 2018 by M. Reichmann (remichae@phys.ethz.ch)
# --------------------------------------------------------

from json import load, loads
from os.path import *
from argparse import ArgumentParser

p = ArgumentParser()
p.add_argument('batches')
args = p.parse_args()

batches = [args.batches] if not ',' in args.batches else [b.strip(' ') for b in args.batches.split(',')]
print batches

d = dirname(dirname(realpath(__file__)))

f = open(join(d, 'runlog.json'))
info = load(f)
f.close()

nrs = []
for key, dic in sorted(info.iteritems()):
    if dic['Batch'] in batches:
        if dic['Status'] == 'green':
            nrs.append(dic['Run Number']['CMS'].zfill(3))

for nr in sorted(nrs):
    print 'ljutel_{}.root'.format(nr),