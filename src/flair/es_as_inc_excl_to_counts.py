#!/usr/bin/env python3
from __future__ import print_function
import os, sys
import numpy as np

data = dict()
with open(sys.argv[1]) as fin1:
    header = next(fin1).rstrip().split()
    nSamps = len(header) - 1
    for line in fin1:
        cols = line.rstrip().split()
        tid, vals = cols[0] , np.asarray(cols[1:], dtype=np.float32)
        data[tid] = vals

with open(sys.argv[2]) as fin2:
    print('\t'.join(['feature_id', 'coordinate']+header[1:]+['isoform_ids']))
    for line in fin2:
        cols = line.rstrip().split()
        if int(cols[3]) == 0:
            continue
        else:
            exon, strand, inc, exc, incIsos, excIsos = cols
        incVals = np.asarray([data.get(x,np.zeros(nSamps)) for x in incIsos.split(",")])
        excVals = np.asarray([data.get(x,np.zeros(nSamps)) for x in excIsos.split(",")])

        incVals = np.sum(incVals, axis=0)
        excVals = np.sum(excVals, axis=0)
        totVals = incVals + excVals

        # must have at least 1 count to support the inc and exc of this exon
        # if incVals.all() < 1 or excVals.all() < 1:
        #     continue

        print("inclusion_%s" % exon,exon,"\t".join(str(x) for x in incVals),incIsos,sep="\t")
        print("exclusion_%s" % exon,exon,"\t".join(str(x) for x in excVals),excIsos,sep="\t")

        #print(exon,"\t".join("%.2f" % x for x in incVals/totVals))

