#!/usr/bin/env python3

import sys, csv, os
import scipy.stats as sps

try:
	events_quant = open(sys.argv[1])
	colname1 = sys.argv[2]
	colname2 = sys.argv[3]
	outfilename = sys.argv[4]
except:
	print('usage: script.py events.quant.tsv colname1 colname2 out.fishers.tsv')
	sys.exit()

header = events_quant.readline().rstrip().split('\t')

if colname1 in header:
	col1 = header.index(colname1)
else:
	sys.stderr.write('Could not find {} in {}\n'.format(colname1, ' '.join(header)))
	sys.exit(1)

if colname2 in header:
	col2 = header.index(colname2)
else:
	sys.stderr.write('Could not find {} in {}\n'.format(colname2, ' '.join(header)))
	sys.exit(1)

events = {}
for line in events_quant:
	line = line.rstrip().split('\t')
	feature = line[0][line[0].find('_')+1:]
	if feature not in events:
		events[feature] = {}
		events[feature]['entries'] = []
		events[feature]['counts'] = []
	events[feature]['entries'] += [line]
	events[feature]['counts'] += [[float(line[col1]), float(line[col2])]]

features_sorted = sorted(events.keys())
with open(outfilename, 'wt') as outfile:
	writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
	writer.writerow(header+[colname1+'-'+colname2+'_pval'])
	for feature in features_sorted:
		for line in events[feature]['entries']:
			writer.writerow(line + [sps.fisher_exact(events[feature]['counts'])[1]])
