#!/usr/bin/env python3
import sys, csv, os, argparse

parser = argparse.ArgumentParser(description='script for collapse-range',
	usage='match_counts.py countsfile psl min_read_threshold outfilename [options]')
required = parser.add_argument_group('required named arguments')
required.add_argument('counts_file', action='store',
	type=str, help='counts file from count_sam_transcripts')
required.add_argument('psl', action='store',
	type=str, help='psl or bed file to filter by counts')
parser.add_argument('min_reads', action='store',
	type=float, help='minimum number of supporting reads')
parser.add_argument('out_psl', action='store',
	type=str, help='output file name, bed or psl')
parser.add_argument('--generate_map', action='store', dest='generate_map',
	required=False, help='''name of read-isoform mapping file to be filtered''')
parser.add_argument('--append_counts', action='store_true', dest='append_counts',
	required=False, help='''specify this argument to append the supporting
	read counts as a new column in the psl file''')
args = parser.parse_args()

isbed = args.psl.lower()[-3:] != 'psl'

counts = {}
for line in open(args.counts_file):
	line = line.rstrip().split('\t')
	counts[line[0]] = float(line[1])


if args.generate_map:
	map_file = {}
	for line in open(args.generate_map):
		line = line.rstrip().split('\t')
		map_file[line[0]] = line[1]
	out_map = open(args.generate_map, 'wt')  # writing to a file of the same name
	writer_map = csv.writer(out_map, delimiter='\t', lineterminator=os.linesep)

out_psl = open(args.out_psl, 'wt')
writer = csv.writer(out_psl, delimiter='\t', lineterminator=os.linesep)
for line in open(args.psl):
	line = line.rstrip().split('\t')
	name = line[3] if isbed else line[9]

	if name in counts:
		count = counts[name]
	else:
		count = 0
	if count >= args.min_reads:
		if args.append_counts:
			writer.writerow(line + [count])
		else:
			writer.writerow(line)
		if args.generate_map:
			writer_map.writerow([name, map_file[name]])
