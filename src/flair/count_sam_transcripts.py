#!/usr/bin/env python3

import sys
import argparse
import re
import csv
import math
import os
from collections import Counter
from collections import namedtuple


def parseargs():
	parser = argparse.ArgumentParser(description='''for counting transcript abundances after
		aligning reads to transcripts; for multiple mappers, only the best alignment
		for each read is used, usage=python -s samfile -o outputfile''')
	required = parser.add_argument_group('required named arguments')
	required.add_argument('-s', '--sam', type=argparse.FileType('r'), help='sam file or - for STDIN')
	required.add_argument('-o', '--output', default='counts.txt', help='output file name')
	parser.add_argument('-i', '--isoforms',
						help='specify isoforms.bed file if --stringent and/or --check_splice is specified')
	parser.add_argument('--stringent', action='store_true',
						help='only count if read alignment passes stringent criteria')
	parser.add_argument('--check_splice', action='store_true',
						help='''enforce coverage of 4 out of 6 bp around each splice site and no
		insertions greater than 3 bp at the splice site''')
	parser.add_argument('--trust_ends', action='store_true',
						help='specify if reads are generated from a long read method with minimal fragmentation')
	parser.add_argument('-t', '--threads', default=4, type=int,
						help='number of threads to use')
	parser.add_argument('-w', '--window', type=int, default=10,
						help='number of bases for determining which end is best match (10)')
	parser.add_argument('--quality', default=1, type=int,
						help='minimum quality threshold to consider if ends are to be trusted (1)')
	parser.add_argument('--generate_map',
						help='''specify an output path for a txt file of which isoform each read is assigned to''')
	parser.add_argument('--fusion_dist',
						help='''minimium distance between separate read alignments on the same chromosome to be
		considered a fusion, otherwise no reads will be assumed to be fusions''')
	# parser.add_argument('--minimal_input',
	# 	help='''input file is not actually a sam file, but only the info necessary''')
	args = parser.parse_args()
	return args


def checkargs(args):
	if args.stringent or args.fusion_dist or args.check_splice:
		if not os.path.exists(args.isoforms):
			sys.stderr.write('A valid isoforms bed file needs to be specified: {}\n'.format(args.isoforms))
			sys.exit(1)
	if args.fusion_dist:
		args.trust_ends = True
	return args


def getannotinfo(args):
	transcripttoexons = {}
	if args.stringent or args.check_splice or args.fusion_dist:
		for line in open(args.isoforms):
			line = line.rstrip().split('\t')
			name, left, right, chrom = line[3], int(line[1]), int(line[2]), line[0]
			blocksizes = [int(n) for n in line[10].rstrip(',').split(',')]

			if line[5] == '+': transcripttoexons[name] = blocksizes
			else: transcripttoexons[name] = blocksizes[::-1]

	return transcripttoexons

def check_singleexon(read_start, read_end, tlen):
	if read_start < 25 and read_end > tlen - 25:
		return True
	else:
		return False

def check_exonenddist(blocksize, disttoend, trust_ends, disttoblock):
	if blocksize < 25: # if the terminal exon is sufficiently short, relax by 5 bp bc drna misses bases on 5' end
		return disttoend < 5
	elif trust_ends:
		return disttoend <= trust_ends_window
	else:
		return disttoblock > 25

def check_firstlastexon(first_blocksize, last_blocksize, read_start, read_end, tlen, trust_ends):
	left_coverage = check_exonenddist(first_blocksize, read_start, trust_ends, first_blocksize-read_start)
	right_coverage = check_exonenddist(last_blocksize, tlen-read_end, trust_ends, read_end - (tlen - last_blocksize))
	return right_coverage and left_coverage

def check_stringent(coveredpos, exonpos, tlen, blockstarts, blocksizes, trust_ends):
	matchpos = len([x for x in coveredpos if x == 1])
	if matchpos/tlen < 0.8:
		return False
	else:
		read_start, read_end = blockstarts[0], blockstarts[-1] + blocksizes[-1]
		first_blocksize, last_blocksize = exonpos[0], exonpos[-1]
		# covers enough bases into the first and last exons
		if len(blocksizes) == 1:  # single exon transcript
			return check_singleexon(read_start, read_end, tlen)
		else:
			return check_firstlastexon(first_blocksize, last_blocksize, read_start, read_end, tlen, trust_ends)

def check_splicesites(coveredpos, exonpos, tstart, tend):
	currpos = 0
	for i in range(len(exonpos)-1):
		elen = exonpos[i]
		currpos += elen
		if tstart < currpos < tend:
			ssvals = coveredpos[currpos - 3:currpos + 3]
			# print(currpos, ssvals)
			totinsert = sum([x for x in ssvals if x > 1])
			totmatch = sum([x for x in ssvals if x == 1])
			if totmatch - totinsert <= num_match_in_ss_window:
				return False
	return True

def get_matchvals(args, md):
	matchvals = []
	if args.stringent or args.check_splice:
		mdblocks = re.findall('\d+|\D+', md)  # ['531', '^CCAGGTGAGCCGCCCGCG', '50', 'G', '2031']
		for b in mdblocks:
			if b[0] != '^':
				if b.isnumeric():
					matchvals.extend([1] * int(b))
				else:
					matchvals.append(0)
	return matchvals

def process_cigar(args, matchvals, cigar, startpos):
	matchpos = 0
	cigarblocks = re.findall('([0-9]+)([A-Z])', cigar)
	coveredpos = [0] * (startpos - 1)
	queryclipping = []
	tendpos = startpos - 1
	# indel_detected = False
	blockstarts, blocksizes = [], []
	for blen, btype in cigarblocks:
		blen = int(blen)
		if btype == 'S' or btype == 'H':
			queryclipping.append(blen)
		elif btype == 'M' and (args.stringent or args.check_splice):
			# coveredpos.extend([1] * blen)
			coveredpos.extend(matchvals[matchpos:matchpos + blen])
			blockstarts.append(tendpos)
			blocksizes.append(blen)
			matchpos += blen
			tendpos += blen
		elif btype == 'D' or btype == 'N':
			if args.stringent or args.check_splice:
				coveredpos.extend([0] * blen)
				tendpos += blen
			if blen > large_indel_tolerance: return True, None, None, None, None, None #
		elif btype == 'I':
			if args.stringent or args.check_splice:
				coveredpos[-1] += blen
			if blen > large_indel_tolerance: return True, None, None, None, None, None
	return False, coveredpos, queryclipping, blockstarts, blocksizes, tendpos

def checktranscriptinannot(exondict, tname):
	try:
		exoninfo = exondict[tname]
	except KeyError:
		raise Exception(
			"The transcript names in the annotation fasta do not appear to match the ones in the isoforms file. You may be able to fix this by using gtf_to_bed and bed_to_sequence on your annotation gtf and using the resulting file as your annotation fasta input to this program")
	except Exception as ex:
		raise Exception("** check_splice FAILED for %s" % (tname)) from ex
	return exoninfo

def check_stringentandsplice(args, transcripttoexons, tname, coveredpos, tlen, blockstarts, blocksizes, tstart, tend):
	passesstringent, passessplice = True, True
	if args.stringent or args.check_splice:
		exoninfo = checktranscriptinannot(transcripttoexons, tname)
		passesstringent = check_stringent(coveredpos, exoninfo, tlen, blockstarts, blocksizes,
										  args.trust_ends) if args.stringent else True
		passessplice = check_splicesites(coveredpos, exoninfo, tstart, tend) if args.check_splice else True
	return passesstringent and passessplice



def getbesttranscript(tinfo, args, transcripttoexons):
	###parse CIGAR + MD tag to ID transcript pos covered by alignment
	###get start + end of transcript on read, alignment block positions
	###also save soft/hard clipping at ends of read
	###not positions of insertions larger than min_insertion_len, apply those to check_splice
	##filter out reads with long indels
	###generate list of 0s and 1s - transcript pos with match to query, val > 1 = insertion
	passingtranscripts = []
	for tname in tinfo:
		thist = tinfo[tname]
		###process MD tag here to query positions with mismatches
		##for MD tag, keep track of position of mismatch in all match positions
		matchvals = get_matchvals(args, thist.md)

		indel_detected, coveredpos, queryclipping, blockstarts, blocksizes, tendpos = process_cigar(args, matchvals, thist.cigar, thist.startpos)

		if not indel_detected:
			if check_stringentandsplice(args, transcripttoexons, thist.name, coveredpos, thist.tlen, blockstarts, blocksizes, thist.startpos, tendpos):
				passingtranscripts.append([-1 * thist.alignscore, sum(queryclipping), thist.tlen, tname])
	###order passing transcripts by alignment score
	###then order by amount of query covered
	##then order by amount of transcript covered
	if len(passingtranscripts) > 0:
		passingtranscripts.sort()
		return passingtranscripts[0][-1]
	else: return None

class IsoAln(object):
	def __init__(self, name=None, q=None, p=None, cigar=None, tlen=None, als=None, md=None):
		self.name = name
		self.q = q
		self.startpos = p
		self.cigar = cigar
		self.tlen = tlen
		self.alignscore = als
		self.md = md

def parsesamheaderline(line):
	name = line[line.find('SN:') + 3:].split()[0]
	if name.count('|') > 2:
		name = name[:name.find('|')]
	length = float(line[line.find('LN:') + 3:].split()[0])
	return name, length

def gettranscriptinfo(line):
	line = line.rstrip().split('\t')
	read, transcript, cigar, quality, pos = line[0], line[2], line[5], int(line[4]), int(line[3]) - 1
	if transcript.count('|') > 2: transcript = transcript[:transcript.find('|')]
	tags = ' '.join(line[11:])
	if transcript == '*': alignscore, mdtag = None, None
	else:
		alignscore = int(tags.split('AS:i:')[1].split(' ')[0])
		mdtag = tags.split('MD:Z:')[1].split(' ')[0]
	return read, transcript, cigar, quality, pos, alignscore, mdtag


def parsesam(args, transcripttoexons):
	transcript_lengths = {}
	lastread = None
	curr_transcripts = {}
	transcripttoreads = {}
	for line in args.sam:
		if line.startswith('@SQ'):
			name, length = parsesamheaderline(line)
			transcript_lengths[name] = length
		elif not line[0] == '@':
			read, transcript, cigar, quality, pos, alignscore, mdtag = gettranscriptinfo(line)
			if quality >= args.quality and transcript != '*':
				if lastread and read != lastread:
					assignedt = getbesttranscript(curr_transcripts, args, transcripttoexons)
					if assignedt:
						if assignedt not in transcripttoreads: transcripttoreads[assignedt] = []
						transcripttoreads[assignedt].append(lastread)
					curr_transcripts = {}
				tlen = transcript_lengths[transcript]
				curr_transcripts[transcript] = IsoAln(transcript, quality, pos, cigar, tlen, alignscore, mdtag)
				lastread = read
	if lastread:
		assignedt = getbesttranscript(curr_transcripts, args, transcripttoexons)
		if assignedt:
			if assignedt not in transcripttoreads: transcripttoreads[assignedt] = []
			transcripttoreads[assignedt].append(lastread)
	return transcripttoreads

def write_output(args, transcripttoreads):
	if args.generate_map: mapout = open(args.generate_map, 'w')
	countout = open(args.output, 'wt')
	for t in transcripttoreads:
		if args.generate_map:
			mapout.write(t + '\t' + ','.join(transcripttoreads[t]) + '\n')
		countout.write(t + '\t' + str(len(transcripttoreads[t])) + '\n')


min_insertion_len = 3
# ss_window = 3 unused
num_match_in_ss_window = 4
trust_ends_window = 50
large_indel_tolerance = 25

if __name__ == '__main__':
	args = parseargs()
	args = checkargs(args)
	transcripttoexons = getannotinfo(args)
	transcripttoreads = parsesam(args, transcripttoexons)
	write_output(args, transcripttoreads)
