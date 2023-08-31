#!/usr/bin/env python3

import sys
from tqdm import tqdm
import pybedtools
import re

def addOtherJuncs(juncs, bedJuncs, chromosomes, fa, printErrFname, known, verbose, printErr):
	verbose = False
	if verbose: sys.stderr.write("Step 2/5: Processing additional junction file  %s ..." % (bedJuncs))
	cols = None

	with open(bedJuncs) as l:
		for num,ll in enumerate(l,0):
			cols = ll.rstrip().split()
			if num > 10:
				break

	# guess what kind of bedFile
	if cols is None:
		raise Exception("Empty junctions BED file, not supported")

	if cols[-1] == "+" or cols[-1] == "-":
		# normal bed
		strandCol = -1
		starOffset = 0

	elif len(cols) == 12:
		# bed12
		raise Exception("Bed12 not currently supported for other_juncs.bed. Please convert to bed6.")

	elif cols[3] == "0" or cols[3] == "1" or cols[3] == "2":
		# star junc.tab
		strandCol = 3
		starOffset = 1

	else:
		raise Exception("Cannot find strand info for %s. Is this bed6 or STAR_juncs.tab file?" % bedJuncs)

	if printErr:
		with open(printErrFname,'a+') as fo:
			print("** Adding other juncs, assuming file is %s" % "bed6" if strandCol == -1 else "STAR", file=fo)

	tempJuncs = list()
	addedFlag = False
	with open(bedJuncs,'r') as bedLines:
		for line in bedLines:
			cols = line.rstrip().split()
			chrom, c1, c2, strand = cols[0], int(cols[1])-starOffset, int(cols[2]), cols[strandCol]

			if chrom not in juncs:
				juncs[chrom] = dict()

			if c2-c1 < 5:
				continue

			if starOffset:
				if strand == "1": strand = "+"
				elif strand == "2": strand = "-"
				else: continue

			chromosomes.add(chrom)
			key = (c1, c2, strand)
			if key in juncs[chrom]:
				juncs[chrom][key] = "both"
				continue
			tempJuncs.append((chrom,c1,c2,"%s,%s,%s,%s" % (chrom,c1,c2,strand),0,strand))
			addedFlag = True
	if addedFlag == False:
		return juncs, chromosomes, addedFlag

	try:
		btJuncs = pybedtools.BedTool(tempJuncs)
		dinucSeq = btJuncs.sequence(fi=fa, s=True, tab=True, name=True)
		with open(dinucSeq.seqfn) as fileObj:
			for i in fileObj:
				header,seq = i.rstrip().split()
				chrom,c1,c2,strand = header.split(",")
				c1,c2 = int(c1),int(c2)
				if "+" in strand:
					strand = strand[strand.rfind('+')]
				elif "-" in strand:
					strand = strand[strand.rfind('-')]
				key = (c1,c2, strand)
				known1,known2 = known.get((chrom,c1),None),known.get((chrom,c2),None)

				if known1 is not None:
					if known1 != strand:
						continue
					else:
						pass
				else:
					pass

				if known2 is not None:
					if known2 != strand:
						continue
					else:
						pass
				else:
					pass

				fivePrime = seq[:2]
				if key not in juncs[chrom]:
					juncs[chrom][key] = "sr"
				elif fivePrime == "GT":
					juncs[chrom][key] = "sr"

	except Exception as e:
		print(e,"Splice site motif filtering failed. Check pybedtools and bedtools is properly install and in $PATH",file=sys.stderr)
		sys.exit(1)

	if printErr:
		with open(printErrFname,'a+') as fo:
			print("** GTF Juncs + other juncs now total %s juncs from %s chromosomes." % (sum([len(x)for x in juncs.values()]), len(list(juncs.keys()))), file=fo)

	return juncs, chromosomes, addedFlag


def gtfToSSBed(file, knownSS, printErr, printErrFname, verbose):
	''' Convenience function, reformats GTF to bed'''

	# First: get all exons per transcript.
	exons = dict()
	chromosomes = set()
	with open(file,'r') as lines:
		for l in lines:
			if l[0] == "#": # skip header lines
				continue

			cols = l.split("\t")

			if "exon" == cols[2]:

				# -1 for 1 to 0 based conversion
				chrom, c1, c2, strand = cols[0], int(cols[3])-1, int(cols[4]), cols[6]
				chromosomes.add(chrom)
				#txn info is in the SECOND position of the shoutout column
				txn = re.search('transcript_id "([^\"]+)"', l).group(1)#
				#cols[-1].split(";")[1].split()[-1].replace('"','')

				key = (chrom, txn, strand)

				if key not in exons:
					exons[key] = list()
				exons[key].append(c1)
				exons[key].append(c2)

	if printErr:
		with open(printErrFname,'a+') as fo:
			print("** Read GTF. Got %s transcripts" % len(list(exons.keys())), file=fo)
			print("** Getting introns...Read GTF", file=fo)

	# Second: get junction and splice sites from transcript exons.
	txnList = list(exons.keys())
	juncs = dict()

	for exonInfo in tqdm(txnList, total=len(txnList), desc="Step 1/5: Splitting junctions from GTF by chromosome", dynamic_ncols=True, position=1) if verbose else txnList:
		chrom, txn, strand = exonInfo

		if chrom not in juncs:
			juncs[chrom] = dict()

		coords = list(exons[exonInfo])

		# assume lowest and highest as TSS and TES, and remove them
		coords.sort()
		coords = coords[1:-1]

		# Coords is list of exons, so a list less than 2 is a single exon gene.
		if len(coords) < 2: continue

		for pos in range(0,len(coords)-1,2):
			c1 = coords[pos]
			c2 = coords[pos+1]

			if abs(c2 - c1) <= 5:
				continue

			juncs[chrom][(c1,c2,strand)] = "gtf"
			knownSS[(chrom, c1)] = strand
			knownSS[(chrom, c2)] = strand

	if printErr:
		with open(printErrFname,'a+') as fo:
			print("** Created %s juncs from %s chromosomes." % (sum([len(x)for x in juncs.values()]), len(list(juncs.keys()))), file=fo)
	return juncs, chromosomes, knownSS





