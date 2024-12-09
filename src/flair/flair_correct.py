#!/usr/bin/env python3

import argparse
import sys
from multiprocessing import Pool
from tqdm import tqdm
import os
import pybedtools
import shutil
import uuid
from ssUtils import addOtherJuncs, gtfToSSBed
from ssPrep import ssPrep

def parseargs(aligned_reads=''):
	parser = argparse.ArgumentParser(description='flair-correct parse options',
									 usage='flair correct -q query.bed12 [-f annotation.gtf]v[-j introns.tab] -g genome.fa [options]')
	required = parser.add_argument_group('required named arguments')
	atleastone = parser.add_argument_group('at least one of the following arguments is required')
	if not aligned_reads:
		required.add_argument('-q', '--query', type=str, required=True,
							  help='uncorrected bed12 file')
	required.add_argument('-g', '--genome', type=str, required=True,
						  help='FastA of reference genome')
	atleastone.add_argument('-j', '--shortread', type=str, default='',
							help='bed format splice junctions from short-read sequencing')
	atleastone.add_argument('-f', '--gtf', default='',
							help='GTF annotation file')
	parser.add_argument('-o', '--output', default='flair',
						help='output name base (default: flair)')
	parser.add_argument('-t', '--threads', type=int, default=4,
						help='number of threads (4)')
	parser.add_argument('--nvrna', action='store_true', default=False,
						help='''specify this flag to make the strand of a read consistent with the annotation during correction''')
	parser.add_argument('-w', '--ss_window', type=int, default=15,
						help='window size for correcting splice sites (15)')
	# parser.add_argument('--print_check', action='store_true', default=False,
	# 					help='Print err.txt with step checking.')
	no_arguments_passed = len(sys.argv) == 1
	if no_arguments_passed:
		parser.print_help()
		sys.exit(1)

	args, unknown = parser.parse_known_args()
	if unknown:
		sys.stderr.write('Correct unrecognized arguments: {}\n'.format(' '.join(unknown)))
	return args

def correct(aligned_reads='', args=None):
	if not args: args = parseargs(aligned_reads)

	if aligned_reads:
		query = aligned_reads
	else:
		query = args.query

# TODO:This seems opposite the intended use, see what happens
	resolveStrand = False
	if not args.nvrna:
		resolveStrand = True

	if os.path.isfile(args.genome+".fai"):
		pass
	else:
		testString = """
			chrX 1	100   feature1  0 +
		"""
		# TODO: Does this create .fai?
		pybedtools.BedTool(testString, from_string=True)

	# make temp dir for dumping
	tempDirName = str(uuid.uuid4())
	try:
		current_directory = os.getcwd()
		tempDir = os.path.join(current_directory, tempDirName)
		os.mkdir(tempDir)
	except OSError:
		print("Creation of the directory %s failed" % tempDirName)
		sys.exit(1)

	# There are a few functions that evaluate what verbose is defined as.
	# Instead of passing it around, just global it.
	global verbose
	global printErrFname
	global printErr
	verbose  = False # TODO
	# printErr = args.print_check
	printErr = False
	printErrFname = False
	if printErr:
		printErrFname = os.path.join(tempDirName, 'ssCorrect.err')

	# Convert gtf to bed and split by chromosome.
	juncs, chromosomes, knownSS  = dict(), set(), dict() # initialize juncs for adding to db
	if args.gtf: 
		juncs, chromosomes, knownSS = gtfToSSBed(args.gtf, knownSS, printErr, printErrFname, verbose)

	# Do the same for the other juncs file.
	if args.shortread: 
		juncs, chromosomes, addFlag = addOtherJuncs(juncs, args.shortread, chromosomes, args.genome, 
			printErrFname, knownSS, verbose, printErr)
		if addFlag == False:
			sys.stderr.write('\nERROR Added no extra junctions from {}\n\n'.format(args.shortread))  
			sys.exit(1)
	knownSS = dict()

	# added to allow annotations not to be used.
	if len(list(juncs.keys())) < 1:
		print("No junctions from GTF or junctionsBed to correct with. Exiting...", file=sys.stderr)
		sys.exit(1)

	annotations = dict()
	for chrom, data in tqdm(juncs.items(), desc="Step 3/5: Preparing annotated junctions to use for correction", total=len(list(juncs.keys())), dynamic_ncols=True, position=1) if verbose else juncs.items():
		annotations[chrom] = os.path.join(tempDir,"%s_known_juncs.bed" % chrom)
		with open(os.path.join(tempDir,"%s_known_juncs.bed" % chrom),"w") as out:
			sortedData = sorted(list(data.keys()), key=lambda item: item[0])
			for k in sortedData:
				annotation = data[k]
				c1, c2, strand = k
				print(chrom,c1,c2,annotation,".",strand, sep="\t", file=out)

	sortedData = None
	skippedChroms = set()
	readDict = dict()
	notfound = False
	prevchrom = False
	with open(query) as lines, open("%s_cannot_verify.bed" % args.output,'w') as nochrom:
		outDict = dict()
		for line in tqdm(lines, desc="Step 4/5: Preparing reads for correction", dynamic_ncols=True, position=1) if verbose else lines:
			cols  = line.rstrip().split()
			chrom = cols[0]
			if chrom not in chromosomes:
				notfound = True
				nochrom.write(line)
				if chrom not in skippedChroms:
					skippedChroms.add(chrom)
					if verbose: tqdm.write("Reference sequence %s not found in annotations, skipping" % (chrom), file=sys.stdout)
					continue
			else:
				if chrom not in outDict:
					readDict[chrom] = os.path.join(tempDir,"%s_temp_reads.bed" % chrom)
					outDict[chrom] = os.path.join(tempDir,"%s_temp_reads.bed" % chrom)
					#outDict[chrom] = open(os.path.join(tempDir,"%s_temp_reads.bed" % chrom),'w')
				if chrom != prevchrom:
					if prevchrom is not False:
						tempoutfile.close()
					tempoutfile = open(outDict[chrom], 'a')
					prevchrom = chrom
				print(line.rstrip(),file=tempoutfile)
				#with open(outDict[chrom], 'a') as tempoutfile:
				#print(line.rstrip(),file=outDict[chrom])
	nochrom.close()
	tempoutfile.close()
	if notfound is False:
		os.remove(f'{args.output}_cannot_verify.bed')

	cmds = list()
	for chrom in readDict:
		juncs = annotations[chrom]
		reads = readDict[chrom]

#		outDict[chrom].close()

		cmds.append([reads, juncs, args.genome, args.ss_window, chrom, resolveStrand,
	      tempDir, printErrFname])

	if printErr:
		with open(printErrFname,'a+') as fo:
			print("** Prepared correct commands for %s read files" % len(cmds), file=fo)

	juncs = None
	annotations = None
	p = Pool(args.threads)
	childErrs = set()
	for i in tqdm(p.imap(ssPrep, cmds), total=len(cmds), desc="Step 5/5: Correcting Splice Sites", 
	       dynamic_ncols=True,position=1) if verbose else p.imap(ssPrep,cmds):
		childErrs.add(i)
	p.close()
	p.join()
	if len(childErrs) > 1:
		print(childErrs,file=sys.stderr)
		sys.exit(1)

	with open("%s_all_inconsistent.bed" % args.output,'wb') as inconsistent:
		for chrom in readDict:
			with open(os.path.join(tempDir, "%s_inconsistent.bed" % chrom),'rb') as fd:
				shutil.copyfileobj(fd, inconsistent, 1024*1024*10)

	correct_bed = args.output + '_all_corrected.bed' 
	with open(correct_bed,'wb') as corrected:
		for chrom in readDict:
			with open(os.path.join(tempDir, "%s_corrected.bed" % chrom),'rb') as fd:
				shutil.copyfileobj(fd, corrected, 1024*1024*10)
	if printErr: 
		shutil.move(printErrFname, f'{args.output}.err') 
	shutil.rmtree(tempDir)

	return correct_bed

if __name__ == '__main__':
	correct()
