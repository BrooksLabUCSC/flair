#! /usr/bin/env python3

import sys, argparse, os, pipettor, glob
os.environ['OPENBLAS_NUM_THREADS'] = '1'
from bed_to_sequence import bed_to_sequence
from gtf_to_bed import gtf_to_bed
from match_counts import match_counts
import mappy as mm
from flair_align import doalignment, dofiltering
import tempfile, pysam
from flair_correct import correct
from filter_collapsed_isoforms import filter_collapsed_isoforms
from identify_gene_isoform import identify_gene_isoform
from filter_collapsed_isoforms_from_annotation import filter_collapsed_isoforms_from_annotation
from select_from_bed import select_from_bed
from pull_starts import pull_starts
from bed_to_gtf import bed_to_gtf

def getargs():
    parser = argparse.ArgumentParser(description='flair-collapse parse options',
                                     usage='''flair collapse -g genome.fa -q query.bed
    		-r <reads.fq>/<reads.fa> [options]''')
    required = parser.add_argument_group('required named arguments')
    required.add_argument('-r', '--reads', nargs='+',
                          type=str, required=True, help='FastA/FastQ files of raw reads, can specify multiple files')
    parser.add_argument('-b', '--genomealignedbam', help='Optional: sorted and indexed bam file (or files) aligned to the genome. Only use this if you have already aligned reads to the genome for some other purpose')
    parser.add_argument('-o', '--output', default='flair.collapse',
                        help='output file name base for FLAIR isoforms (default: flair.collapse)')
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='minimap2 number of threads (4)')
    parser.add_argument('-f', '--gtf', default='',
                        help='GTF annotation file, used for renaming FLAIR isoforms to annotated isoforms and adjusting TSS/TESs')
    parser.add_argument('--generate_map', default=False, action='store_true',
                        help='''specify this argument to generate a txt file of read-isoform assignments
    		(default: not specified)''')
    parser.add_argument('--transcriptfasta', default=False,
                        help='''specify transcript fasta that corresponds to transcripts in the gtf; to ask flair to make transcript sequences given the gtf and genome fa,
    		type --transcriptfasta generate''')
    # supporting read assignment options
    parser.add_argument('-s', '--support', type=float, default=3.0,
                        help='''minimum number of supporting reads for an isoform;
    		if s < 1, it will be treated as a percentage of expression of the gene (3)''')
    parser.add_argument('--stringent', default=False, action='store_true',
                        help='''specify if all supporting reads need to be full-length
    		(80%% coverage and spanning 25 bp of the first and last exons)''')
    parser.add_argument('--check_splice', default=False, action='store_true',
                        help='''enforce coverage of 4 out of 6 bp around each splice site and no
    		insertions greater than 3 bp at the splice site''')
    parser.add_argument('--trust_ends', default=False, action='store_true',
                        help='specify if reads are generated from a long read method with minimal fragmentation')
    parser.add_argument('--quality', type=int, default=1,
                        help='minimum MAPQ of read assignment to an isoform (1)')
    # variant options
    parser.add_argument('--longshot_bam', default='',
                        help='bam from longshot containing haplotype information for each read')
    parser.add_argument('--longshot_vcf', default='',
                        help='vcf from longshot')
    ###genome alignment options
    parser.add_argument('--minfragmentsize', type=int, default=80,
                        help='minimum size of alignment kept, used in minimap -s. More important when doing downstream fusion detection')
    parser.add_argument('--maxintronlen', default='200k',
                        help='maximum intron length in genomic alignment. Longer can help recover more novel isoforms with long introns')
    parser.add_argument('--filtertype', type=str, default='keepsup',
                        help='method of filtering chimeric alignments (potential fusion reads). Options: removesup (default), separate (required for downstream work with fusions), keepsup (keeps supplementary alignments for isoform detection, does not allow gene fusion detection)')
    atleastone = parser.add_argument_group('Either one of the following arguments is required')
    atleastone.add_argument('-g', '--genome', type=str,
                            help='FastA of reference genome, can be minimap2 indexed')
    atleastone.add_argument('--mm_index', type=str, default='',
                            help='minimap2 index .mmi file')
    parser.add_argument('--nvrna', action='store_true', default=False,
                        help='specify this flag to use native-RNA specific alignment parameters for minimap2')
    parser.add_argument('--junction_bed',
                        help='annotated isoforms/junctions bed file for splice site-guided minimap2 genomic alignment')
    # ends options
    parser.add_argument('-w', '--end_window', type=int, default=100,
                        help='window size for comparing TSS/TES (100)')
    parser.add_argument('-p', '--promoters', default='',
                        help='promoter regions bed file to identify full-length reads')
    parser.add_argument('--3prime_regions', dest='threeprime', default='',
                        help='TES regions bed file to identify full-length reads')
    parser.add_argument('-n', '--no_redundant', default='none',
                        help='''For each unique splice junction chain, report options include:
    		none--best TSSs/TESs chosen for each unique set of splice junctions;
    		longest--single TSS/TES chosen to maximize length;
    		best_only--single most supported TSS/TES used in conjunction chosen (none)''')
    parser.add_argument('-i', '--isoformtss', default=False, action='store_true',
                        help='''when specified, TSS/TES for each isoform will be determined from supporting reads
    		for individual isoforms (default: not specified, determined at the gene level)''')
    parser.add_argument('--no_gtf_end_adjustment', default=False, action='store_true',
                        dest='no_end_adjustment',
                        help='''when specified, TSS/TES from the gtf provided with -f will not be used to adjust
    		isoform TSSs/TESs each isoform will be determined from supporting reads''')
    parser.add_argument('--max_ends', type=int, default=2,
                        help='maximum number of TSS/TES picked per isoform (2)')
    parser.add_argument('--filter', default='default',
                        help='''Report options include:
    		nosubset--any isoforms that are a proper set of another isoform are removed;
    		default--subset isoforms are removed based on support;
    		comprehensive--default set + all subset isoforms;
    		ginormous--comprehensive set + single exon subset isoforms''')
    ##correct options
    parser.add_argument('-j', '--shortread', type=str, default='',
							help='bed format splice junctions from short-read sequencing')
    parser.add_argument('--ss_window', type=int, default=15,
                        help='window size for correcting splice sites (15)')
    # other
    parser.add_argument('--annotated_bed', default=False,
                        help='''annotation_reliant also requires a bedfile of annotated isoforms; if this isn't provided,
    		flair collapse will generate the bedfile from the gtf. eventually this argument will be removed''')
    parser.add_argument('--mm2_args', type=str, default=[],
                        help='''additional minimap2 arguments when aligning reads first-pass transcripts;
    		separate args by commas, e.g. --mm2_args=-I8g,--MD ''')
    parser.add_argument('--temp_dir', default='',
                        help='directory for temporary files. use "./" to indicate current directory (default: python tempfile directory)')
    parser.add_argument('--keep_intermediate', default=False, action='store_true',
                        help='''specify if intermediate and temporary files are to be kept for debugging.
    		Intermediate files include: promoter-supported reads file,
    		read assignments to firstpass isoforms''')

    no_arguments_passed = len(sys.argv) == 1
    if no_arguments_passed:
        parser.print_help()
        sys.exit(1)

    args, unknown = parser.parse_known_args()
    if unknown:
        sys.stderr.write('unrecognized arguments: {}\n'.format(' '.join(unknown)))

    args = checkfilepaths(args)
    args.quality = '0' if args.trust_ends else args.quality
    if args.mm2_args:
        args.mm2_args = args.mm2_args.split(',')
    return args

def maketempdir(args):
    tempfile_dir = tempfile.NamedTemporaryFile().name
    tempfile_name = tempfile_dir[tempfile_dir.rfind('/') + 1:] + '.'
    if args.temp_dir == '':
        args.temp_dir = tempfile_dir + '/'
        sys.stderr.write('Writing temporary files to {}\t\n'.format(args.temp_dir))
    if not os.path.isdir(args.temp_dir):  # make temporary directory
        pipettor.run([('mkdir', args.temp_dir)])
    if args.temp_dir[-1] != '/':
        args.temp_dir += '/'
    tempfileprefix = args.temp_dir + tempfile_name
    return tempfileprefix

def checkfilepaths(args):
    if ',' in args.reads[0]:
        args.reads = args.reads[0].split(',')
    for rfile in args.reads:
        if not os.path.exists(rfile):
            sys.stderr.write(f'Read file path does not exist: {rfile}\n')
            sys.exit(1)
    if not os.path.exists(args.genome):
        sys.stderr.write('Genome file path does not exist: {}\n'.format(args.genome))
        sys.exit(1)
    if float(args.support) < 1 and not args.gtf:
        sys.stderr.write('Provide gtf for gene grouping if -s is percentage of total gene expression\n')
        sys.exit(1)
    return args

def getannotatedseq(args):
    if args.transcriptfasta == 'generate':
        if not args.gtf:
            sys.stderr.write('Please specify annotated gtf with -f for --transcriptfasta generate\n')
            sys.exit(1)
        elif not os.path.exists(args.gtf):
            sys.stderr.write('GTF file path does not exist\n')
            sys.exit(1)

        sys.stderr.write('Making transcript fasta using annotated gtf and genome sequence\n')
        if not args.annotated_bed:
            args.annotated_bed = args.output + '.annotated_transcripts.bed'
            if not os.path.exists(args.output + '.annotated_transcripts.bed'):
                gtf_to_bed(args.annotated_bed, args.gtf, include_gene=True)
        # get transcript sequences
        args.transcriptfasta = args.output + '.annotated_transcripts.fa'
        if not os.path.exists(args.output + '.annotated_transcripts.fa'):
            bed_to_sequence(query=args.output + '.annotated_transcripts.bed', genome=args.genome,
                            outfilename=args.transcriptfasta)
    return args

def getcountsamcommand(args, outputname, mapfile, isannot):
    # count sam transcripts ; the dash at the end means STDIN
    count_cmd = ['count_sam_transcripts.py', '--sam', '-',
                 '-o', outputname, '-t', str(args.threads),
                 '--quality', str(args.quality), '-w', str(args.end_window)]
    if mapfile:
        count_cmd += ['--generate_map', mapfile]
    if args.stringent or isannot:
        count_cmd += ['--stringent']
    if args.check_splice: #or isannot:
        count_cmd += ['--check_splice']
    if args.check_splice or args.stringent or isannot:
        count_cmd += ['-i', args.annotated_bed]  # annotated isoform bed file
    if args.trust_ends:
        count_cmd += ['--trust_ends']
    return tuple(count_cmd)



def transcriptomealignandcount(args, outputname, mapfile, isannot):
    # minimap (results are piped into count_sam_transcripts.py)
    ##'--split-prefix', 'minimap2transcriptomeindex', doesn't work with MD tag
    mm2_cmd = tuple(['minimap2', '-a', '-t', str(args.threads), '-N', '4', '--MD'] + args.mm2_args + [args.transcriptfasta] + args.reads)
    ###FIXME add in step to filter out chimeric reads here
    ###FIXME really need to go in and check on how count_sam_transcripts is working
    count_cmd = getcountsamcommand(args, outputname, mapfile, isannot)
    print(' '.join(mm2_cmd))
    print(' '.join(count_cmd))
    sys.stderr.write('Aligning and counting supporting reads for transcripts\n')
    pipettor.run([mm2_cmd, count_cmd])

def subset_unassigned_reads(readmap, fastx, outfa):
    sys.stderr.write('Setting up unassigned reads for novel isoform detection\n')
    assigned_names = set()
    for line in open(readmap):  # map
        iso, reads = line.rstrip().split('\t')
        reads = reads.split(',')
        ###reads that are assigned to annotated transcripts in any counts probably won't make new transcripts, don't use support
        # if len(reads) < support:
        #   continue
        for r in reads:
            assigned_names.add(r)
    with open(outfa, 'w') as outf:
        for fle in fastx:
            for read in mm.fastx_read(fle):
                header, seq, qual = read
                if header not in assigned_names:
                    print('>' + header, file=outf)
                    print(seq, file=outf)
    outf.close()

def runreadcollapse(args):
    # TODO: collapse_isoforms_precise uses pool and map, which makes it difficult to capture in a function
    sys.stderr.write('Collapsing reads into firstpass isoforms\n')
    collapse_cmd = ['collapse_isoforms_precise.py', '-q', args.query, '-t', str(args.threads),
                    '-m', str(args.max_ends), '-w', str(args.end_window), '-n', args.no_redundant,
                    '-o', args.output + '.firstpass.unfiltered.bed']
    if args.gtf and not args.no_end_adjustment:
        collapse_cmd += ['-f', args.gtf]
    if args.isoformtss:
        collapse_cmd += ['-i']
    collapse_cmd = tuple(collapse_cmd)
    pipettor.run([collapse_cmd])

def filterlowsupisos(args):
    # filtering out subset isoforms with insufficient support
    keep_extra_column = float(args.support) < 1
    filter_collapsed_isoforms(queryfile=args.output + '.firstpass.unfiltered.bed', mode=args.filter,
                              tol=args.end_window, outfile=args.output + '.firstpass.bed',
                              keep_extra_column=keep_extra_column)

def renameisoforms(args):
    # rename first-pass isoforms to annotated transcript IDs if they match
    if args.gtf:
        sys.stderr.write('Renaming isoforms using gtf\n')
        identify_gene_isoform(query=args.output + '.firstpass.bed', gtf=args.gtf,
                              outfilename=args.output + '.firstpass.named.bed',
                              annotation_reliant=args.transcriptfasta)

        os.rename(args.output + '.firstpass.named.bed', args.output + '.firstpass.bed')

        # if we want a certain read fraction to support the isoform (instead of a number of reads)
        ###FIXME This is definitely broken, need to fix
        if float(args.support) < 1:
            filter_isoforms_by_proportion_of_gene_expr(isoforms=args.output + '.firstpass.bed',
                                                       outfilename=args.output + '.firstpass.filtered.bed',
                                                       support=args.support)
            os.rename(args.output + '.firstpass.filtered.bed', args.output + '.firstpass.bed')

def aligntofirstpasstranscripts(args):
    # get the isoform sequences for the first pass we just did
    bed_to_sequence(query=args.output + '.firstpass.bed', genome=args.genome, outfilename=args.output + '.firstpass.fa')
    sys.stderr.write('Aligning reads to first-pass isoform reference\n')
    args.transcriptfasta = args.output + '.firstpass.fa'
    args.annotated_bed = args.output + '.firstpass.bed'
    mapfile = args.output + '.novel.isoform.read.map.txt' if args.generate_map or args.transcriptfasta else None
    transcriptomealignandcount(args, args.output + '.firstpass.q.counts', mapfile, False)
    return args

def filterbasedonbed(args, filterfile, outsuffix):
    if filterfile:
        sys.stderr.write('Filtering out reads without support in ' + filterfile + '\n')
        # get the ends
        tes_bedfile = args.temp_dir + outsuffix + '.bed'
        pull_starts(args.query, tes_bedfile, reverse=True)

        # intersect with known gene ends
        bedtools_output = args.temp_dir + outsuffix + '_intersect.bed'
        bedtools_cmd = ('bedtools', 'intersect', '-a', tes_bedfile, '-b', filterfile)
        pipettor.Popen([bedtools_cmd], 'w', stdout=bedtools_output)

        # select reads that contain ends
        select_from_bed(bedtools_output, args.query, args.output + outsuffix + '_supported.bed')
        args.query = args.output + outsuffix + '_supported.bed'
        return args, [args.output + outsuffix + '_supported.bed']
    else: return args, []

def doprefiltering(args, intermediate_files):
    args.query = args.output + '.bed'
    args, newint = filterbasedonbed(args, args.promoters, 'tss')
    intermediate_files.extend(newint)
    args, newint = filterbasedonbed(args, args.threeprime, 'tes')
    intermediate_files.extend(newint)
    return args, intermediate_files

def combineannotandnovel(args, intermediate_files, min_reads, didannotreliant):
    if didannotreliant:
        filter_collapsed_isoforms_from_annotation(support=min_reads,queryfile=args.output + '.isoforms.bed',
                                                  map_i=args.output + '.novel.isoform.read.map.txt',
                                                  annotation=args.output + '.annotated_transcripts.supported.bed',
                                                  map_a=args.output + '.annotated_transcripts.isoform.read.map.txt',
                                                  outputfile=args.output + '.isoforms.filter.bed',
                                                  new_map=args.output + '.isoform.read.map.txt')
        os.rename(args.output + '.isoforms.filter.bed', args.output + '.isoforms.bed')
        intermediate_files.extend([args.output + '.novel.isoform.read.map.txt', args.output + '.annotated_transcripts.isoform.read.map.txt',
                                   args.output + '.annotated_transcripts.supported.bed'])
        if not args.generate_map: intermediate_files.extend([args.output + '.isoform.read.map.txt'])
    else: os.rename(args.output + '.novel.isoform.read.map.txt', args.output + '.isoform.read.map.txt')
    return intermediate_files

def removeintermediatefiles(args, intermediate_files):
    intermediate_files.extend([args.output + '.firstpass.unfiltered.bed', args.output + '.firstpass.bed',
                               args.output + '.firstpass.q.counts', args.output + '.firstpass.fa'])
    intermediate_files.extend(
        [args.output + '.bed', args.output + '_all_corrected.bed', args.output + '_all_inconsistent.bed'])
    if not args.keep_intermediate:
        intermediate_files.extend(glob.glob(args.temp_dir + '*'))  # TODO: CHECK
        pipettor.run([tuple(['rm'] + intermediate_files)])

def matchcountsonfirstpass(args, min_reads):
    if args.generate_map or args.transcriptfasta:
        isoform_file = args.output + '.novel.isoform.read.map.txt'
    else:
        isoform_file = False
    match_counts(counts_file=args.output + '.firstpass.q.counts', output_file=args.output + '.isoforms.bed',
                 bed=args.output + '.firstpass.bed', min_reads=min_reads, isoform_file=isoform_file)

def collapse(args):
    args.temp_dir = maketempdir(args)
    ###FIXME this invalidates goal of allowing fractional support - tbh maybe fractional support should be a separate option?
    min_reads = float(args.support) if float(args.support) >= 1 else 3
    intermediate_files = []
    #check whether any annotation
    didannotreliant = True if args.transcriptfasta else False
    if args.transcriptfasta:
        args = getannotatedseq(args)
        transcriptomealignandcount(args, args.output + '.annotated_transcripts.alignment.counts', args.output + '.annotated_transcripts.isoform.read.map.txt', True)
        ###outputs only supported annotated transcripts to a new file
        ###FIXME I feel like it's unnecessary to have both the counts file and the read map file - remove the counts file
        match_counts(counts_file=args.output + '.annotated_transcripts.alignment.counts',
                     output_file=args.output + '.annotated_transcripts.supported.bed',
                     bed=args.annotated_bed, min_reads=args.support)
        intermediate_files.append(args.output + '.annotated_transcripts.alignment.counts')

    if os.stat(args.reads[0]).st_size > 0: ###check that there are remaining reads to align to the genome
        if args.genomealignedbam:
            if args.transcriptfasta:
                dofiltering(args, args.genomealignedbam, args.output + '.annotated_transcripts.isoform.read.map.txt')
            else:
                dofiltering(args, args.genomealignedbam)
        else:
            if args.transcriptfasta:
                # get the unassigned reads separately
                subset_reads = args.output + '.unassigned.fasta'
                subset_unassigned_reads(readmap=args.output + '.annotated_transcripts.isoform.read.map.txt',
                                        fastx=args.reads, outfa=subset_reads)
                args.reads = [subset_reads]
                intermediate_files.append(args.output + '.unassigned.fasta')
            ##do alignment and filtering
            sys.stderr.write('Aligning to genome\n')
            doalignment(args)
            dofiltering(args, args.output + '_unfiltered.bam')
            pipettor.run([('rm', args.output + '_unfiltered.bam', args.output + '_unfiltered.bam.bai')])

        args, intermediate_files = doprefiltering(args, intermediate_files)
        if args.shortread or args.gtf: ##only do correction if have gtf or shortread
            sys.stderr.write('Correcting splice sites\n')
            correct(args=args)
            args.query = args.output + '_all_corrected.bed'
        ###FIXME we must load the gtf file in every step here, is there a way to just keep it loaded into memory and pass it to different steps?
        runreadcollapse(args)
        filterlowsupisos(args)
        renameisoforms(args)
        args = aligntofirstpasstranscripts(args)
        matchcountsonfirstpass(args, min_reads)

    intermediate_files = combineannotandnovel(args, intermediate_files, min_reads, didannotreliant)

    sys.stderr.write('Generating final transcriptome fasta and gtf\n')
    bed_to_sequence(query=args.output + '.isoforms.bed', genome=args.genome, outfilename=args.output + '.isoforms.fa')
    if args.gtf:
        bed_to_gtf(query=args.output + '.isoforms.bed', outputfile=args.output + '.isoforms.gtf')
    removeintermediatefiles(args, intermediate_files)

if __name__ == "__main__":
    args = getargs()
    collapse(args)
