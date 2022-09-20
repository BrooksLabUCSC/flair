flair collapse
=======

usage: python flair.py collapse -g genome.fa -q <query.psl>|<query.bed>
		-r <reads.fq>/<reads.fa> [options]

flair-collapse parse options

positional arguments:
  collapse

optional arguments:
  -h, --help            show this help message and exit
  -f F, --gtf F         GTF annotation file, used for renaming FLAIR isoforms
                        to annotated isoforms and adjusting TSS/TESs
  -t T, --threads T     minimap2 number of threads (4)
  -o O, --output O      output file name base for FLAIR isoforms (default:
                        flair.collapse)
  --generate_map        specify this argument to generate a txt file of read-
                        isoform assignments note: only works if the
                        quantification method is not using salmon (default:
                        not specified)
  --annotation_reliant ANNOTATION_RELIANT
                        specify transcript fasta that corresponds to
                        transcripts in the gtf to run annotation- reliant
                        flair collapse; to ask flair to make transcript
                        sequences given the gtf and genome fa, type
                        --annotation_reliant generate
  -s S, --support S     minimum number of supporting reads for an isoform; if
                        s < 1, it will be treated as a percentage of
                        expression of the gene (3)
  --stringent           specify if all supporting reads need to be full-length
                        (80% coverage and spanning 25 bp of the first and last
                        exons)
  --check_splice        enforce coverage of 4 out of 6 bp around each splice
                        site and no insertions greater than 3 bp at the splice
                        site
  --trust_ends          specify if reads are generated from a long read method
                        with minimal fragmentation
  --quality QUALITY     minimum MAPQ of read assignment to an isoform (1)
  --longshot_bam LONGSHOT_BAM
                        bam from longshot containing haplotype information for
                        each read
  --longshot_vcf LONGSHOT_VCF
                        vcf from longshot
  -w W, --end_window W  window size for comparing TSS/TES (100)
  -p P, --promoters P   promoter regions bed file to identify full-length
                        reads
  --3prime_regions THREEPRIME
                        TES regions bed file to identify full-length reads
  -n N, --no_redundant N
                        For each unique splice junction chain, report options
                        include: none--best TSSs/TESs chosen for each unique
                        set of splice junctions; longest--single TSS/TES
                        chosen to maximize length; best_only--single most
                        supported TSS/TES used in conjunction chosen (none)
  -i, --isoformtss      when specified, TSS/TES for each isoform will be
                        determined from supporting reads for individual
                        isoforms (default: not specified, determined at the
                        gene level)
  --no_gtf_end_adjustment
                        when specified, TSS/TES from the gtf provided with -f
                        will not be used to adjust isoform TSSs/TESs each
                        isoform will be determined from supporting reads
  --max_ends MAX_ENDS   maximum number of TSS/TES picked per isoform (2)
  --filter FILTER       Report options include: nosubset--any isoforms that
                        are a proper set of another isoform are removed;
                        default--subset isoforms are removed based on support;
                        comprehensive--default set + all subset isoforms;
                        ginormous--comprehensive set + single exon subset
                        isoforms
  --temp_dir TEMP_DIR   directory for temporary files. use "./" to indicate
                        current directory (default: python tempfile directory)
  --keep_intermediate   specify if intermediate and temporary files are to be
                        kept for debugging. Intermediate files include:
                        promoter-supported reads file, read assignments to
                        firstpass isoforms
  -m M, --minimap2 M    path to minimap2 if not in $PATH
  -b B, --bedtools B    bedtools executable path, provide if TSS/TES regions
                        specified and bedtools is not in $PATH
  -sam SAM, --samtools SAM
                        samtools executable path if not in $PATH
  --salmon SALMON       Path to salmon executable, specify if salmon
                        quantification is desired
  --fusion_dist FUSION_DIST
                        minimium distance between separate read alignments on
                        the same chromosome to be considered a fusion,
                        otherwise no reads will be assumed to be fusions
  --mm2_args MM2_ARGS   additional minimap2 arguments when aligning reads
                        first-pass transcripts; separate args by commas, e.g.
                        --mm2_args=-I8g,--MD
  --quiet               Suppress progress statements from being printed
  --annotated_bed ANNOTATED_BED
                        annotation_reliant also requires a bedfile of
                        annotated isoforms; if this isn't provided, flair
                        collapse will generate the bedfile from the gtf.
                        eventually this argument will be removed
  --range RANGE         interval for which to collapse isoforms for, formatted
                        chromosome:coord1-coord2 or tab-delimited; if a range
                        is specified, then the aligned reads bam must be
                        specified with -r and the query must be a sorted,
                        bgzip-ed bed file

required named arguments:
  -q Q, --query Q       bed or psl file of aligned/corrected reads
  -g G, --genome G      FastA of reference genome
  -r R [R ...], --reads R [R ...]
                        FastA/FastQ files of raw reads, can specify multiple
                        files
