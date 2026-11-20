root = ..
include ${root}/defs.mk

##
# the test-installed targets set
#   use_installed_flair=yes
# to not add local targets
##
ifneq (${use_installed_flair},yes)
    export PATH := ${root}/bin:${PATH}
    export PYTHONPATH := ${root}/src:${PYTHONPATH}
endif

ifeq (${PYTHONWARNINGS},)
   # get resource and other warnings
   #   ResourceWarning: Enable tracemalloc to get the object allocation traceback
   # use
   #    PYTHONTRACEMALLOC=20 make
   # to find out were ResourceWarning are comming from
   # WANRING tracemalloc is really slow
   #
   # FIXME: gets lots of warnings by not closing files, and these break some
   # expect results.  Disable for now.
   #export PYTHONWARNINGS=always
endif

ifneq (${PYTHONTRACEMALLOC},)
    msg = WARNING: PYTHONTRACEMALLOC is set test will be very slow
    $(warning ${msg})
    $(warning ${msg})
    $(warning ${msg})
    $(warning ${msg})
endif



##
# all targets
##

# A test target says how it ended on stderr, so that the tail of a
# `make -j ... test -O >&log' says whether the run passed, rather than burying a
# failure in the middle of the log behind the tests that finished after it.
define run_tests
	${MAKE} $(1) && echo "===== TESTS PASSED =====" >&2 || { echo "===== TESTS FAILED =====" >&2; exit 1; }
endef

default:
	@echo "====="
	@echo "===== to see all targets: make help"
	@echo "====="
	${MAKE} test


help:
	@echo To run all tests, use: >&2
	@echo "   " make test >&2
	@echo To run tests without R depedencies, use: >&2
	@echo "   " make test-base >&2
	@echo To run tests using installed FLAIR rather than source tree, use: >&2
	@echo "   " make test-installed >&2
	@echo To run tests without R depedencies using installed FLAIR, use: >&2
	@echo "   " make test-base-installed >&2
	@echo Otherwise run individual tests by specifying test target, such as: >&2
	@echo "   " make transcriptome-tests >&2
	@echo to installed FLAIR for individual tests use the use_installed_flair option, such as: >&2
	@echo "   " make transcriptome-tests use_installed_flair=yes >&2
	@echo "Run only library tests" >&2
	@echo "   " make lib-tests >&2
	@echo "Write to an output directory of this run's own, so that runs for" >&2
	@echo "different python versions can go at the same time:" >&2
	@echo "   " make test outsuffix=-py3.14 >&2
	@echo "Build the generated inputs such runs share, before starting them:" >&2
	@echo "   " make test-inputs >&2
	@echo "Regenerate the committed whatshap inputs of the isoalleles tests:" >&2
	@echo "   " make regen-isoalleles-inputs >&2
	@echo "Remove output files:" >&2
	@echo "   " make clean >&2
	@echo Remove output and downloaded genomes files: >&2
	@echo "   " make real-clean >&2
	@echo run flake8: >&2
	@echo "   " make flake8 >&2
	@echo "build test case hub" >&2
	@echo "   " make test-hub >&2
	@echo "   " make make hub-check huburl=https://public.gi.ucsc.edu/~you/flair/test/hub/hub.txt >&2

# the two halves run one after the other: both include the help tests, and a
# parallel run of the two had the same target writing an output file while the
# other half's run of it was comparing that file
test:
	$(call run_tests, test-base)
	$(call run_tests, test-with-R)

test-base:
	${MAKE} test-alignments quantify-tests bed-to-gtf-tests \
		combine-tests diff-iso-usage-tests \
		plot-usage-tests transcriptome-tests fusion-tests variants-tests \
		spliceevents-tests partition-tests lib-tests test-help \
		test-protocol-unit scripts-tests

test-no-align:
	${MAKE} test-alignments quantify-tests bed-to-gtf-tests \
		combine-tests diff-iso-usage-tests \
		plot-usage-tests transcriptome-tests fusion-tests variants-tests \
		spliceevents-tests partition-tests lib-tests test-diffsplice-strands \
		scripts-tests

test-with-R:
	$(call run_tests, diffexp-tests diffsplice-tests)

test-installed:
	${MAKE} test use_installed_flair=yes

test-base-installed:
	${MAKE} test-base use_installed_flair=yes

# the help of each program, which test-base runs: test-main-help had listed a
# subcommand that no longer exists since it was only run on its own
test-help:
	${MAKE} test-main-help \
		test-quantify-help test-diffexp-help test-diffsplice-help \
		test-diff-iso-usage-help \
		test-plot-usage-help test-transcriptome-help test-fusion-help \
		test-spliceevents-help test-isoalleles-help test-doc-cli

# the committed option pages must match what the parsers currently define
test-doc-cli: mkdirs
	rm -rf $O/cli
	${root}/dev/bin/flair-gen-cli-doc $O/cli >/dev/null
	$(if ${HELP_FORMAT_313},diff -r ${root}/docs/source/cli $O/cli,$(call skip_help_cmp))

test-fast:
	$(call run_tests, test-combine-bed-only \
		test-diff-iso-usage test-plot-usage test-transcriptome \
		test-fusion)

##
# Non-python programs FLAIR runs.  test-env-only reruns the tests with only the
# conda environment and system directories in PATH; a developer's personal PATH
# will otherwise hide a program missing from misc/flair_dev_conda_env.yaml.
# Used for release testing.
##
external_progs = minimap2 samtools bedtools intron-prospector bedToBigBed longshot Rscript
env_only_path = ${CONDA_PREFIX}/bin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin

check-externals:
	@for prog in ${external_progs} ; do \
	    command -v $$prog >/dev/null || \
	        { echo "Error: required program not in PATH: $$prog; add it to misc/flair_dev_conda_env.yaml and misc/Dockerfile, then update the environment" >&2 ; exit 1 ; } ; \
	done

test-env-only:
	@[ -n "${CONDA_PREFIX}" ] || { echo "Error: activate the flair conda environment before running $@" >&2 ; exit 1 ; }
	PATH=${env_only_path} ${MAKE} check-externals test use_installed_flair=${use_installed_flair}



lint: flake8
flake8:
	cd .. && ${MAKE} flake8

##
# options for pytest
##
# keep pytest xdist from complaining about the terminal
export TERM=dumb

###
# pytest config
###
# override number of process with nproc=0 to run
# in main process.
# stdout not echo with mulitple processes
# keep small for now due to startup time with only a few tests
nproc = 0

# logging level
loglevel = INFO

# override with pytest_extra=-v or what ever
pytest_extra =

# Sometimes standard traceback (--tb=native) helts debug, less decoration,
# but lots pytest internal
# traceback = --tb=native

# --numprocesses is provided by pytest-xdist
PYTEST_CMD = ${PYTEST} --capture=no ${traceback} --log-cli-level=${loglevel} ${pytest_extra}

##
# Note: expected and output are sorted in cases where the order is not stable
##

##
# Output directory.  A run that must not disturb another one sets outsuffix, as
# the python version tests do with `outsuffix=-py3.14', giving each version its
# own output/ and letting the versions run at the same time.  ${TI} holds
# generated inputs that do not vary by python version, so it stays shared; `make
# test-inputs' builds them before parallel runs start.
##
outsuffix =

I = input
E = expected
O = output${outsuffix}
TI = tmp-input

# pytest tests that read a file make built get the directory from here
export FLAIR_TEST_OUTPUT_DIR = $O

# .pytest_cache is in the test directory, which a suffixed run shares with the
# runs going on beside it, so it keeps no cache
ifneq (${outsuffix},)
    pytest_extra += -p no:cacheprovider
endif

# reference data
GENOME_FA = ${TI}/genome.fa
# the test genome with all chromosomes concatenated into one, chr12_17_20
ONECHROM_GENOME_FA = ${TI}/genome.onechrom.fa
# the test genome with RefSeq's sequence names, for the RefSeq annotation test
REFSEQ_GENOME_FA = ${TI}/genome.refseq.fa

##
# basic: general tests that programs seem to work
##
BASIC_READS_FA = $I/basic.reads.fa
BASIC_ANNOTATION = $I/basic.annotation.gtf
BASIC_SHORTREAD_JUNCTIONS = $I/basic.shortread_junctions.tab

# FIXME: test are based on old test data, however flair quantify
# does not produce sufficient pvalues to work diffexpr
BASIC_COUNTS_MATRIX = $I/basic.counts_matrix.diffexp.tsv
BASIC_DIFFSPLICE_MATRIX = $I/basic.counts_matrix.diffsplice.tsv
# isoforms matching BASIC_DIFFSPLICE_MATRIX; test-collapse-annot.isoforms.bed
# is a different generation of the data and does not match the counts
BASIC_DIFFSPLICE_BED = $I/basic.isoforms.diffsplice.bed
STRANDS_DIFFSPLICE_BED = $I/strands.isoforms.diffsplice.bed
STRANDS_DIFFSPLICE_MATRIX = $I/strands.counts_matrix.diffsplice.tsv
COLLAPSE_ANNOT_ISOFORMS_BED = $I/test-collapse-annot.isoforms.bed

TOME_OUT_FA = $O/test-transcriptome.isoforms.fa
TOME_OUT_FAI = $O/test-transcriptome.isoforms.fa.fai

QUANTIFY_INPUT_READS = basic.reads.1 basic.reads.2 basic.reads.3 \
                       basic.reads.4 basic.reads.5 basic.reads.6
QUANTIFY_INPUT_BAMS = ${QUANTIFY_INPUT_READS:%=$O/%.genome_aligned.bam}

##
# Files generate by tests that are used by other tests.
# These enfoce files being create in the right order.
## 
BASIC_ALIGN_BAM = $O/basic.reads.bam
# basic reads with poly(A) tails removed, so single-exon reads have no poly(A)
# evidence for strand correction
BASIC_READS_NOPOLYA_FA = $O/basic.reads.nopolya.fa
NOPOLYA_ALIGN_BAM = $O/basic.reads.nopolya.bam
# poly(A)-trimmed basic reads with a poly(T) head added to every read, so single-exon
# reads from minus strand genes look like plus strand reads with a poly(A) tail
BASIC_READS_POLYT_FA = $O/basic.reads.polyt.fa
POLYT_ALIGN_BAM = $O/basic.reads.polyt.bam
# basic annotation plus a synthetic gene that shares one DDX5 intron
SHARED_JUNC_GENE_GTF = $I/shared_junction_gene.gtf
SHARED_JUNC_ANNOTATION = $O/basic.annotation.shared_junction.gtf
# basic.reads.bam with only one chr20 read, so the chr20 partition has reads but
# too few to support any isoform
ONECHR20_ALIGN_BAM = $O/basic.reads.onechr20.bam
# basic reads aligned to the one-chromosome genome
ONECHROM_ALIGN_BAM = $O/basic.reads.onechrom.bam

###
# seg1: single exon tests
###
SEG1_READS_FQ = $I/seg1.cdna_ont.fastq
SEG1_ANNOTATION = $I/seg1.gencode_v47.gtf

SEG1_ALIGN_BAM = $O/seg1.reads.bam

##
# mutant splice site, BED is truth
##
DDX5_READTHROUGH_FQ = $I/DDX5_readthrough.badread_40x.fastq
DDX5_READTHROUGH_BAM = $O/DDX5_readthrough.badread_40x.bam
DDX5_READTHROUGH_ANNOT = $I/DDX5_POLG2.gencode38.gtf
DDX5_READTHROUGH_SJ = $I/DDX5_readthrough.ip_junctions.bed
# DDX5 readthrough reads with all but 3 reads of the two supersets of the
# chr17:64500697-64506289 subset isoform dropped, so that subset outnumbers them
DDX5_BYSUPPORT_DROP_READS = $I/DDX5_readthrough.bysupport_drop_reads.txt
DDX5_BYSUPPORT_BAM = $O/DDX5_readthrough.bysupport.bam
# DDX5/POLG2 annotation plus the readthrough transcript as its own gene
DDX5_READTHROUGH_GENE_GTF = $I/DDX5_POLG2.readthrough_gene.gtf
DDX5_READTHROUGH_GENE_ANNOT = $O/DDX5_POLG2.gencode38.readthrough_gene.gtf

KRAS_MUT_SPLICE_FQ = $I/KRAS_sim_genomic_mut_driven_splicing.badread_sim_pb.fastq
KRAS_MUT_SPLICE_BAM = $O/KRAS_sim_genomic_mut_driven_splicing.badread_sim_pb.bam
KRAS_MUT_SPLICE_SJ_BED = $I/KRAS_sim_genomic_mut_driven_splicing.badread_sim_pb.sj.bed
KRAS_MUT_SPLICE_SJ_BB = $O/KRAS_sim_genomic_mut_driven_splicing.badread_sim_pb.sj.bb
KRAS_MUT_SPLICE_BED = $I/KRAS_sim_genomic_mut_driven_splicing.bed
KRAS_MUT_SPLICE_BB = $O/KRAS_sim_genomic_mut_driven_splicing.bb

##
# fusion tests; these are above ALIGNED_TEST_BAMS, which names them
##
ERGIC3_KRAS_SIM_FA_GZ = $O/ERGIC3_KRAS_sim_fusion.badread_20x.stranded.fastq.gz
FUSION_TEST_BAM = $O/ERGIC3_KRAS_sim_fusion.badread_20x.stranded.bam
# a 3-gene fusion, ERGIC3 (chr20 +) -- DDX5 (chr17 -) -- KRAS (chr12 -), joined at
# annotated exon boundaries and simulated with badread at 30x, reads oriented 5' to 3'
FUSION3_TEST_BAM = $O/ERGIC3_DDX5_KRAS_sim_fusion.badread_30x.stranded.bam
# the fusion test reads aligned to the test genome with RefSeq's sequence names
REFSEQ_FUSION_TEST_BAM = $O/ERGIC3_KRAS_sim_fusion.badread_20x.stranded.refseq.bam

COMBINE_MANIFEST = $I/combine.bed_only.manifest.txt
##
# Manifests name BAMs that make built, so the committed copy names the default
# output directory and a run using another one generates its own copy naming
# that, with the compressed copies following from it.
##
MANIFEST_NAMES = basic.bam_manifest.txt spliceevents.manifest.tsv \
	spliceevents.multi.manifest.tsv spliceevents.conditions.manifest.tsv
OUT_MANIFESTS = ${MANIFEST_NAMES:%=$O/%}
OUT_MANIFEST_GZS = $O/basic.bam_manifest.txt.gz $O/spliceevents.manifest.tsv.gz

${OUT_MANIFESTS}: $O/%: $I/% | mkdirs
	sed 's|output/|$O/|g' $< > $@.$$$$.tmp && mv -f $@.$$$$.tmp $@

${OUT_MANIFEST_GZS}: %.gz: %
	gzip -c $< > $@.$$$$.tmp && mv -f $@.$$$$.tmp $@

QUANTIFY_MANIFEST = $O/basic.bam_manifest.txt
SPLICEEVENTS_MANIFEST = $O/spliceevents.manifest.tsv
# six samples of the basic reads, for outliers: four random 60% samples, one with every
# read including KRAS's skipped exon (chr12:25215436-25215560) and a third of the
# rest, one with every read retaining a chr17 intron (64500326-64500548) and a third
# of the rest
SPLICEEVENTS_MULTI_MANIFEST = $O/spliceevents.multi.manifest.tsv
SPLICEEVENTS_MULTI_BAMS = $O/spliceevents.s1.bam $O/spliceevents.s2.bam $O/spliceevents.s3.bam $O/spliceevents.s4.bam \
	$O/spliceevents.exon_inclusion.bam $O/spliceevents.intron_retention.bam

TOME_ISOFORM_BED = $E/test-transcriptome.isoforms.bed
TOME_READ_MAP = $E/test-transcriptome.isoforms.read.map.tsv

# The help test results vary on screen width.  This forces them to be consistent.
# It might seem weird to see if help works, but some it did break in the past.
export COLUMNS = 80

##
# Compressed copies of inputs, built on demand, so the compression tests read the
# same data as the plain tests.  Every text input flair reads may be gzipped; BAM
# is already compressed and never is.
##
${TI}/%.gz: $I/%
	@mkdir -p ${TI}
	gzip -c $< > $@.$$$$.tmp && mv -f $@.$$$$.tmp $@

##
# A bgzip genome, indexed.  The genome is the one input that cannot be plain gzip:
# the sequence reader needs random access, which only bgzip gives it.
##
BGZIP_GENOME_FA = ${TI}/genome.bgz.fa.gz

${BGZIP_GENOME_FA}: ${GENOME_FA}
	bgzip -c $< > $@.$$$$.tmp && mv -f $@.$$$$.tmp $@
	samtools faidx $@

${TI}/%.gz: $E/%
	@mkdir -p ${TI}
	gzip -c $< > $@.$$$$.tmp && mv -f $@.$$$$.tmp $@

##
# compare output, expected with result
##
define cmp_expect
	diff $E/$(1) $O/$(1)
endef

##
# compare a captured log, whose records are stamped with the time they were logged,
# which is not reproducible: the stamps are dropped first
##
define cmp_expect_log
	sed -i -r -e 's/^[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2} //' $O/$(1)
	diff $E/$(1) $O/$(1)
endef

##
# compare --help output.  argparse writes "-o, --output OUTPUT" from python 3.13
# on, where 3.12 wrote "-o OUTPUT, --output OUTPUT", and wraps usage lines at
# different points.  The expected files are in the 3.13 form, so an older python
# skips the comparison rather than reporting a difference that is not a defect.
##
HELP_FORMAT_313 := $(shell python3 -c 'import sys; sys.stdout.write("yes" if sys.version_info >= (3, 13) else "")')
PYTHON_VERSION := $(shell python3 -c 'import platform; print(platform.python_version())')

define skip_help_cmp
	@echo "Warning: $@ not compared: python ${PYTHON_VERSION} writes argparse help in the pre-3.13 format, the expected output is in the 3.13 format" >&2
endef

define cmp_expect_help
	$(if ${HELP_FORMAT_313},diff $E/$(1) $O/$(1),$(call skip_help_cmp))
endef

##
# compare a compression test's output against the expected file of the plain test,
# which is what it must reproduce from gzipped input
##
define cmp_expect_named
	@test -e $O/${1} || (echo Error: output file missing: $O/${1}; exit 1)
	diff $E/${2} $O/${1}
endef

##
# compare a compressed output against the expected file of the plain test, reading
# it through zcat and leaving it compressed
##
define cmp_expect_gz
	@test -e $O/${1}.gz || (echo Error: output file missing: $O/${1}.gz; exit 1)
	diff $E/${2} <(zcat $O/${1}.gz)
endef

##
# a compressed output that has no expected file of its own: it must be there and
# must decompress
##
define check_gz
	@test -e $O/${1}.gz || (echo Error: output file missing: $O/${1}.gz; exit 1)
	zcat $O/${1}.gz > /dev/null
endef

##
# a compressed FASTA is bgzip rather than plain gzip so that it can be indexed,
# which is what this checks; samtools faidx refuses a plain gzip file
##
define check_fasta_indexable
	samtools faidx $O/${1}.gz
endef

##
# compare an output file against another output file, for tests that must produce
# the same result as an earlier test from differently arranged input
##
# a read map, merged by isoform and sorted (bin/read-map-sort): quantify writes a row
# when its batch finishes, so the row order follows the batching that --threads sets
define cmp_expect_read_map
	@test -e $O/${1} || (echo Error: output file missing: $O/${1}; exit 1)
	diff <(bin/read-map-sort < $E/${1}) <(bin/read-map-sort < $O/${1})
endef

# a read map that must reproduce another test's, sorted the same way
define cmp_expect_read_map_named
	@test -e $O/${1} || (echo Error: output file missing: $O/${1}; exit 1)
	diff <(bin/read-map-sort < $E/${2}) <(bin/read-map-sort < $O/${1})
endef

define cmp_expect_sort_as
	@test -e $O/${1} || (echo Error: output file missing: $O/${1}; exit 1)
	diff <(sort $O/${2} | bin/tsv-round-numeric) <(sort $O/${1} | bin/tsv-round-numeric)
endef

##
# as above, but comparing only the leading ids and the trailing statistics.  The
# DRIMSeq result also carries one column per sample, which a test that reorders the
# input columns permutes by construction.
##
define cmp_expect_sort_stats_as
	@test -e $O/${1} || (echo Error: output file missing: $O/${1}; exit 1)
	diff <(sort $O/${2} | awk -F'\t' -v OFS='\t' '{print $$1, $$2, $$(NF-1), $$NF}' | bin/tsv-round-numeric) \
	     <(sort $O/${1} | awk -F'\t' -v OFS='\t' '{print $$1, $$2, $$(NF-1), $$NF}' | bin/tsv-round-numeric)
endef

###
# MAIN COMMAND
###

test-main-help: mkdirs
	flair --help | fgrep -v Version> $O/$@.out
	$(call cmp_expect_help,$@.out)

###
# reads aligned to the genome, the input of the tests
###

# align reads $(1) to genome $(2) with minimap2, as users are asked to (spliced, MD
# tags, no secondary alignments), with further minimap2 options $(3), to a sorted,
# indexed $@
define minimap2_genome_align
	@mkdir -p $(dir $@)
	minimap2 -ax splice --MD --secondary=no -t 4 $(3) $(2) $(1) | samtools sort -o $@.tmp.bam -
	mv -f $@.tmp.bam $@
	samtools index $@
endef

# The test BAMs' alignments, checked against expected/<bam>.alignments.tsv first, so
# a minimap2 that aligns the test reads differently fails here rather than as changed
# output of the tests using them.  A mismatching minimap2 version alone is only
# reported.  update-alignment-tests takes the current alignments as expected, after
# checking that the changes are wanted.
ALIGNED_TEST_BAMS = ${BASIC_ALIGN_BAM} ${SEG1_ALIGN_BAM} ${NOPOLYA_ALIGN_BAM} ${POLYT_ALIGN_BAM} \
	${ONECHROM_ALIGN_BAM} ${FUSION_TEST_BAM} ${FUSION3_TEST_BAM} ${REFSEQ_FUSION_TEST_BAM} \
	${KRAS_MUT_SPLICE_BAM} ${DDX5_READTHROUGH_BAM}

test-alignments: mkdirs ${ALIGNED_TEST_BAMS}
	@echo "minimap2 $$(minimap2 --version), expected alignments are from $$(cat $E/minimap2.version)" >&2
	@status=0; \
	for bam in ${ALIGNED_TEST_BAMS}; do \
	    name=$$(basename $$bam .bam); \
	    bin/bam-alignment-summary $$bam > $O/$$name.alignments.tsv; \
	    if ! cmp -s $E/$$name.alignments.tsv $O/$$name.alignments.tsv; then \
	        echo "Error: alignments of $$name differ: diff $E/$$name.alignments.tsv $O/$$name.alignments.tsv" >&2; \
	        status=1; \
	    fi; \
	done; \
	exit $$status

update-alignment-tests: ${ALIGNED_TEST_BAMS}
	for bam in ${ALIGNED_TEST_BAMS}; do \
	    bin/bam-alignment-summary $$bam > $E/$$(basename $$bam .bam).alignments.tsv; \
	done
	minimap2 --version > $E/minimap2.version

${BASIC_ALIGN_BAM}: ${BASIC_READS_FA} ${GENOME_FA}
	$(call minimap2_genome_align,${BASIC_READS_FA},${GENOME_FA},-s 80 -G 200k)

${SEG1_ALIGN_BAM}: ${SEG1_READS_FQ} ${GENOME_FA}
	$(call minimap2_genome_align,${SEG1_READS_FQ},${GENOME_FA},-s 80 -G 200k)

${BASIC_READS_NOPOLYA_FA}: ${BASIC_READS_FA} bin/trim_polya_reads
	@mkdir -p $O
	bin/trim_polya_reads ${BASIC_READS_FA} $@.tmp
	mv -f $@.tmp $@

# basic reads with poly(A) tails trimmed
${NOPOLYA_ALIGN_BAM}: ${BASIC_READS_NOPOLYA_FA} ${GENOME_FA}
	$(call minimap2_genome_align,${BASIC_READS_NOPOLYA_FA},${GENOME_FA},-s 80 -G 200k)

${BASIC_READS_POLYT_FA}: ${BASIC_READS_NOPOLYA_FA}
	awk '/^>/ {print; next} {print "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" $$0}' ${BASIC_READS_NOPOLYA_FA} > $@.tmp
	mv -f $@.tmp $@

${POLYT_ALIGN_BAM}: ${BASIC_READS_POLYT_FA} ${GENOME_FA}
	$(call minimap2_genome_align,${BASIC_READS_POLYT_FA},${GENOME_FA},-s 80 -G 200k)

${ONECHR20_ALIGN_BAM}: ${BASIC_ALIGN_BAM}
	(samtools view -h ${BASIC_ALIGN_BAM} chr12 chr17 && samtools view -F 0x900 ${BASIC_ALIGN_BAM} chr20 | head -n 1) | samtools view -b -o $@.tmp.bam -
	mv -f $@.tmp.bam $@
	samtools index $@

${SHARED_JUNC_ANNOTATION}: ${BASIC_ANNOTATION} ${SHARED_JUNC_GENE_GTF}
	@mkdir -p $O
	cat ${BASIC_ANNOTATION} ${SHARED_JUNC_GENE_GTF} > $@.tmp
	mv -f $@.tmp $@

# regenerate the expected files a family already has from the output of a run: a
# name no expected file has is not created, since what a test compares is in its recipe
define update_expected
	for exp in $E/$(1)*; do \
	    out=$O/$$(basename $$exp); \
	    if test -f $$exp -a -f $$out; then cp -f $$out $$exp; fi; \
	done
endef

# macros to do comparions with expected results
define diff_bed_gtf
	$(call cmp_expect,$@.isoforms.bed)
	$(call cmp_expect,$@.isoforms.gtf)
endef

###
# FLAIR COMBINE
###

define combine_diff
	$(call cmp_expect,$@.combined.isoforms.bed)
	$(call cmp_expect,$@.combined.isoforms.counts.tsv)
	$(call cmp_expect,$@.combined.isoforms.map.tsv)
	$(call cmp_expect,$@.combined.genes.map.tsv)
endef


combine-tests: test-combine test-combine-files test-combine-prefixes test-combine-nose test-combine-10pct test-combine-multends test-combine-multends-10pct \
	test-combine-conflict test-combine-conflict-gtf test-combine-isoalleles test-combine-compress \
	test-combine-combined

update-combine-tests:
	$(call update_expected,test-combine)

test-combine: mkdirs ${GENOME_FA}
	flair combine --manifest $(COMBINE_MANIFEST) --genome ${GENOME_FA} -o $O/$@
	$(combine_diff)
	$(call cmp_expect,$@.combined.isoforms.fa)
	$(call cmp_expect,$@.combined.isoforms.gtf)

# the same eight BEDs test-combine-files uses, gzipped
COMBINE_BEDS_GZ = ${TI}/test-transcriptome.isoforms.bed.gz ${TI}/test-transcriptome-mult-ends.isoforms.bed.gz \
	${TI}/test-transcriptome-mutant-splice.isoforms.bed.gz ${TI}/test-transcriptome-readthrough.isoforms.bed.gz \
	${TI}/test-fusion.fusions.isoforms.bed.gz ${TI}/test-transcriptome-seg1-badannot.isoforms.bed.gz \
	${TI}/test-transcriptome-seg1-badannot-multends.isoforms.bed.gz \
	${TI}/test-transcriptome-seg1-badannot-multends-chromshifted.isoforms.bed.gz

test-combine-compress: mkdirs ${GENOME_FA} ${COMBINE_BEDS_GZ}
	flair -z combine --genome ${GENOME_FA} -o $O/$@ --isoform_beds $(shell echo $(COMBINE_BEDS_GZ) | tr ' ' ',')
	$(call cmp_expect_gz,$@.combined.isoforms.counts.tsv,test-combine-files.combined.isoforms.counts.tsv)
	$(call cmp_expect_gz,$@.combined.isoforms.bed,test-combine-files.combined.isoforms.bed)
	$(call check_gz,$@.combined.isoforms.gtf)
	$(call check_gz,$@.combined.isoforms.fa)
	$(call check_fasta_indexable,$@.combined.isoforms.fa)

test-combine-files: mkdirs ${GENOME_FA}
	flair combine --genome ${GENOME_FA} -o $O/$@ --isoform_beds expected/test-transcriptome.isoforms.bed,expected/test-transcriptome-mult-ends.isoforms.bed,expected/test-transcriptome-mutant-splice.isoforms.bed,expected/test-transcriptome-readthrough.isoforms.bed,expected/test-fusion.fusions.isoforms.bed,expected/test-transcriptome-seg1-badannot.isoforms.bed,expected/test-transcriptome-seg1-badannot-multends.isoforms.bed,input/test-transcriptome-seg1-badannot-multends-chromshifted.isoforms.bed
	$(combine_diff)
	$(call cmp_expect_named,$@.combined.isoforms.fa,test-combine.combined.isoforms.fa)
	$(call cmp_expect_named,$@.combined.isoforms.gtf,test-combine.combined.isoforms.gtf)

# the same eight transcriptomes as test-combine-files, named by prefix with their
# directories, comma separated and by repeating the option
test-combine-prefixes: mkdirs ${GENOME_FA}
	flair combine --genome ${GENOME_FA} -o $O/$@ \
		--prefixes expected/test-transcriptome,expected/test-transcriptome-mult-ends \
		--prefixes expected/test-transcriptome-mutant-splice,expected/test-transcriptome-readthrough, \
		--prefixes expected/test-fusion.fusions \
		--prefixes expected/test-transcriptome-seg1-badannot,expected/test-transcriptome-seg1-badannot-multends \
		--prefixes input/test-transcriptome-seg1-badannot-multends-chromshifted
	diff $E/test-combine-files.combined.isoforms.bed $O/$@.combined.isoforms.bed
	diff $E/test-combine-files.combined.isoforms.counts.tsv $O/$@.combined.isoforms.counts.tsv
	$(call cmp_expect_named,$@.combined.genes.map.tsv,test-combine-files.combined.genes.map.tsv)
	$(call cmp_expect_named,$@.combined.isoforms.fa,test-combine.combined.isoforms.fa)
	$(call cmp_expect_named,$@.combined.isoforms.gtf,test-combine.combined.isoforms.gtf)
	$(call cmp_expect_named,$@.combined.isoforms.map.tsv,test-combine-files.combined.isoforms.map.tsv)

test-combine-nose: mkdirs ${GENOME_FA}
	flair combine --manifest $(COMBINE_MANIFEST) --genome ${GENOME_FA} -o $O/$@ --remove_single_exon
	$(combine_diff)
	$(call cmp_expect,$@.combined.isoforms.fa)
	$(call cmp_expect,$@.combined.isoforms.gtf)

test-combine-10pct: mkdirs ${GENOME_FA}
	flair combine --manifest $(COMBINE_MANIFEST) --genome ${GENOME_FA} -o $O/$@ --min_percent_usage 10
	$(combine_diff)
	$(call cmp_expect,$@.combined.isoforms.fa)
	$(call cmp_expect,$@.combined.isoforms.gtf)

test-combine-multends: mkdirs ${GENOME_FA}
	flair combine --manifest $(COMBINE_MANIFEST) --genome ${GENOME_FA} -o $O/$@ --max_ends 10
	$(combine_diff)
	$(call cmp_expect,$@.combined.isoforms.fa)
	$(call cmp_expect,$@.combined.isoforms.gtf)

test-combine-multends-10pct: mkdirs ${GENOME_FA}
	flair combine --manifest $(COMBINE_MANIFEST) --genome ${GENOME_FA} -o $O/$@ --max_ends 10 --min_percent_usage 10
	$(combine_diff)
	$(call cmp_expect,$@.combined.isoforms.fa)
	$(call cmp_expect,$@.combined.isoforms.gtf)

# a sample whose isoform of ENST00000348547 was given KRAS: without the annotation,
# an error naming the samples' genes, and with it, the isoform gets its own gene
test-combine-conflict: mkdirs ${GENOME_FA}
	! flair combine --manifest $I/combine.conflict.manifest.txt --genome ${GENOME_FA} -o $O/$@ 2> $O/$@.err
	grep -q 'test-combine-conflict: ENSG00000133703.12; test-transcriptome: ENSG00000125991.19' $O/$@.err

test-combine-conflict-gtf: mkdirs ${GENOME_FA}
	flair combine --manifest $I/combine.conflict.manifest.txt --genome ${GENOME_FA} -f $(BASIC_ANNOTATION) -o $O/$@
	$(combine_diff)
	$(call cmp_expect,$@.combined.aaseq.tsv)
	$(call cmp_expect,$@.combined.isoforms.fa)
	$(call cmp_expect,$@.combined.isoforms.gtf)

# isoforms of alleles, from flair isoalleles: combined into their isoforms,
# without allele information, with a warning.  read_support is per allele,
# frac_support the isoform's, so --min_percent_usage treats an isoform's alleles alike
test-combine-isoalleles: mkdirs ${GENOME_FA}
	flair combine --isoform_beds $I/combine.isoalleles.bed --genome ${GENOME_FA} --min_percent_usage 1 -o $O/$@ 2> $O/$@.err
	grep -q 'flair combine drops the allele information' $O/$@.err
	$(combine_diff)
	$(call cmp_expect,$@.combined.isoforms.fa)
	$(call cmp_expect,$@.combined.isoforms.gtf)

# combine's own output names every sample an isoform came from and contains their
# summed reads: combining it again would give each of those samples all the reads
test-combine-combined: mkdirs ${GENOME_FA}
	! flair combine --isoform_beds $E/test-combine.combined.isoforms.bed --genome ${GENOME_FA} -o $O/$@ 2> $O/$@.err
	grep -q 'isoform FLT00000001 names 3 samples' $O/$@.err

###
# FLAIR QUANTIFY
###
quantify-tests: test-quantify test-quantify-normends test-quantify-onethread-keep test-quantify-compress \
	test-quantify-bgzip-genome test-quantify-fusion test-quantify-fusion-3gene test-quantify-help

test-quantify: mkdirs ${TOME_OUT_FAI} ${QUANTIFY_INPUT_BAMS} ${GENOME_FA} ${QUANTIFY_MANIFEST}
	flair quantify --manifest $(QUANTIFY_MANIFEST) --genome ${GENOME_FA} --isoform_bed $(TOME_ISOFORM_BED) --generate_map -o $O/$@
	$(call cmp_expect,$@.isoforms.counts.tsv)
	$(call cmp_expect,$@.isoforms.sample_info.tsv)
	$(call cmp_expect,$@.isoforms.cpm.tsv)
	$(call cmp_expect_read_map,$@.isoforms.A1.read.map.tsv)
	$(call cmp_expect_read_map,$@.isoforms.B3.read.map.tsv)
	$(call cmp_expect_read_map,$@.isoforms.A2.read.map.tsv)
	$(call cmp_expect_read_map,$@.isoforms.A3.read.map.tsv)
	$(call cmp_expect_read_map,$@.isoforms.B1.read.map.tsv)
	$(call cmp_expect_read_map,$@.isoforms.B2.read.map.tsv)

test-quantify-normends: mkdirs ${TOME_OUT_FAI} ${QUANTIFY_INPUT_BAMS} ${GENOME_FA} ${QUANTIFY_MANIFEST}
	flair quantify --manifest $(QUANTIFY_MANIFEST) --genome ${GENOME_FA} --isoform_bed $(TOME_ISOFORM_BED) --generate_map -o $O/$@ --normalize_ends
	$(call cmp_expect,$@.isoforms.counts.tsv)
	$(call cmp_expect,$@.isoforms.sample_info.tsv)
	$(call cmp_expect,$@.isoforms.cpm.tsv)
	$(call cmp_expect_read_map,$@.isoforms.A1.read.map.tsv)
	$(call cmp_expect_read_map_named,$@.isoforms.A2.read.map.tsv,test-quantify.isoforms.A2.read.map.tsv)
	$(call cmp_expect_read_map_named,$@.isoforms.A3.read.map.tsv,test-quantify.isoforms.A3.read.map.tsv)
	$(call cmp_expect_read_map_named,$@.isoforms.B1.read.map.tsv,test-quantify.isoforms.B1.read.map.tsv)
	$(call cmp_expect_read_map_named,$@.isoforms.B2.read.map.tsv,test-quantify.isoforms.B2.read.map.tsv)
	$(call cmp_expect_read_map,$@.isoforms.B3.read.map.tsv)

# one thread, without the pool, and the intermediate files kept: the same counts as test-quantify
test-quantify-onethread-keep: mkdirs ${TOME_OUT_FAI} ${QUANTIFY_INPUT_BAMS} ${GENOME_FA} ${QUANTIFY_MANIFEST}
	rm -rf $O/$@.tmp
	flair quantify --manifest $(QUANTIFY_MANIFEST) --genome ${GENOME_FA} --isoform_bed $(TOME_ISOFORM_BED) --generate_map -o $O/$@ \
	   -t 1 --temp_dir $O/$@.tmp --keep_intermediate
	diff $E/test-quantify.isoforms.counts.tsv $O/$@.isoforms.counts.tsv
	diff <(bin/read-map-sort < $E/test-quantify.isoforms.A1.read.map.tsv) <(bin/read-map-sort < $O/$@.isoforms.A1.read.map.tsv)
	ls $O/$@.tmp/*.intermediate/*/isoforms.fa > /dev/null
	$(call cmp_expect_read_map_named,$@.isoforms.A2.read.map.tsv,test-quantify.isoforms.A2.read.map.tsv)
	$(call cmp_expect_read_map_named,$@.isoforms.A3.read.map.tsv,test-quantify.isoforms.A3.read.map.tsv)
	$(call cmp_expect_read_map_named,$@.isoforms.B1.read.map.tsv,test-quantify.isoforms.B1.read.map.tsv)
	$(call cmp_expect_read_map_named,$@.isoforms.B2.read.map.tsv,test-quantify.isoforms.B2.read.map.tsv)
	$(call cmp_expect_read_map,$@.isoforms.B3.read.map.tsv)
	$(call cmp_expect_named,$@.isoforms.cpm.tsv,test-quantify.isoforms.cpm.tsv)
	$(call cmp_expect_named,$@.isoforms.sample_info.tsv,test-quantify.isoforms.sample_info.tsv)

# combine's isoforms with a two-piece fusion, on opposite strands, with the fusion's simulated reads as a
# sample: the fusion gets flair fusion's counts (test-fusion), and its reads no other isoform's
test-quantify-fusion: mkdirs ${TOME_OUT_FAI} ${QUANTIFY_INPUT_BAMS} ${GENOME_FA} ${FUSION_TEST_BAM} ${QUANTIFY_MANIFEST}
	( cat $(QUANTIFY_MANIFEST) && printf 'fusion\tF\tb0\t%s\n' ${FUSION_TEST_BAM} ) > $O/$@.manifest.tsv
	flair quantify --manifest $O/$@.manifest.tsv --genome ${GENOME_FA} --isoform_bed $E/test-combine-files.combined.isoforms.bed \
	   --generate_map -o $O/$@
	$(call cmp_expect,$@.isoforms.counts.tsv)
	$(call cmp_expect,$@.isoforms.cpm.tsv)
	$(call cmp_expect_read_map,$@.isoforms.fusion.read.map.tsv)
	$(call cmp_expect_read_map,$@.isoforms.A1.read.map.tsv)
	$(call cmp_expect_read_map,$@.isoforms.A2.read.map.tsv)
	$(call cmp_expect_read_map,$@.isoforms.A3.read.map.tsv)
	$(call cmp_expect_read_map,$@.isoforms.B1.read.map.tsv)
	$(call cmp_expect_read_map,$@.isoforms.B2.read.map.tsv)
	$(call cmp_expect_read_map,$@.isoforms.B3.read.map.tsv)
	$(call cmp_expect,$@.isoforms.sample_info.tsv)

# a three-piece fusion: flair fusion's counts (test-fusion-3gene)
test-quantify-fusion-3gene: mkdirs ${GENOME_FA} ${FUSION3_TEST_BAM}
	printf 'fusion3\tF\tb0\t%s\n' ${FUSION3_TEST_BAM} > $O/$@.manifest.tsv
	flair quantify --manifest $O/$@.manifest.tsv --genome ${GENOME_FA} --isoform_bed $E/test-fusion-3gene.fusions.isoforms.bed \
	   -o $O/$@
	$(call cmp_expect,$@.isoforms.counts.tsv)
	$(call cmp_expect,$@.isoforms.cpm.tsv)
	$(call cmp_expect,$@.isoforms.sample_info.tsv)

# every text input gzipped: the isoform BED and the manifest
test-quantify-compress: mkdirs ${TOME_OUT_FAI} ${QUANTIFY_INPUT_BAMS} ${GENOME_FA} \
		${TI}/test-transcriptome.isoforms.bed.gz $O/basic.bam_manifest.txt.gz
	flair --compress quantify --manifest $O/basic.bam_manifest.txt.gz --genome ${GENOME_FA} \
	   --isoform_bed ${TI}/test-transcriptome.isoforms.bed.gz --generate_map -o $O/$@
	$(call cmp_expect_gz,$@.isoforms.counts.tsv,test-quantify.isoforms.counts.tsv)
	$(call cmp_expect_gz,$@.isoforms.sample_info.tsv,test-quantify.isoforms.sample_info.tsv)
	$(call cmp_expect_gz,$@.isoforms.cpm.tsv,test-quantify.isoforms.cpm.tsv)
	@# read maps are bgzip, which flair isoalleles reads by offset, holding the same reads
	htsfile $O/$@.isoforms.A1.read.map.tsv.gz | grep -q BGZF
	diff <(bin/read-map-sort < $E/test-quantify.isoforms.A1.read.map.tsv) <(gunzip -c $O/$@.isoforms.A1.read.map.tsv.gz | bin/read-map-sort)

# the genome bgzipped and indexed, every other input plain
test-quantify-bgzip-genome: mkdirs ${TOME_OUT_FAI} ${QUANTIFY_INPUT_BAMS} ${BGZIP_GENOME_FA} ${QUANTIFY_MANIFEST}
	flair quantify --manifest $(QUANTIFY_MANIFEST) --genome ${BGZIP_GENOME_FA} \
	   --isoform_bed $(TOME_ISOFORM_BED) --generate_map -o $O/$@
	$(call cmp_expect_named,$@.isoforms.counts.tsv,test-quantify.isoforms.counts.tsv)
	$(call cmp_expect_named,$@.isoforms.cpm.tsv,test-quantify.isoforms.cpm.tsv)
	$(call cmp_expect_read_map,$@.isoforms.A1.read.map.tsv)
	$(call cmp_expect_read_map_named,$@.isoforms.A2.read.map.tsv,test-quantify.isoforms.A2.read.map.tsv)
	$(call cmp_expect_read_map_named,$@.isoforms.A3.read.map.tsv,test-quantify.isoforms.A3.read.map.tsv)
	$(call cmp_expect_read_map_named,$@.isoforms.B1.read.map.tsv,test-quantify.isoforms.B1.read.map.tsv)
	$(call cmp_expect_read_map_named,$@.isoforms.B2.read.map.tsv,test-quantify.isoforms.B2.read.map.tsv)
	$(call cmp_expect_read_map,$@.isoforms.B3.read.map.tsv)
	$(call cmp_expect_named,$@.isoforms.sample_info.tsv,test-quantify.isoforms.sample_info.tsv)

test-quantify-help: mkdirs
	flair quantify --help > $O/$@.out
	$(call cmp_expect_help,$@.out)

###
# FLAIR DIFFEXP
###
diffexp-tests: test-diffexp test-diffexp-compress test-diffexp-colorder test-diffexp-spliceevents test-diffexp-help

# the conditions are named, so which one is the reference does not depend on where
# their columns fall in the counts matrix
test-diffexp : mkdirs
	flair diffexp --counts_matrix $(BASIC_COUNTS_MATRIX) -o $O/$@ --min_expression 1 \
	   --condition_a MCF7 --condition_b A549
	$(call cmp_expect,$@.genes_deseq2_MCF7_v_A549.tsv)
	$(call cmp_expect,$@.isoforms_deseq2_MCF7_v_A549.tsv)
	$(call cmp_expect,$@.isoforms_drimseq_MCF7_v_A549.tsv)
	$(call cmp_expect,$@.isoforms_sig_exp_change_norm_by_gene.tsv)

# same counts, columns interleaved rather than grouped by condition, compared against
# the same expected results.  Taking the conditions from the first and last column
# got this wrong, and with an order like A,B,B,A filtered one condition twice.
test-diffexp-compress: mkdirs ${TI}/basic.counts_matrix.diffexp.tsv.gz
	flair -z diffexp --counts_matrix ${TI}/basic.counts_matrix.diffexp.tsv.gz -o $O/$@ --min_expression 1 \
	   --condition_a MCF7 --condition_b A549
	$(call cmp_expect_gz,$@.genes_deseq2_MCF7_v_A549.tsv,test-diffexp.genes_deseq2_MCF7_v_A549.tsv)
	$(call cmp_expect_gz,$@.isoforms_deseq2_MCF7_v_A549.tsv,test-diffexp.isoforms_deseq2_MCF7_v_A549.tsv)

test-diffexp-colorder : mkdirs test-diffexp
	${PYTHON} bin/tsv-condition-interleave $(BASIC_COUNTS_MATRIX) $O/$@.counts_matrix.tsv
	flair diffexp --counts_matrix $O/$@.counts_matrix.tsv -o $O/$@ --min_expression 1 \
	   --condition_a MCF7 --condition_b A549
	$(call cmp_expect_sort_as,$@.genes_deseq2_MCF7_v_A549.tsv,test-diffexp.genes_deseq2_MCF7_v_A549.tsv)
	$(call cmp_expect_sort_as,$@.isoforms_deseq2_MCF7_v_A549.tsv,test-diffexp.isoforms_deseq2_MCF7_v_A549.tsv)
	$(call cmp_expect_sort_stats_as,$@.isoforms_drimseq_MCF7_v_A549.tsv,test-diffexp.isoforms_drimseq_MCF7_v_A549.tsv)
	$(call cmp_expect_named,$@.isoforms_sig_exp_change_norm_by_gene.tsv,test-diffexp.isoforms_sig_exp_change_norm_by_gene.tsv)

# the basic counts as a flair spliceevents table, each gene an event, the last isoform of
# a gene of three or more its other row: usage alone is tested, no DESeq2, and other rows
# count toward their event but are not reported.  The other features' DRIMSeq results
# are those of test-diffexp
test-diffexp-spliceevents: mkdirs
	flair diffexp --counts_matrix $I/basic.diffexp.spliceevents.counts.tsv -o $O/$@ --min_expression 1 \
	   --condition_a MCF7 --condition_b A549
	test ! -e $O/$@.genes_deseq2_MCF7_v_A549.tsv
	! grep -q '|other' $O/$@.isoforms_drimseq_MCF7_v_A549.tsv $O/$@.isoforms_sig_exp_change_norm_by_gene.tsv
	$(call cmp_expect,$@.isoforms_drimseq_MCF7_v_A549.tsv)
	$(call cmp_expect,$@.isoforms_sig_exp_change_norm_by_gene.tsv)

test-diffexp-help: mkdirs
	flair diffexp --help > $O/$@.out
	$(call cmp_expect_help,$@.out)

###
# FLAIR DIFFSPLICE
###
diffsplice-tests: test-diffsplice test-diffsplice-cond test-diffsplice-strands test-diffsplice-help test-diffsplice-compress

test-diffsplice: mkdirs
	flair diffsplice --isoform_bed $(BASIC_DIFFSPLICE_BED) --counts_matrix $(BASIC_DIFFSPLICE_MATRIX) --test -o $O/$@
	$(call cmp_expect,$@.alt3.events.quant.tsv)
	$(call cmp_expect,$@.alt5.events.quant.tsv)
	$(call cmp_expect,$@.es.events.quant.tsv)
	$(call cmp_expect,$@.ir.events.quant.tsv)
	$(call cmp_expect,$@.drimseq_alt3_A_v_B.tsv)
	$(call cmp_expect,$@.drimseq_alt5_A_v_B.tsv)
	$(call cmp_expect,$@.drimseq_es_A_v_B.tsv)
	$(call cmp_expect,$@.drimseq_ir_A_v_B.tsv)

test-diffsplice-compress: mkdirs ${TI}/basic.isoforms.diffsplice.bed.gz ${TI}/basic.counts_matrix.diffsplice.tsv.gz
	flair -z diffsplice --isoform_bed ${TI}/basic.isoforms.diffsplice.bed.gz \
	   --counts_matrix ${TI}/basic.counts_matrix.diffsplice.tsv.gz --test -o $O/$@
	$(call cmp_expect_gz,$@.alt3.events.quant.tsv,test-diffsplice.alt3.events.quant.tsv)
	$(call cmp_expect_gz,$@.es.events.quant.tsv,test-diffsplice.es.events.quant.tsv)
	$(call cmp_expect_gz,$@.drimseq_es_A_v_B.tsv,test-diffsplice.drimseq_es_A_v_B.tsv)
	$(call cmp_expect_gz,$@.drimseq_ir_A_v_B.tsv,test-diffsplice.drimseq_ir_A_v_B.tsv)

# naming the two conditions gives what taking them in sorted order gives.  The
# comparisons used to name test-diffsplice on both sides, so they re-checked that
# test's output and never looked at this one's
test-diffsplice-cond: mkdirs test-diffsplice
	flair diffsplice --isoform_bed $(BASIC_DIFFSPLICE_BED) --counts_matrix $(BASIC_DIFFSPLICE_MATRIX) --test -o $O/$@ --condition_a A --condition_b B
	$(call cmp_expect_sort_as,$@.alt3.events.quant.tsv,test-diffsplice.alt3.events.quant.tsv)
	$(call cmp_expect_sort_as,$@.alt5.events.quant.tsv,test-diffsplice.alt5.events.quant.tsv)
	$(call cmp_expect_sort_as,$@.es.events.quant.tsv,test-diffsplice.es.events.quant.tsv)
	$(call cmp_expect_sort_as,$@.ir.events.quant.tsv,test-diffsplice.ir.events.quant.tsv)
	$(call cmp_expect_sort_as,$@.drimseq_alt3_A_v_B.tsv,test-diffsplice.drimseq_alt3_A_v_B.tsv)
	$(call cmp_expect_sort_as,$@.drimseq_es_A_v_B.tsv,test-diffsplice.drimseq_es_A_v_B.tsv)

##
# events at the same coordinates on both strands.  The strand is part of the event
# id; without it the two strands produced one feature_id twice and DRIMSeq refused
# the file.  No --test, so this runs without DRIMSeq installed
##
test-diffsplice-strands: mkdirs
	flair diffsplice --isoform_bed $(STRANDS_DIFFSPLICE_BED) --counts_matrix $(STRANDS_DIFFSPLICE_MATRIX) -o $O/$@
	$(call cmp_expect,$@.alt3.events.quant.tsv)
	$(call cmp_expect,$@.alt5.events.quant.tsv)
	$(call cmp_expect,$@.es.events.quant.tsv)
	$(call cmp_expect,$@.ir.events.quant.tsv)

test-diffsplice-help: mkdirs
	flair diffsplice --help > $O/$@.out
	$(call cmp_expect_help,$@.out)


bed-to-gtf-tests: test-bed-to-gtf test-bed-to-gtf-compress test-bed-to-gtf-extracols test-bed-to-gtf-no-genes

# a bed not from flair, of names transcript_gene, which are not split: each name is
# its gene id too.  Its thick range covers each transcript, which is a coding region
# over every exon; a bed from a FLAIR before 3.1, which wrote that for a transcript
# with no coding region, is converted with --no_cds
test-bed-to-gtf : mkdirs
	bed_to_gtf --force $(COLLAPSE_ANNOT_ISOFORMS_BED) > $O/$@.gtf
	$(call cmp_expect,$@.gtf)

test-bed-to-gtf-compress : mkdirs ${TI}/test-collapse-annot.isoforms.bed.gz
	bed_to_gtf --force ${TI}/test-collapse-annot.isoforms.bed.gz > $O/$@.gtf
	$(call cmp_expect_named,$@.gtf,test-bed-to-gtf.gtf)

test-bed-to-gtf-extracols : mkdirs
	bed_to_gtf --is_flair_bed $(TOME_ISOFORM_BED) > $O/$@.gtf
	$(call cmp_expect,$@.gtf)

# without --force, a bed not from flair, of no gene ids, is refused
test-bed-to-gtf-no-genes : mkdirs
	! bed_to_gtf $(COLLAPSE_ANNOT_ISOFORMS_BED) > $O/$@.gtf 2> $O/$@.err
	grep -q 'has no gene ids' $O/$@.err


##
# the other installed programs, each run on the test data
##
scripts-tests: test-gtf-to-bed test-mark-intron-retention test-fasta-seq-lengths \
	test-diffsplice-fishers-exact

test-gtf-to-bed: mkdirs
	gtf_to_bed --include_gene $(BASIC_ANNOTATION) $O/$@.bed
	$(call cmp_expect,$@.bed)

# the isoforms of the basic transcriptome, with retained introns marked
test-mark-intron-retention: mkdirs
	mark_intron_retention $(TOME_ISOFORM_BED) $O/$@.marked.bed $O/$@.introns.tsv
	$(call cmp_expect,$@.marked.bed)
	$(call cmp_expect,$@.introns.tsv)

test-fasta-seq-lengths: mkdirs
	fasta_seq_lengths $E/test-transcriptome.isoforms.fa $O/$@.lengths.tsv $O/$@.histogram.tsv
	$(call cmp_expect,$@.lengths.tsv)
	$(call cmp_expect,$@.histogram.tsv)

# two samples of the committed diffsplice counts, so no R is needed
test-diffsplice-fishers-exact: mkdirs
	diffsplice_fishers_exact $E/test-diffsplice.es.events.quant.tsv A1_A_b0 B1_B_b0 $O/$@.tsv
	$(call cmp_expect,$@.tsv)


##
# DIFF ISO USAGE
##
diff-iso-usage-tests: test-diff-iso-usage test-diff-iso-usage-compress test-diff-iso-usage-help

test-diff-iso-usage: mkdirs
	diff_iso_usage $(BASIC_COUNTS_MATRIX) MCF7cdnar1r3_MCF7_batch1 MCF7cdnastrandr3r2_MCF7_batch1 $O/$@.tsv
	$(call cmp_expect,$@.tsv)

test-diff-iso-usage-compress: mkdirs ${TI}/basic.counts_matrix.diffexp.tsv.gz
	diff_iso_usage ${TI}/basic.counts_matrix.diffexp.tsv.gz MCF7cdnar1r3_MCF7_batch1 MCF7cdnastrandr3r2_MCF7_batch1 $O/$@.tsv
	$(call cmp_expect_named,$@.tsv,test-diff-iso-usage.tsv)

test-diff-iso-usage-help: mkdirs
	diff_iso_usage --help > $O/$@.out
	$(call cmp_expect_help,$@.out)

# ##
# # plot isoform usage
# ##
plot-usage-tests: test-plot-usage test-plot-usage-compress test-plot-usage-help

test-plot-usage: mkdirs
	plot_isoform_usage $(TOME_ISOFORM_BED) $E/test-quantify.isoforms.counts.tsv FLG00000002 --min_reads 1 -o $O/$@
	file $O/test-plot-usage_FLG00000002_usage.pdf | grep 'PDF document' >/dev/null
	file $O/test-plot-usage_FLG00000002_isoforms.pdf | grep 'PDF document' >/dev/null

test-plot-usage-compress: mkdirs ${TI}/test-transcriptome.isoforms.bed.gz ${TI}/test-quantify.isoforms.counts.tsv.gz
	plot_isoform_usage ${TI}/test-transcriptome.isoforms.bed.gz ${TI}/test-quantify.isoforms.counts.tsv.gz \
	   FLG00000002 --min_reads 1 -o $O/$@
	file $O/$@_FLG00000002_usage.pdf | grep 'PDF document' >/dev/null
	file $O/$@_FLG00000002_isoforms.pdf | grep 'PDF document' >/dev/null

test-plot-usage-help: mkdirs
	plot_isoform_usage --help > $O/$@.out
	$(call cmp_expect_help,$@.out)

###
# flair transcriptome
###

transcriptome-tests: test-transcriptome test-transcriptome-directrna test-transcriptome-region test-transcriptome-readthrough-region test-transcriptome-lowannot test-transcriptome-noannot test-transcriptome-noannot-trust-junctions test-transcriptome-normends test-transcriptome-comprehensive test-transcriptome-comprehensive-basic \
                    test-transcriptome-seg1 test-transcriptome-mutant-splice test-transcriptome-readthrough test-transcriptome-readthrough-gene test-transcriptome-nosubset test-transcriptome-filter-1.2X test-transcriptome-nopolya test-transcriptome-nopolya-trust-strand test-transcriptome-keep-intprimed test-transcriptome-keep-intprimed-no-trust-strand test-transcriptome-polyt test-transcriptome-shared-junc test-transcriptome-no-firstpass test-transcriptome-onechrom test-transcriptome-mult-ends \
					test-transcriptome-seg1-badannot test-transcriptome-seg1-badannot-multends test-transcriptome-seg1-badannot-comprehensive test-transcriptome-help \
					test-transcriptome-keep-intermediate-no-temp-dir test-transcriptome-trust-ends-removed \
					test-transcriptome-no-stringent-removed test-transcriptome-no-check-splice-removed \
					test-transcriptome-frac-support-0 test-transcriptome-frac-support-0-basic test-transcriptome-compress

update-transcriptome-tests:
	$(call update_expected,test-transcriptome)


test-transcriptome: mkdirs ${BASIC_ALIGN_BAM}
	rm -rf $O/$@.tmp
	flair transcriptome -b ${BASIC_ALIGN_BAM} -g ${GENOME_FA} --junction_star ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --keep_intermediate --temp_dir $O/$@.tmp --generate_map
	$(diff_bed_gtf)
	$(call cmp_expect,$@.isoforms.fa)
	$(call cmp_expect,$@.isoforms.counts.tsv)
	$(call cmp_expect,$@.aaseq.tsv)
	$(call cmp_expect,$@.isoforms.read.map.tsv)
	$(call cmp_expect,$@.regions.bed)
	$(call cmp_expect,$@.firstpass.bed)
	$(call cmp_expect,$@.firstpass.reallyunfiltered.bed)
	$(call cmp_expect,$@.firstpass.unfiltered.bed)

# --directRNA: trusts read strand and does not tolerate large indels in terminal exons
test-transcriptome-directrna: mkdirs ${BASIC_ALIGN_BAM}
	flair transcriptome -b ${BASIC_ALIGN_BAM} -g ${GENOME_FA} --junction_star ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --directRNA
	$(diff_bed_gtf)
	$(call cmp_expect_named,$@.aaseq.tsv,test-transcriptome.aaseq.tsv)
	$(call cmp_expect,$@.isoforms.counts.tsv)
	$(call cmp_expect_named,$@.isoforms.fa,test-transcriptome.isoforms.fa)
	$(call cmp_expect_named,$@.regions.bed,test-transcriptome.regions.bed)

# --normalize_ends: spliced isoforms report the ends shared by the isoforms of
# their gene with similar terminal splice sites, rather than their reads' ends
test-transcriptome-normends: mkdirs ${BASIC_ALIGN_BAM}
	flair transcriptome -b ${BASIC_ALIGN_BAM} -g ${GENOME_FA} --junction_star ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --generate_map --normalize_ends
	$(diff_bed_gtf)
	$(call cmp_expect,$@.isoforms.read.map.tsv)
	$(call cmp_expect_named,$@.aaseq.tsv,test-transcriptome.aaseq.tsv)
	$(call cmp_expect_named,$@.isoforms.counts.tsv,test-transcriptome.isoforms.counts.tsv)
	$(call cmp_expect,$@.isoforms.fa)
	$(call cmp_expect_named,$@.regions.bed,test-transcriptome.regions.bed)

# test of --filter comprehensive, which keeps the subset isoforms that the
# default nosubset filter drops from test-transcriptome-mutant-splice
test-transcriptome-comprehensive: mkdirs ${KRAS_MUT_SPLICE_BAM}
	flair transcriptome -b ${KRAS_MUT_SPLICE_BAM} -g ${GENOME_FA} \
	   --junction_bed=${KRAS_MUT_SPLICE_SJ_BED} \
	   -t 4 -f ${BASIC_ANNOTATION} -o $O/$@ --sample_name $@ --sjc_support 2 --filter comprehensive
	$(diff_bed_gtf)
	$(call cmp_expect,$@.aaseq.tsv)
	$(call cmp_expect,$@.isoforms.counts.tsv)
	$(call cmp_expect,$@.isoforms.fa)
	$(call cmp_expect,$@.regions.bed)

# --filter comprehensive on the basic data, which has single-exon isoforms
# contained in spliced exons
test-transcriptome-comprehensive-basic: mkdirs ${BASIC_ALIGN_BAM}
	flair transcriptome -b ${BASIC_ALIGN_BAM} -g ${GENOME_FA} --junction_star ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --filter comprehensive
	$(diff_bed_gtf)
	$(call cmp_expect_named,$@.aaseq.tsv,test-transcriptome.aaseq.tsv)
	$(call cmp_expect_named,$@.isoforms.counts.tsv,test-transcriptome.isoforms.counts.tsv)
	$(call cmp_expect_named,$@.isoforms.fa,test-transcriptome.isoforms.fa)
	$(call cmp_expect_named,$@.regions.bed,test-transcriptome.regions.bed)

test-transcriptome-mult-ends: mkdirs ${BASIC_ALIGN_BAM}
	rm -rf $O/$@.tmp
	flair transcriptome  -b ${BASIC_ALIGN_BAM} -g ${GENOME_FA} --junction_star ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --max_ends 5 --keep_intermediate --temp_dir $O/$@.tmp
	$(diff_bed_gtf)
	$(call cmp_expect_named,$@.aaseq.tsv,test-transcriptome.aaseq.tsv)
	$(call cmp_expect_named,$@.firstpass.bed,test-transcriptome.firstpass.bed)
	$(call cmp_expect_named,$@.firstpass.reallyunfiltered.bed,test-transcriptome.firstpass.reallyunfiltered.bed)
	$(call cmp_expect_named,$@.firstpass.unfiltered.bed,test-transcriptome.firstpass.unfiltered.bed)
	$(call cmp_expect,$@.isoforms.counts.tsv)
	$(call cmp_expect,$@.isoforms.fa)
	$(call cmp_expect_named,$@.regions.bed,test-transcriptome.regions.bed)

# --region_bed builds transcripts in one region, grown to the partition holding every
# read and annotation that reaches it, so they are the whole-genome run's transcripts
# for that gene; the sample name is that run's, as the output is compared to its
# expected files, and the second region of the BED has no reads, which is a warning.
# --region repeats the first region of the BED, in one-based coordinates, so the
# regions of the two options merge back to the regions of the BED
test-transcriptome-region: mkdirs ${BASIC_ALIGN_BAM}
	flair transcriptome -b ${BASIC_ALIGN_BAM} -g ${GENOME_FA} --junction_star ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name test-transcriptome --region_bed $I/transcriptome.region.bed --region=chr12:25205001-25210000
	diff <(awk '$$1 == "chr12"' $E/test-transcriptome.isoforms.bed) $O/$@.isoforms.bed
	diff <(awk '$$1 == "chr12"' $E/test-transcriptome.isoforms.gtf) $O/$@.isoforms.gtf
	$(call cmp_expect,$@.regions.bed)
	$(call cmp_expect,$@.aaseq.tsv)
	$(call cmp_expect,$@.isoforms.counts.tsv)
	$(call cmp_expect,$@.isoforms.fa)

# a 1kb region inside DDX5 is expanded over the readthrough reads into POLG2, the gene
# next to it, so the region reproduces a whole-genome run, POLG2's transcripts included
test-transcriptome-readthrough-region: mkdirs ${DDX5_READTHROUGH_BAM}
	flair transcriptome -b ${DDX5_READTHROUGH_BAM} -g ${GENOME_FA} -f ${DDX5_READTHROUGH_ANNOT} -o $O/$@ --sample_name test-transcriptome-readthrough --junction_bed=${DDX5_READTHROUGH_SJ} --sjc_support 2 --region chr17:64500001-64501000
	diff $E/test-transcriptome-readthrough.isoforms.bed $O/$@.isoforms.bed
	diff $E/test-transcriptome-readthrough.isoforms.gtf $O/$@.isoforms.gtf
	$(call cmp_expect,$@.regions.bed)
	$(call cmp_expect_named,$@.aaseq.tsv,test-transcriptome-readthrough.aaseq.tsv)
	$(call cmp_expect_named,$@.isoforms.counts.tsv,test-transcriptome-readthrough.isoforms.counts.tsv)
	$(call cmp_expect_named,$@.isoforms.fa,test-transcriptome-readthrough.isoforms.fa)

test-transcriptome-lowannot: mkdirs ${BASIC_ALIGN_BAM}
	intron-prospector --genome-fasta=${GENOME_FA} --intron-bed6=$O/$@.ip_junctions.bed -C 0.0 ${BASIC_ALIGN_BAM}
	flair transcriptome  -b ${BASIC_ALIGN_BAM} -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) --junction_bed $O/$@.ip_junctions.bed --junction_support 2 -o $O/$@ --sample_name $@ --no_align_to_annot --sjc_support 3
	$(diff_bed_gtf)
	$(call cmp_expect_named,$@.aaseq.tsv,test-transcriptome.aaseq.tsv)
	$(call cmp_expect_named,$@.isoforms.counts.tsv,test-transcriptome.isoforms.counts.tsv)
	$(call cmp_expect_named,$@.isoforms.fa,test-transcriptome.isoforms.fa)
	$(call cmp_expect_named,$@.regions.bed,test-transcriptome.regions.bed)

test-transcriptome-noannot: mkdirs ${BASIC_ALIGN_BAM}
	intron-prospector --genome-fasta=${GENOME_FA} --intron-bed6=$O/$@.ip_junctions.bed -C 0.0 ${BASIC_ALIGN_BAM}
	flair transcriptome  -b ${BASIC_ALIGN_BAM} -g ${GENOME_FA} -t 4 -o $O/$@ --sample_name $@ --junction_bed $O/$@.ip_junctions.bed --junction_support 2 --sjc_support 3
	$(diff_bed_gtf)
	$(call cmp_expect,$@.aaseq.tsv)
	$(call cmp_expect_named,$@.isoforms.counts.tsv,test-transcriptome.isoforms.counts.tsv)
	$(call cmp_expect_named,$@.isoforms.fa,test-transcriptome.isoforms.fa)
	$(call cmp_expect,$@.regions.bed)

# no annotation, so junction strands come only from intron-prospector; --trust_junctions
# takes every junction's strand as given, without checking for a GT-AG motif
test-transcriptome-noannot-trust-junctions: mkdirs ${BASIC_ALIGN_BAM}
	intron-prospector --genome-fasta=${GENOME_FA} --intron-bed6=$O/$@.ip_junctions.bed -C 0.0 ${BASIC_ALIGN_BAM}
	flair transcriptome  -b ${BASIC_ALIGN_BAM} -g ${GENOME_FA} -t 4 -o $O/$@ --sample_name $@ --junction_bed $O/$@.ip_junctions.bed --junction_support 2 --sjc_support 3 \
	   --trust_junctions
	$(diff_bed_gtf)
	$(call cmp_expect,$@.aaseq.tsv)
	$(call cmp_expect_named,$@.isoforms.counts.tsv,test-transcriptome.isoforms.counts.tsv)
	$(call cmp_expect_named,$@.isoforms.fa,test-transcriptome.isoforms.fa)
	$(call cmp_expect,$@.regions.bed)

test-transcriptome-seg1: mkdirs  ${SEG1_ALIGN_BAM}
	flair transcriptome -b ${SEG1_ALIGN_BAM} -g ${GENOME_FA} -t 4 -f $(SEG1_ANNOTATION) -o $O/$@ --sample_name $@
	$(diff_bed_gtf)
	$(call cmp_expect,$@.aaseq.tsv)
	$(call cmp_expect,$@.isoforms.counts.tsv)
	$(call cmp_expect,$@.isoforms.fa)
	$(call cmp_expect,$@.regions.bed)

test-transcriptome-seg1-badannot: mkdirs  ${SEG1_ALIGN_BAM}
	flair transcriptome -b ${SEG1_ALIGN_BAM} -g ${GENOME_FA} -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@
	$(diff_bed_gtf)
	$(call cmp_expect,$@.aaseq.tsv)
	$(call cmp_expect,$@.isoforms.counts.tsv)
	$(call cmp_expect,$@.isoforms.fa)
	$(call cmp_expect,$@.regions.bed)

test-transcriptome-seg1-badannot-multends: mkdirs  ${SEG1_ALIGN_BAM}
	flair transcriptome -b ${SEG1_ALIGN_BAM} -g ${GENOME_FA} -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --max_ends 4 --single_exon_support 2 --filter ginormous
	$(diff_bed_gtf)
	$(call cmp_expect,$@.aaseq.tsv)
	$(call cmp_expect,$@.isoforms.counts.tsv)
	$(call cmp_expect,$@.isoforms.fa)
	$(call cmp_expect,$@.regions.bed)

# as seg1-badannot-multends, but comprehensive compares the overlapping single-exon
# isoforms by expression instead of keeping them all
test-transcriptome-seg1-badannot-comprehensive: mkdirs  ${SEG1_ALIGN_BAM}
	flair transcriptome -b ${SEG1_ALIGN_BAM} -g ${GENOME_FA} -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --max_ends 4 --single_exon_support 2 --filter comprehensive
	$(diff_bed_gtf)
	$(call cmp_expect,$@.aaseq.tsv)
	$(call cmp_expect,$@.isoforms.counts.tsv)
	$(call cmp_expect,$@.isoforms.fa)
	$(call cmp_expect,$@.regions.bed)

# test of an artifical mutated splice site.
test-transcriptome-mutant-splice: mkdirs ${KRAS_MUT_SPLICE_BAM}
	flair transcriptome -b ${KRAS_MUT_SPLICE_BAM} -g ${GENOME_FA} \
	   --junction_bed=${KRAS_MUT_SPLICE_SJ_BED} \
	   -t 4 -f ${BASIC_ANNOTATION} -o $O/$@ --sample_name $@ --sjc_support 2
	$(diff_bed_gtf)
	$(call cmp_expect,$@.aaseq.tsv)
	$(call cmp_expect,$@.isoforms.counts.tsv)
	$(call cmp_expect,$@.isoforms.fa)
	$(call cmp_expect,$@.regions.bed)

#test of an artificial read through transcript
test-transcriptome-readthrough: mkdirs ${DDX5_READTHROUGH_BAM}
	flair transcriptome -b ${DDX5_READTHROUGH_BAM} -g ${GENOME_FA} -f ${DDX5_READTHROUGH_ANNOT} -o $O/$@ --sample_name $@ --junction_bed=${DDX5_READTHROUGH_SJ} --sjc_support 2
	$(diff_bed_gtf)
	$(call cmp_expect,$@.aaseq.tsv)
	$(call cmp_expect,$@.isoforms.counts.tsv)
	$(call cmp_expect,$@.isoforms.fa)
	$(call cmp_expect,$@.regions.bed)


# the readthrough transcript annotated as its own gene: the readthrough isoform
# matches it, instead of being a novel isoform assigned to both DDX5 and POLG2
test-transcriptome-readthrough-gene: mkdirs ${DDX5_READTHROUGH_BAM} ${DDX5_READTHROUGH_GENE_ANNOT}
	flair transcriptome -b ${DDX5_READTHROUGH_BAM} -g ${GENOME_FA} -f ${DDX5_READTHROUGH_GENE_ANNOT} -o $O/$@ --sample_name $@ --junction_bed=${DDX5_READTHROUGH_SJ} --sjc_support 2
	$(diff_bed_gtf)
	$(call cmp_expect,$@.aaseq.tsv)
	$(call cmp_expect,$@.isoforms.counts.tsv)
	$(call cmp_expect,$@.isoforms.fa)
	$(call cmp_expect,$@.regions.bed)

${DDX5_READTHROUGH_GENE_ANNOT}: ${DDX5_READTHROUGH_ANNOT} ${DDX5_READTHROUGH_GENE_GTF}
	@mkdir -p $O
	cat ${DDX5_READTHROUGH_ANNOT} ${DDX5_READTHROUGH_GENE_GTF} > $@.tmp
	mv -f $@.tmp $@

# test of --filter nosubset (the default), which removes all subset isoforms, on
# reads chosen so subset isoforms have more reads than their supersets
test-transcriptome-nosubset: mkdirs ${DDX5_BYSUPPORT_BAM}
	flair transcriptome -b ${DDX5_BYSUPPORT_BAM} -g ${GENOME_FA} -f ${DDX5_READTHROUGH_ANNOT} -o $O/$@ --sample_name $@ --junction_bed=${DDX5_READTHROUGH_SJ} --sjc_support 2 \
	   --filter nosubset
	$(diff_bed_gtf)
	$(call cmp_expect,$@.aaseq.tsv)
	$(call cmp_expect,$@.isoforms.counts.tsv)
	$(call cmp_expect,$@.isoforms.fa)
	$(call cmp_expect,$@.regions.bed)

# test of --filter 1.2X on the same reads, which keeps a subset isoform with more
# than 1.2x the reads of the isoforms it is a subset of
test-transcriptome-filter-1.2X: mkdirs ${DDX5_BYSUPPORT_BAM}
	flair transcriptome -b ${DDX5_BYSUPPORT_BAM} -g ${GENOME_FA} -f ${DDX5_READTHROUGH_ANNOT} -o $O/$@ --sample_name $@ --junction_bed=${DDX5_READTHROUGH_SJ} --sjc_support 2 \
	   --filter 1.2X
	$(diff_bed_gtf)
	$(call cmp_expect,$@.aaseq.tsv)
	$(call cmp_expect,$@.isoforms.counts.tsv)
	$(call cmp_expect,$@.isoforms.fa)
	$(call cmp_expect,$@.regions.bed)

${DDX5_BYSUPPORT_BAM}: ${DDX5_READTHROUGH_BAM} ${DDX5_BYSUPPORT_DROP_READS}
	samtools view -h ${DDX5_READTHROUGH_BAM} | grep -v -F -w -f ${DDX5_BYSUPPORT_DROP_READS} | samtools view -b -o $@.tmp.bam -
	mv -f $@.tmp.bam $@
	samtools index $@

# basic reads with poly(A) trimmed, keeping single-exon isoforms.  Without poly(A),
# single-exon reads can only be stranded with --trust_strand, so compare these two
test-transcriptome-nopolya: mkdirs ${NOPOLYA_ALIGN_BAM}
	flair transcriptome -b ${NOPOLYA_ALIGN_BAM} -g ${GENOME_FA} --junction_star ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ \
	   --filter ginormous --single_exon_support 2
	$(diff_bed_gtf)
	$(call cmp_expect_named,$@.aaseq.tsv,test-transcriptome.aaseq.tsv)
	$(call cmp_expect,$@.isoforms.counts.tsv)
	$(call cmp_expect,$@.isoforms.fa)
	$(call cmp_expect,$@.regions.bed)

test-transcriptome-nopolya-trust-strand: mkdirs ${NOPOLYA_ALIGN_BAM}
	flair transcriptome -b ${NOPOLYA_ALIGN_BAM} -g ${GENOME_FA} --junction_star ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ \
	   --filter ginormous --single_exon_support 2 --trust_strand
	$(diff_bed_gtf)
	$(call cmp_expect_named,$@.aaseq.tsv,test-transcriptome.aaseq.tsv)
	$(call cmp_expect,$@.isoforms.counts.tsv)
	$(call cmp_expect,$@.isoforms.fa)
	$(call cmp_expect,$@.regions.bed)

# poly(T) heads on all reads, no annotation: the single-exon KRAS isoform gets re-assigned to the + strand
# all spliced reads strand as normal due to splice junction strand correction, only single exon reads get re-stranded
test-transcriptome-polyt: mkdirs ${POLYT_ALIGN_BAM}
	intron-prospector --genome-fasta=${GENOME_FA} --intron-bed6=$O/$@.ip_junctions.bed -C 0.0 ${POLYT_ALIGN_BAM}
	flair transcriptome -b ${POLYT_ALIGN_BAM} -g ${GENOME_FA} -t 4 -o $O/$@ --sample_name $@ --junction_bed $O/$@.ip_junctions.bed --junction_support 2 --sjc_support 3 \
	   --filter ginormous --single_exon_support 2
	$(diff_bed_gtf)
	$(call cmp_expect,$@.aaseq.tsv)
	$(call cmp_expect,$@.isoforms.counts.tsv)
	$(call cmp_expect,$@.isoforms.fa)
	$(call cmp_expect,$@.regions.bed)

# novel DDX5 isoforms use the intron shared with the synthetic SHAREDJ1 gene; SHAREDJ1
# has fewer of their junctions, so gene assignment drops it and the results match
# test-transcriptome
# deliberately uses the deprecated --junction_tab, so that spelling stays covered
test-transcriptome-shared-junc: mkdirs ${BASIC_ALIGN_BAM} ${SHARED_JUNC_ANNOTATION}
	flair transcriptome -b ${BASIC_ALIGN_BAM} -g ${GENOME_FA} --junction_tab ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f ${SHARED_JUNC_ANNOTATION} -o $O/$@ --sample_name $@
	$(diff_bed_gtf)
	$(call cmp_expect_named,$@.aaseq.tsv,test-transcriptome.aaseq.tsv)
	$(call cmp_expect_named,$@.isoforms.counts.tsv,test-transcriptome.isoforms.counts.tsv)
	$(call cmp_expect_named,$@.isoforms.fa,test-transcriptome.isoforms.fa)
	$(call cmp_expect_named,$@.regions.bed,test-transcriptome.regions.bed)

# the chr20 partition has one read, below --sjc_support 2, so it has no firstpass
# isoforms; results should be the chr12 and chr17 isoforms of test-transcriptome
test-transcriptome-no-firstpass: mkdirs ${ONECHR20_ALIGN_BAM}
	flair transcriptome -b ${ONECHR20_ALIGN_BAM} -g ${GENOME_FA} --junction_star ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --generate_map --sjc_support 2
	$(diff_bed_gtf)
	$(call cmp_expect,$@.isoforms.read.map.tsv)
	$(call cmp_expect,$@.aaseq.tsv)
	$(call cmp_expect,$@.isoforms.counts.tsv)
	$(call cmp_expect,$@.isoforms.fa)
	$(call cmp_expect,$@.regions.bed)

# test-transcriptome-noannot with all loci on one chromosome, so KRAS and DDX5 (both
# minus strand) are separate novel genes on the same strand of one partition
test-transcriptome-onechrom: mkdirs ${ONECHROM_ALIGN_BAM}
	intron-prospector --genome-fasta=${ONECHROM_GENOME_FA} --intron-bed6=$O/$@.ip_junctions.bed -C 0.0 ${ONECHROM_ALIGN_BAM}
	flair transcriptome -b ${ONECHROM_ALIGN_BAM} -g ${ONECHROM_GENOME_FA} -t 4 -o $O/$@ --sample_name $@ --junction_bed $O/$@.ip_junctions.bed --junction_support 2 --sjc_support 3
	$(diff_bed_gtf)
	$(call cmp_expect,$@.aaseq.tsv)
	$(call cmp_expect_named,$@.isoforms.counts.tsv,test-transcriptome.isoforms.counts.tsv)
	$(call cmp_expect_named,$@.isoforms.fa,test-transcriptome.isoforms.fa)
	$(call cmp_expect,$@.regions.bed)

${ONECHROM_ALIGN_BAM}: ${BASIC_READS_FA} ${ONECHROM_GENOME_FA}
	$(call minimap2_genome_align,${BASIC_READS_FA},${ONECHROM_GENOME_FA},-s 80 -G 200k)

${KRAS_MUT_SPLICE_BAM}: ${KRAS_MUT_SPLICE_FQ}  ${GENOME_FA}
	$(call minimap2_genome_align,${KRAS_MUT_SPLICE_FQ},${GENOME_FA},-s 40 -G 350k)

#DDX5_readthrough.badread_40x.fastq
${DDX5_READTHROUGH_BAM}: ${DDX5_READTHROUGH_FQ}  ${GENOME_FA}
	$(call minimap2_genome_align,${DDX5_READTHROUGH_FQ},${GENOME_FA},-s 40 -G 350k)

# --keep_intermediate without --temp_dir keeps the run's own directory, in $TMPDIR,
# and reports where it is
test-transcriptome-keep-intermediate-no-temp-dir: mkdirs ${BASIC_ALIGN_BAM}
	TMPDIR=$O/$@.tmpdir; export TMPDIR; mkdir -p $$TMPDIR; \
	  flair transcriptome -b ${BASIC_ALIGN_BAM} -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --keep_intermediate > $O/$@.out 2>&1
	grep -q -- 'intermediate files kept in' $O/$@.out
	test -n "$$(ls $O/$@.tmpdir)" 
	$(call cmp_expect_named,$@.aaseq.tsv,test-transcriptome.aaseq.tsv)
	$(call cmp_expect_named,$@.firstpass.bed,test-transcriptome.firstpass.bed)
	$(call cmp_expect_named,$@.firstpass.reallyunfiltered.bed,test-transcriptome.firstpass.reallyunfiltered.bed)
	$(call cmp_expect_named,$@.firstpass.unfiltered.bed,test-transcriptome.firstpass.unfiltered.bed)
	$(call cmp_expect,$@.isoforms.bed)
	$(call cmp_expect_named,$@.isoforms.counts.tsv,test-transcriptome.isoforms.counts.tsv)
	$(call cmp_expect_named,$@.isoforms.fa,test-transcriptome.isoforms.fa)
	$(call cmp_expect,$@.isoforms.gtf)
	$(call cmp_expect_named,$@.regions.bed,test-transcriptome.regions.bed)

# --keep_internally_primed: single-exon reads ending in or before genomic A runs are
# kept rather than removed as oligo-dT mispriming, for reads not primed on poly(A)
test-transcriptome-keep-intprimed: mkdirs ${BASIC_ALIGN_BAM}
	flair transcriptome -b ${BASIC_ALIGN_BAM} -g ${GENOME_FA} --junction_star ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --trust_strand --keep_internally_primed
	$(diff_bed_gtf)
	$(call cmp_expect_named,$@.aaseq.tsv,test-transcriptome.aaseq.tsv)
	$(call cmp_expect_named,$@.isoforms.counts.tsv,test-transcriptome.isoforms.counts.tsv)
	$(call cmp_expect_named,$@.isoforms.fa,test-transcriptome.isoforms.fa)
	$(call cmp_expect_named,$@.regions.bed,test-transcriptome.regions.bed)

# --keep_internally_primed without --trust_strand: single-exon groups are still
# stranded by their poly(A) tails
test-transcriptome-keep-intprimed-no-trust-strand: mkdirs ${BASIC_ALIGN_BAM}
	flair transcriptome -b ${BASIC_ALIGN_BAM} -g ${GENOME_FA} --junction_star ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --keep_internally_primed
	$(diff_bed_gtf)
	$(call cmp_expect_named,$@.aaseq.tsv,test-transcriptome.aaseq.tsv)
	$(call cmp_expect_named,$@.isoforms.counts.tsv,test-transcriptome.isoforms.counts.tsv)
	$(call cmp_expect_named,$@.isoforms.fa,test-transcriptome.isoforms.fa)
	$(call cmp_expect_named,$@.regions.bed,test-transcriptome.regions.bed)

# test of --frac_support 0: isoforms need only sjc_support and single_exon_support
# reads, though single-exon end variants still need SINGLE_EXON_END_MIN_FRAC of
# their cluster's reads
test-transcriptome-frac-support-0: mkdirs ${NOPOLYA_ALIGN_BAM}
	flair transcriptome -b ${NOPOLYA_ALIGN_BAM} -g ${GENOME_FA} --junction_star ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ \
	   --filter ginormous --single_exon_support 2 --trust_strand --frac_support 0
	$(diff_bed_gtf)
	$(call cmp_expect,$@.aaseq.tsv)
	$(call cmp_expect,$@.isoforms.counts.tsv)
	$(call cmp_expect,$@.isoforms.fa)
	$(call cmp_expect,$@.regions.bed)

# the basic test with --frac_support 0: spliced isoforms need only sjc_support reads,
# 3 here, as the default 1 would keep every corrected junction chain
test-transcriptome-frac-support-0-basic: mkdirs ${BASIC_ALIGN_BAM}
	flair transcriptome -b ${BASIC_ALIGN_BAM} -g ${GENOME_FA} --junction_star ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ \
	   --frac_support 0 --sjc_support 3
	$(diff_bed_gtf)
	$(call cmp_expect,$@.aaseq.tsv)
	$(call cmp_expect,$@.isoforms.counts.tsv)
	$(call cmp_expect,$@.isoforms.fa)
	$(call cmp_expect_named,$@.regions.bed,test-transcriptome.regions.bed)


test-transcriptome-trust-ends-removed: mkdirs ${BASIC_ALIGN_BAM}
	! flair transcriptome -b ${BASIC_ALIGN_BAM} -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --trust_ends > $O/$@.out 2>&1
	grep -q -- '--trust_ends has been removed' $O/$@.out

test-transcriptome-no-stringent-removed: mkdirs ${BASIC_ALIGN_BAM}
	! flair transcriptome -b ${BASIC_ALIGN_BAM} -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --no_stringent > $O/$@.out 2>&1
	grep -q -- '--no_stringent has been removed' $O/$@.out

test-transcriptome-no-check-splice-removed: mkdirs ${BASIC_ALIGN_BAM}
	! flair transcriptome -b ${BASIC_ALIGN_BAM} -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --no_check_splice > $O/$@.out 2>&1
	grep -q -- '--no_check_splice has been removed' $O/$@.out

test-transcriptome-help: mkdirs
	flair transcriptome --help > $O/$@.out
	$(call cmp_expect_help,$@.out)


###
# flair variants
# FIXME need to add better checks for AAseq results, however the inconsistent sorting of the variants leads to inconsistent AAseq IDs which are hard to check
###

BASIC_LONGSHOT = $O/basic.reads.longshot.vcf

${BASIC_LONGSHOT}: ${GENOME_FA} ${BASIC_ALIGN_BAM}
	@mkdir -p $O
	longshot --force_overwrite --bam ${BASIC_ALIGN_BAM} --ref ${GENOME_FA} --out ${BASIC_LONGSHOT} --min_cov 3 --min_alt_count 3 --strand_bias_pvalue_cutoff 0.000001

isoalleles_tests = test-isoalleles test-isoalleles-pair test-isoalleles-compress test-isoalleles-help
variants-tests: test-identifyvars-snvs test-identifyvars-indels test-identifyvars-snvsindels \
	test-identifyvars-thresholds ${isoalleles_tests}

##
# flair isoalleles: the test reads' longshot variants phased by whatshap, which
# haplotags the reads, and flair quantify's read map of them, as a user prepares
# them.
#
# The whatshap output is committed, as bioconda builds whatshap through python 3.13
# and these tests must run on every version `requires-python' spans.  Two
# files are committed: the phased VCF, and whatshap's read-to-haplotype table,
# which bin/bam-add-haplotags applies to the alignments in place of running
# whatshap haplotag.  Regenerate both with `make regen-isoalleles-inputs' when
# the reads, the aligner or longshot change; that needs whatshap, installed as
# docs/source/installing.rst describes.
##
BASIC_PHASED_VCF_IN = $I/basic.reads.phased.vcf
BASIC_HAPLOTAG_LIST = $I/basic.reads.haplotag_list.tsv

BASIC_PHASED_VCF = $O/basic.reads.phased.vcf.gz
BASIC_HAPLOTAGGED_BAM = $O/basic.reads.haplotagged.bam
ISOALLELES_QUANTIFY = $O/basic.reads.isoalleles_quantify
ISOALLELES_READ_MAP = ${ISOALLELES_QUANTIFY}.isoforms.basic.read.map.tsv

# flair isoalleles wants it bgzipped, indexed with pysam, as no tabix is in the
# docker image
${BASIC_PHASED_VCF}: ${BASIC_PHASED_VCF_IN} | mkdirs
	bgzip -c $< > $@.$$$$.tmp && mv -f $@.$$$$.tmp $@
	python3 -c 'import pysam, sys; pysam.tabix_index(sys.argv[1], preset="vcf", force=True)' $@

${BASIC_HAPLOTAGGED_BAM}: ${BASIC_ALIGN_BAM} ${BASIC_HAPLOTAG_LIST} bin/bam-add-haplotags
	bin/bam-add-haplotags ${BASIC_ALIGN_BAM} ${BASIC_HAPLOTAG_LIST} $@.$$$$.tmp && mv -f $@.$$$$.tmp $@
	samtools index $@

# whatshap haplotag reads the VCF through tabix, so it is bgzipped and indexed
# first, in a copy that is thrown away: only the two $I files are kept
regen-isoalleles-inputs: ${GENOME_FA} ${BASIC_ALIGN_BAM} ${BASIC_LONGSHOT}
	whatshap phase --ignore-read-groups --reference ${GENOME_FA} -o ${BASIC_PHASED_VCF_IN} \
	    ${BASIC_LONGSHOT} ${BASIC_ALIGN_BAM}
	bgzip -c ${BASIC_PHASED_VCF_IN} > $O/regen.phased.vcf.gz
	python3 -c 'import pysam, sys; pysam.tabix_index(sys.argv[1], preset="vcf", force=True)' $O/regen.phased.vcf.gz
	whatshap haplotag --ignore-read-groups --reference ${GENOME_FA} \
	    --output-haplotag-list ${BASIC_HAPLOTAG_LIST} -o $O/regen.haplotagged.bam \
	    $O/regen.phased.vcf.gz ${BASIC_ALIGN_BAM}
	rm -f $O/regen.phased.vcf.gz $O/regen.phased.vcf.gz.tbi $O/regen.haplotagged.bam

${ISOALLELES_READ_MAP}: ${BASIC_HAPLOTAGGED_BAM}
	printf 'basic\tA\tb1\t%s\n' ${BASIC_HAPLOTAGGED_BAM} > ${ISOALLELES_QUANTIFY}.manifest.tsv
	flair quantify --manifest ${ISOALLELES_QUANTIFY}.manifest.tsv -g ${GENOME_FA} --isoform_bed $(TOME_ISOFORM_BED) \
	   --generate_map -o ${ISOALLELES_QUANTIFY}

ISOALLELES_INPUTS = --isoform_bed $(TOME_ISOFORM_BED) -g ${GENOME_FA} --phased_vcf ${BASIC_PHASED_VCF}

# one sample: each isoform's reads by haplotype, with each haplotype's sequence and protein
test-isoalleles: mkdirs ${ISOALLELES_READ_MAP} ${BASIC_PHASED_VCF}
	flair isoalleles ${ISOALLELES_INPUTS} --bam ${BASIC_HAPLOTAGGED_BAM} --read_map ${ISOALLELES_READ_MAP} \
	   --sample_name basic -o $O/$@
	$(call cmp_expect,$@.isoalleles.counts.tsv)
	$(call cmp_expect,$@.isoalleles.sequences.tsv)
	$(call cmp_expect,$@.isoalleles.variants.tsv)
	$(call cmp_expect,$@.isoalleles.bed)
	$(call cmp_expect,$@.isoalleles.genes.tsv)
	$(call cmp_expect,$@.isoalleles.excluded.tsv)
	$(call cmp_expect,$@.isoalleles.terminal_sites.tsv)
	$(call cmp_expect,$@.isoalleles.warnings.tsv)
	$(call cmp_expect,$@.isoalleles.aaseq.tsv)
	$(call cmp_expect,$@.isoalleles.fa)
	$(call cmp_expect,$@.isoalleles.somatic_sets.tsv)

# a pair, of the reads as tumor and normal both, the chr17 het given as a somatic
# variant: its reads, on H1, are a somatic set in both
test-isoalleles-pair: mkdirs ${ISOALLELES_READ_MAP} ${BASIC_PHASED_VCF}
	flair isoalleles ${ISOALLELES_INPUTS} --bam ${BASIC_HAPLOTAGGED_BAM} --read_map ${ISOALLELES_READ_MAP} \
	   --sample_name tumor --normal_bam ${BASIC_HAPLOTAGGED_BAM} --normal_read_map ${ISOALLELES_READ_MAP} \
	   --normal_name normal --variant 'chr17:64500016:C>T' -o $O/$@
	$(call cmp_expect,$@.isoalleles.counts.tsv)
	$(call cmp_expect,$@.isoalleles.sequences.tsv)
	$(call cmp_expect,$@.isoalleles.somatic_sets.tsv)
	$(call cmp_expect,$@.isoalleles.variants.tsv)
	$(call cmp_expect,$@.isoalleles.genes.tsv)
	$(call cmp_expect,$@.isoalleles.excluded.tsv)
	$(call cmp_expect,$@.isoalleles.terminal_sites.tsv)
	$(call cmp_expect,$@.isoalleles.warnings.tsv)
	$(call cmp_expect,$@.isoalleles.aaseq.tsv)
	$(call cmp_expect,$@.isoalleles.bed)
	$(call cmp_expect,$@.isoalleles.fa)

# compressed, in two threads, with read maps and labelled bams: the same tables
test-isoalleles-compress: mkdirs ${ISOALLELES_READ_MAP} ${BASIC_PHASED_VCF}
	flair -z isoalleles ${ISOALLELES_INPUTS} --bam ${BASIC_HAPLOTAGGED_BAM} --read_map ${ISOALLELES_READ_MAP} \
	   --sample_name basic -t 2 --generate_map --label_bams -o $O/$@
	$(call cmp_expect_gz,$@.isoalleles.counts.tsv,test-isoalleles.isoalleles.counts.tsv)
	$(call cmp_expect_gz,$@.isoalleles.sequences.tsv,test-isoalleles.isoalleles.sequences.tsv)
	$(call check_fasta_indexable,$@.isoalleles.fa)
	samtools view $O/$@.isoalleles.basic.bam | grep -q 'IA:Z:FLT00000003|H2'

##
# misc/flair_protocol.py, the protocol from one manifest: a prototype, so the runs
# of it are asked for (make test-protocol) rather than in a default suite.  Its
# unit tests, which read the planned steps and run no pipeline, are in test-base.
# The first run of each is --force, so earlier outputs can't hide a failure, and a
# rerun must find every step done
##
PROTOCOL = ../misc/flair_protocol.py

test-protocol: test-protocol-unit test-protocol-patients test-protocol-haplotype

test-protocol-unit:
	FLAIR_PROTOCOL_TESTS=1 ${PYTEST_CMD} test_flair_protocol.py

# a patient of a normal and a tumor, and a tumor-only patient, with splice events:
# junctions from intron-prospector, and only the pair tested
test-protocol-patients: mkdirs ${GENOME_FA} ${QUANTIFY_INPUT_BAMS}
	printf 'sample\tbam\tpatient\ttissue\nP1N\t%s\tP1\tnormal\nP1T\t%s\tP1\ttumor\nP2T\t%s\tP2\ttumor\n' \
	   $O/basic.reads.1.genome_aligned.bam $O/basic.reads.2.genome_aligned.bam $O/basic.reads.4.genome_aligned.bam \
	   > $O/$@.manifest.tsv
	python3 ${PROTOCOL} -m $O/$@.manifest.tsv -g ${GENOME_FA} -f $(BASIC_ANNOTATION) -o $O/$@ -t 2 \
	   --spliceevents --stat_test patients --force > $O/$@.out
	test -s $O/$@.quantify/flair.isoforms.counts.tsv
	test -s $O/$@.spliceevents/flair.spliceevents.counts.tsv
	test -s $O/$@.stats/P1.P1T.isoforms.usage.tsv
	test -s $O/$@.stats/P1.P1T.spliceevents.usage.tsv
	test ! -e $O/$@.stats/P2.P2T.isoforms.usage.tsv
	grep -q 'tumor-only patients, reported but not tested: P2' $O/$@.out
	python3 ${PROTOCOL} -m $O/$@.manifest.tsv -g ${GENOME_FA} -f $(BASIC_ANNOTATION) -o $O/$@ -t 2 \
	   --spliceevents --stat_test patients > $O/$@.rerun.out
	! grep -q '^\[run\]' $O/$@.rerun.out

# --haplotype: the test reads, called by longshot, as a pair and as a tumor-only
# patient: phasing, haplotagging, and isoalleles of each
test-protocol-haplotype: mkdirs ${GENOME_FA} ${BASIC_ALIGN_BAM} ${BASIC_LONGSHOT}
	printf 'sample\tbam\tpatient\ttissue\tvcf\nP1N\t%s\tP1\tnormal\t%s\nP1T\t%s\tP1\ttumor\t%s\nP2T\t%s\tP2\ttumor\t%s\n' \
	   ${BASIC_ALIGN_BAM} ${BASIC_LONGSHOT} ${BASIC_ALIGN_BAM} ${BASIC_LONGSHOT} ${BASIC_ALIGN_BAM} ${BASIC_LONGSHOT} \
	   > $O/$@.manifest.tsv
	python3 ${PROTOCOL} -m $O/$@.manifest.tsv -g ${GENOME_FA} -f $(BASIC_ANNOTATION) -o $O/$@ -t 2 \
	   --stat_test patients --haplotype --force > $O/$@.out
	test -s $O/$@.isoalleles/P1T.isoalleles.counts.tsv
	test -s $O/$@.isoalleles/P2T.isoalleles.counts.tsv
	test -s $O/$@.stats/P1.P1T.isoalleles.usage.tsv
	head -1 $O/$@.isoalleles/P1T.isoalleles.counts.tsv | grep -q 'P1T	P1N'

test-identifyvars-snvs: mkdirs ${GENOME_FA} ${BASIC_ALIGN_BAM}
	python3 ../src/flair/identify_vars.py --genome ${GENOME_FA} --bam ${BASIC_ALIGN_BAM} --identify_snvs yes --identify_indels no --output $O/$@
	$(call cmp_expect,$@.vcf)

test-identifyvars-indels: mkdirs ${GENOME_FA} ${BASIC_ALIGN_BAM}
	python3 ../src/flair/identify_vars.py --genome ${GENOME_FA} --bam ${BASIC_ALIGN_BAM} --identify_snvs no --identify_indels yes --output $O/$@
	$(call cmp_expect,$@.vcf)

test-identifyvars-snvsindels: mkdirs ${GENOME_FA} ${BASIC_ALIGN_BAM}
	python3 ../src/flair/identify_vars.py --genome ${GENOME_FA} --bam ${BASIC_ALIGN_BAM} --identify_snvs yes --identify_indels yes --output $O/$@
	$(call cmp_expect,$@.vcf)

# every call filter at a value of its own, which keeps 23 of the 51 calls of
# test-identifyvars-snvsindels
test-identifyvars-thresholds: mkdirs ${GENOME_FA} ${BASIC_ALIGN_BAM}
	python3 ../src/flair/identify_vars.py --genome ${GENOME_FA} --bam ${BASIC_ALIGN_BAM} --identify_snvs yes --identify_indels yes --output $O/$@ \
	   --min_var_reads 6 --min_vaf_pct 25 --min_sj_dist 15 --min_read_end_dist 20 --density_dist 50 --density_count 2 --num_repeats_for_filtering 4
	$(call cmp_expect,$@.vcf)

###
# flair fusion
###

fusion-tests: test-fusion test-fusion-pary test-fusion-nongenic test-fusion-no-tome test-fusion-3gene test-fusion-sorted-gtf \
	test-fusion-refseq test-fusion-recall test-fusion-help test-fusion-compress test-fusion-max-partners \
	test-fusion-reuse test-fusion-reuse-other-reads

update-fusion-tests:
	$(call update_expected,test-fusion)

test-fusion: mkdirs ${GENOME_FA} ${ERGIC3_KRAS_SIM_FA_GZ} ${FUSION_TEST_BAM}
	flair fusion -g ${GENOME_FA} -f $(BASIC_ANNOTATION) -b ${FUSION_TEST_BAM} --max_loci 4 -o $O/$@ --sample_name $@ --min_support 2
	$(call cmp_expect,$@.fusions.isoforms.bed)
	$(call cmp_expect,$@.fusions.isoforms.fa)
	$(call cmp_expect,$@.fusions.tsv)
	$(call cmp_expect,$@.fusions.isoforms.counts.tsv)
	$(call cmp_expect,$@.fusions.aaseq.tsv)
	test -s $O/$@.fusions.isoforms.read.map.tsv -a -s $O/$@.fusions.isoforms.gtf -a -s $O/$@.fusions.isoforms.bb
	$(call cmp_expect,$@.fusions.isoforms.gtf)
	$(call cmp_expect_read_map,$@.fusions.isoforms.read.map.tsv)
	$(call cmp_expect,$@.fusions.prelim.bed)
# 	diff <(sort <(cut -f 1,2,3,5,6,7,8,9,10,11,12 $E/$@.fusions.isoforms.bed)) <(sort <(cut -f 1,2,3,5,6,7,8,9,10,11,12 $O/$@.fusions.isoforms.bed))
# 	diff <(sort <(grep -v ">" $E/$@.fusions.isoforms.fa)) <(sort <(grep -v ">" $O/$@.fusions.isoforms.fa))

# test-fusion for the most recall: every spliced subset fusion isoform kept, and every
# isoform whatever its part of the fusion's reads
test-fusion-recall: mkdirs ${GENOME_FA} ${ERGIC3_KRAS_SIM_FA_GZ} ${FUSION_TEST_BAM}
	flair fusion -g ${GENOME_FA} -f $(BASIC_ANNOTATION) -b ${FUSION_TEST_BAM} --max_loci 4 -o $O/$@ --sample_name $@ --min_support 2 \
	    --filter comprehensive --frac_support 0
	$(call cmp_expect,$@.fusions.isoforms.bed)
	$(call cmp_expect,$@.fusions.isoforms.fa)
	$(call cmp_expect,$@.fusions.tsv)
	$(call cmp_expect,$@.fusions.isoforms.counts.tsv)
	$(call cmp_expect_named,$@.fusions.aaseq.tsv,test-fusion.fusions.aaseq.tsv)
	$(call cmp_expect,$@.fusions.isoforms.gtf)
	$(call cmp_expect_read_map_named,$@.fusions.isoforms.read.map.tsv,test-fusion.fusions.isoforms.read.map.tsv)
	$(call cmp_expect_named,$@.fusions.prelim.bed,test-fusion.fusions.prelim.bed)

# the sample name is the one test-fusion used, since it is written into the BED
test-fusion-compress: mkdirs ${GENOME_FA} ${ERGIC3_KRAS_SIM_FA_GZ} ${FUSION_TEST_BAM} ${TI}/basic.annotation.gtf.gz
	flair -z fusion -g ${GENOME_FA} -f ${TI}/basic.annotation.gtf.gz -b ${FUSION_TEST_BAM} --max_loci 4 -o $O/$@ --sample_name test-fusion --min_support 2
	$(call cmp_expect_gz,$@.fusions.isoforms.bed,test-fusion.fusions.isoforms.bed)
	$(call cmp_expect_gz,$@.fusions.isoforms.fa,test-fusion.fusions.isoforms.fa)
	$(call check_gz,$@.fusions.isoforms.gtf)
	$(call check_fasta_indexable,$@.fusions.isoforms.fa)
	$(call cmp_expect_named,$@.fusions.prelim.bed,test-fusion.fusions.prelim.bed)

# --max_partners 0 makes every gene promiscuous, so the fusion test-fusion reports
# is dropped: the option reaches the partner count
test-fusion-max-partners: mkdirs ${GENOME_FA} ${ERGIC3_KRAS_SIM_FA_GZ} ${FUSION_TEST_BAM}
	flair fusion -g ${GENOME_FA} -f $(BASIC_ANNOTATION) -b ${FUSION_TEST_BAM} --max_loci 4 -o $O/$@ --sample_name $@ --min_support 2 --max_partners 0
	$(call cmp_expect,$@.fusions.tsv)
	test ! -s $O/$@.fusions.isoforms.bed
	$(call cmp_expect_named,$@.fusions.aaseq.tsv,test-fusion.fusions.aaseq.tsv)
	$(call cmp_expect,$@.fusions.isoforms.counts.tsv)
	$(call cmp_expect,$@.fusions.isoforms.fa)
	$(call cmp_expect_read_map,$@.fusions.isoforms.read.map.tsv)
	$(call cmp_expect,$@.fusions.prelim.bed)

test-fusion-pary: mkdirs ${GENOME_FA} ${ERGIC3_KRAS_SIM_FA_GZ} ${FUSION_TEST_BAM}
	flair fusion -g ${GENOME_FA} -f $I/parY.annotation.gtf -b ${FUSION_TEST_BAM} --max_loci 4 -o $O/$@ --sample_name $@ --min_support 2
	$(call cmp_expect,$@.fusions.isoforms.bed)
	$(call cmp_expect,$@.fusions.isoforms.fa)
	$(call cmp_expect,$@.fusions.tsv)
	$(call cmp_expect_named,$@.fusions.aaseq.tsv,test-fusion.fusions.aaseq.tsv)
	$(call cmp_expect_named,$@.fusions.isoforms.counts.tsv,test-fusion.fusions.isoforms.counts.tsv)
	$(call cmp_expect,$@.fusions.isoforms.gtf)
	$(call cmp_expect_read_map_named,$@.fusions.isoforms.read.map.tsv,test-fusion.fusions.isoforms.read.map.tsv)
	$(call cmp_expect,$@.fusions.prelim.bed)

test-fusion-nongenic: mkdirs ${GENOME_FA} ${ERGIC3_KRAS_SIM_FA_GZ} ${FUSION_TEST_BAM}
	flair fusion -g ${GENOME_FA} -f $I/basic.annotation.ERGIC3_only.gtf -b ${FUSION_TEST_BAM} --max_loci 4 -o $O/$@ --sample_name $@ --min_support 2
	$(call cmp_expect,$@.fusions.isoforms.bed)
	$(call cmp_expect,$@.fusions.isoforms.fa)
	$(call cmp_expect,$@.fusions.tsv)
	$(call cmp_expect_named,$@.fusions.aaseq.tsv,test-fusion.fusions.aaseq.tsv)
	$(call cmp_expect,$@.fusions.isoforms.counts.tsv)
	$(call cmp_expect,$@.fusions.isoforms.gtf)
	$(call cmp_expect_read_map,$@.fusions.isoforms.read.map.tsv)
	$(call cmp_expect,$@.fusions.prelim.bed)

# a rerun under the same output name with other reads: the chimeric BAMs the first
# run left are made again, since they are not of these reads, and the result is that
# of a run with nothing to reuse.  Reusing them on their existence alone reported the
# first run's fusion for the second run's reads.
test-fusion-reuse-other-reads: mkdirs ${GENOME_FA} ${FUSION_TEST_BAM} ${FUSION3_TEST_BAM}
	for f in genome_aligned transcriptome_aligned; do rm -f $O/$@.$$f.chim.bam $O/$@.$$f.chim.bam.bai $O/$@.$$f.chim.bam.provenance.json; done
	flair fusion -g ${GENOME_FA} -f $(BASIC_ANNOTATION) -b ${FUSION_TEST_BAM} --max_loci 3 -o $O/$@ --sample_name $@ --min_support 2
	flair fusion -g ${GENOME_FA} -f $(BASIC_ANNOTATION) -b ${FUSION3_TEST_BAM} --max_loci 3 -o $O/$@ --sample_name $@ --min_support 2
	diff $E/test-fusion-3gene.fusions.tsv $O/$@.fusions.tsv
	$(call cmp_expect_named,$@.fusions.aaseq.tsv,test-fusion-3gene.fusions.aaseq.tsv)
	$(call cmp_expect,$@.fusions.isoforms.bed)
	$(call cmp_expect_named,$@.fusions.isoforms.counts.tsv,test-fusion-3gene.fusions.isoforms.counts.tsv)
	$(call cmp_expect_named,$@.fusions.isoforms.fa,test-fusion-3gene.fusions.isoforms.fa)
	$(call cmp_expect,$@.fusions.isoforms.gtf)
	$(call cmp_expect_read_map_named,$@.fusions.isoforms.read.map.tsv,test-fusion-3gene.fusions.isoforms.read.map.tsv)
	$(call cmp_expect_named,$@.fusions.prelim.bed,test-fusion-3gene.fusions.prelim.bed)

# the same reads and options again: the chimeric BAMs are reused, and reported
test-fusion-reuse: mkdirs ${GENOME_FA} ${FUSION3_TEST_BAM}
	for f in genome_aligned transcriptome_aligned; do rm -f $O/$@.$$f.chim.bam $O/$@.$$f.chim.bam.bai $O/$@.$$f.chim.bam.provenance.json; done
	flair fusion -g ${GENOME_FA} -f $(BASIC_ANNOTATION) -b ${FUSION3_TEST_BAM} --max_loci 3 -o $O/$@ --sample_name $@ --min_support 2
	flair fusion -g ${GENOME_FA} -f $(BASIC_ANNOTATION) -b ${FUSION3_TEST_BAM} --max_loci 3 -o $O/$@ --sample_name $@ --min_support 2 2> $O/$@.log
	test $$(grep -c "reusing the genome alignment's chimeric reads" $O/$@.log) -eq 1
	test $$(grep -c 'reusing the chimeric alignments to the annotated transcripts' $O/$@.log) -eq 1
	diff $E/test-fusion-3gene.fusions.tsv $O/$@.fusions.tsv
	$(call cmp_expect_named,$@.fusions.aaseq.tsv,test-fusion-3gene.fusions.aaseq.tsv)
	$(call cmp_expect,$@.fusions.isoforms.bed)
	$(call cmp_expect_named,$@.fusions.isoforms.counts.tsv,test-fusion-3gene.fusions.isoforms.counts.tsv)
	$(call cmp_expect_named,$@.fusions.isoforms.fa,test-fusion-3gene.fusions.isoforms.fa)
	$(call cmp_expect,$@.fusions.isoforms.gtf)
	$(call cmp_expect_read_map_named,$@.fusions.isoforms.read.map.tsv,test-fusion-3gene.fusions.isoforms.read.map.tsv)
	$(call cmp_expect_named,$@.fusions.prelim.bed,test-fusion-3gene.fusions.prelim.bed)

test-fusion-3gene: mkdirs ${GENOME_FA} ${FUSION3_TEST_BAM}
	flair fusion -g ${GENOME_FA} -f $(BASIC_ANNOTATION) -b ${FUSION3_TEST_BAM} --max_loci 3 -o $O/$@ --sample_name $@ --min_support 2
	$(call cmp_expect,$@.fusions.isoforms.bed)
	$(call cmp_expect,$@.fusions.isoforms.fa)
	$(call cmp_expect,$@.fusions.tsv)
	$(call cmp_expect,$@.fusions.isoforms.counts.tsv)
	$(call cmp_expect,$@.fusions.aaseq.tsv)
	test -s $O/$@.fusions.isoforms.read.map.tsv -a -s $O/$@.fusions.isoforms.gtf -a -s $O/$@.fusions.isoforms.bb
	$(call cmp_expect,$@.fusions.isoforms.gtf)
	$(call cmp_expect_read_map,$@.fusions.isoforms.read.map.tsv)
	$(call cmp_expect,$@.fusions.prelim.bed)

${FUSION3_TEST_BAM}: $I/ERGIC3_DDX5_KRAS_sim_fusion.badread_30x.stranded.fastq ${GENOME_FA}
	$(call minimap2_genome_align,$I/ERGIC3_DDX5_KRAS_sim_fusion.badread_30x.stranded.fastq,${GENOME_FA},-s 40 -G 350k)

# the test-fusion annotation sorted by position, as minus strand exons are not in
# transcript order: the same fusion isoforms as test-fusion
test-fusion-sorted-gtf: mkdirs ${GENOME_FA} ${ERGIC3_KRAS_SIM_FA_GZ} ${FUSION_TEST_BAM}
	(grep '^#' $(BASIC_ANNOTATION); grep -v '^#' $(BASIC_ANNOTATION) | sort -k1,1 -k4,4n -k5,5nr) > $O/$@.annotation.gtf
	flair fusion -g ${GENOME_FA} -f $O/$@.annotation.gtf -b ${FUSION_TEST_BAM} --max_loci 4 -o $O/$@ --sample_name test-fusion --min_support 2
	diff <(sort $E/test-fusion.fusions.isoforms.bed) <(sort $O/$@.fusions.isoforms.bed)
	diff <(sort $E/test-fusion.fusions.isoforms.fa) <(sort $O/$@.fusions.isoforms.fa)
	diff $E/test-fusion.fusions.tsv $O/$@.fusions.tsv
	$(call cmp_expect_named,$@.fusions.aaseq.tsv,test-fusion.fusions.aaseq.tsv)
	$(call cmp_expect_named,$@.fusions.isoforms.counts.tsv,test-fusion.fusions.isoforms.counts.tsv)
	$(call cmp_expect_named,$@.fusions.isoforms.gtf,test-fusion.fusions.isoforms.gtf)
	$(call cmp_expect_read_map_named,$@.fusions.isoforms.read.map.tsv,test-fusion.fusions.isoforms.read.map.tsv)
	$(call cmp_expect_named,$@.fusions.prelim.bed,test-fusion.fusions.prelim.bed)

# a RefSeq annotation: RefSeq sequence names, gene symbols as gene ids and no
# gene_name attribute, BestRefSeq and Gnomon transcripts
test-fusion-refseq: mkdirs ${REFSEQ_GENOME_FA} ${REFSEQ_FUSION_TEST_BAM}
	flair fusion -g ${REFSEQ_GENOME_FA} -f $I/refseq.test.annotation.gtf -b ${REFSEQ_FUSION_TEST_BAM} --max_loci 4 -o $O/$@ --sample_name $@ --min_support 2
	$(call cmp_expect,$@.fusions.isoforms.bed)
	$(call cmp_expect,$@.fusions.isoforms.fa)
	$(call cmp_expect,$@.fusions.tsv)
	$(call cmp_expect_named,$@.fusions.aaseq.tsv,test-fusion.fusions.aaseq.tsv)
	$(call cmp_expect_named,$@.fusions.isoforms.counts.tsv,test-fusion.fusions.isoforms.counts.tsv)
	$(call cmp_expect,$@.fusions.isoforms.gtf)
	$(call cmp_expect_read_map,$@.fusions.isoforms.read.map.tsv)
	$(call cmp_expect,$@.fusions.prelim.bed)

${REFSEQ_FUSION_TEST_BAM}: ${ERGIC3_KRAS_SIM_FA_GZ} ${REFSEQ_GENOME_FA}
	$(call minimap2_genome_align,$I/ERGIC3_KRAS_sim_fusion.badread_20x.stranded.fastq,${REFSEQ_GENOME_FA},-s 40 -G 350k)

# genome chimeras only: no transcriptome alignment
test-fusion-no-tome: mkdirs ${GENOME_FA} ${ERGIC3_KRAS_SIM_FA_GZ} ${FUSION_TEST_BAM}
	flair fusion -g ${GENOME_FA} -f $(BASIC_ANNOTATION) -b ${FUSION_TEST_BAM} --max_loci 4 -o $O/$@ --sample_name $@ --min_support 2 --no_transcriptome_alignment
	$(call cmp_expect,$@.fusions.isoforms.bed)
	$(call cmp_expect,$@.fusions.isoforms.fa)
	$(call cmp_expect,$@.fusions.tsv)
	$(call cmp_expect_named,$@.fusions.aaseq.tsv,test-fusion.fusions.aaseq.tsv)
	$(call cmp_expect,$@.fusions.isoforms.counts.tsv)
	$(call cmp_expect,$@.fusions.isoforms.gtf)
	$(call cmp_expect_read_map,$@.fusions.isoforms.read.map.tsv)
	$(call cmp_expect,$@.fusions.prelim.bed)

${FUSION_TEST_BAM}: ${ERGIC3_KRAS_SIM_FA_GZ} ${GENOME_FA}
	$(call minimap2_genome_align,$I/ERGIC3_KRAS_sim_fusion.badread_20x.stranded.fastq,${GENOME_FA},-s 40 -G 350k)

test-fusion-help: mkdirs
	flair fusion  --help > $O/$@.out
	$(call cmp_expect_help,$@.out)

${ERGIC3_KRAS_SIM_FA_GZ}: $I/ERGIC3_KRAS_sim_fusion.badread_20x.stranded.fastq
	@mkdir -p $(dir $@)
	gzip -c $< >$@


###
# flair spliceevents
###
spliceevents-tests: test-spliceevents test-spliceevents-extraannot test-spliceevents-outliers test-spliceevents-region test-spliceevents-region-extend \
		test-spliceevents-region-partition test-spliceevents-smallparts test-spliceevents-nopartition test-spliceevents-restart test-spliceevents-compress \
		test-spliceevents-temp-dir test-spliceevents-keep-intermediate-with-temp-dir test-spliceevents-reuse-other-options \
		test-spliceevents-conditions test-spliceevents-help

# a rerun with --keep_intermediate and other correction options corrects the samples
# again rather than reading the first run's corrected reads, which are not what this
# run would correct, and gives what a run with nothing to reuse gives
test-spliceevents-reuse-other-options: mkdirs ${BASIC_ALIGN_BAM} ${SPLICEEVENTS_MANIFEST}
	flair spliceevents -m $(SPLICEEVENTS_MANIFEST) -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) --keep_intermediate -o $O/$@
	flair spliceevents -m $(SPLICEEVENTS_MANIFEST) -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) --keep_intermediate --ss_window 50 -o $O/$@ 2> $O/$@.log
	test $$(grep -c 'from other inputs or options' $O/$@.log) -eq 3
	flair spliceevents -m $(SPLICEEVENTS_MANIFEST) -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) --ss_window 50 -o $O/$@-clean
	diff $O/$@-clean.spliceevents.counts.tsv $O/$@.spliceevents.counts.tsv
	diff $O/$@-clean.spliceevents.psi_junc.tsv $O/$@.spliceevents.psi_junc.tsv
	$(call cmp_expect,$@.spliceevents.bed)
	$(call cmp_expect,$@.spliceevents.psi_tot.tsv)
	$(call cmp_expect_named,$@.spliceevents.regions.bed,test-spliceevents.spliceevents.regions.bed)
	$(call cmp_expect,$@.spliceevents.totals.tsv)

# a manifest of the four flair quantify fields: the conditions are written beside the
# counts as the sample info file flair diffexp reads, and the counts are those of the
# same samples named in the two field form (test-spliceevents-outliers)
test-spliceevents-conditions: mkdirs ${SPLICEEVENTS_MULTI_BAMS} $O/spliceevents.conditions.manifest.tsv
	flair spliceevents -m $O/spliceevents.conditions.manifest.tsv -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@
	$(call cmp_expect,$@.spliceevents.sample_info.tsv)
	diff $E/test-spliceevents-outliers.spliceevents.counts.tsv $O/$@.spliceevents.counts.tsv
	$(call cmp_expect_named,$@.spliceevents.bed,test-spliceevents-outliers.spliceevents.bed)
	$(call cmp_expect_named,$@.spliceevents.psi_junc.tsv,test-spliceevents-outliers.spliceevents.psi_junc.tsv)
	$(call cmp_expect_named,$@.spliceevents.psi_tot.tsv,test-spliceevents-outliers.spliceevents.psi_tot.tsv)
	$(call cmp_expect_named,$@.spliceevents.regions.bed,test-spliceevents-outliers.spliceevents.regions.bed)
	$(call cmp_expect_named,$@.spliceevents.totals.tsv,test-spliceevents-outliers.spliceevents.totals.tsv)

test-spliceevents-help: mkdirs
	flair spliceevents --help > $O/$@.out
	$(call cmp_expect_help,$@.out)

test-isoalleles-help: mkdirs
	flair isoalleles --help > $O/$@.out
	$(call cmp_expect_help,$@.out)

define spliceevents_diff
	$(call cmp_expect,$(1).spliceevents.bed)
	$(call cmp_expect,$(1).spliceevents.counts.tsv)
	$(call cmp_expect,$(1).spliceevents.psi_junc.tsv)
	$(call cmp_expect,$(1).spliceevents.psi_tot.tsv)
	$(call cmp_expect,$(1).spliceevents.totals.tsv)
	$(call cmp_expect,$(1).spliceevents.regions.bed)
endef

test-spliceevents: mkdirs ${BASIC_ALIGN_BAM} ${SPLICEEVENTS_MANIFEST}
	flair spliceevents -m $(SPLICEEVENTS_MANIFEST) -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@
	$(call spliceevents_diff,$@)

test-spliceevents-extraannot: mkdirs ${BASIC_ALIGN_BAM} ${SPLICEEVENTS_MANIFEST}
	flair spliceevents -m $(SPLICEEVENTS_MANIFEST) -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) --annot_basic $(BASIC_ANNOTATION) -o $O/$@
	$(call spliceevents_diff,$@)

test-spliceevents-outliers: mkdirs ${SPLICEEVENTS_MULTI_BAMS} ${SPLICEEVENTS_MULTI_MANIFEST}
	flair spliceevents -m $(SPLICEEVENTS_MULTI_MANIFEST) -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@
	$(call spliceevents_diff,$@)
	python3 ../src/flair/spliceevents_outliers.py -i $O/$@ -o $O/$@
	$(call cmp_expect,$@.spliceevents.outliers.tsv)
	$(call cmp_expect,$@.spliceevents.outliers.filtered.tsv)
	python3 ../src/flair/spliceevents_outliers.py -i $O/$@ -o $O/$@-psitot --psi tot
	$(call cmp_expect,$@-psitot.spliceevents.outliers.tsv)
	$(call cmp_expect,$@-psitot.spliceevents.outliers.filtered.tsv)

# a region cutting KRAS: of its events, only the exon skip lies inside, and its
# rows are those of the whole-genome run, though many of its reads run past the
# region's end; KRAS's end events, of all its reads, are left out.  A second region,
# without genes, adds nothing
test-spliceevents-region: mkdirs ${BASIC_ALIGN_BAM} ${SPLICEEVENTS_MANIFEST}
	flair spliceevents -m $(SPLICEEVENTS_MANIFEST) -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) --region_bed $I/spliceevents.kras_part.bed -o $O/$@
	diff <(grep -F 'es|chr12:25215437-25215560(-)' $E/test-spliceevents.spliceevents.counts.tsv) <(sed 1d $O/$@.spliceevents.counts.tsv)
	diff <(grep -F 'es|chr12:25215437-25215560(-)' $E/test-spliceevents.spliceevents.psi_junc.tsv) <(sed 1d $O/$@.spliceevents.psi_junc.tsv)
	$(call cmp_expect,$@.spliceevents.regions.bed)
	$(call cmp_expect,$@.spliceevents.bed)
	$(call cmp_expect,$@.spliceevents.psi_tot.tsv)
	$(call cmp_expect,$@.spliceevents.totals.tsv)

# a region at KRAS's 5' end: expanding it once over the reads that overlap it reaches
# 25251072, and over the reads that overlap those, 25251120, the whole gene
test-spliceevents-region-extend: mkdirs ${BASIC_ALIGN_BAM} ${SPLICEEVENTS_MANIFEST}
	flair spliceevents -m $(SPLICEEVENTS_MANIFEST) -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) --region=chr12:25205001-25205500 -o $O/$@
	$(call cmp_expect,$@.spliceevents.regions.bed)
	$(call cmp_expect,$@.spliceevents.bed)
	$(call cmp_expect,$@.spliceevents.counts.tsv)
	$(call cmp_expect,$@.spliceevents.psi_junc.tsv)
	$(call cmp_expect,$@.spliceevents.psi_tot.tsv)
	$(call cmp_expect,$@.spliceevents.totals.tsv)

# --region regions, the same as the --region_bed ones of test-spliceevents-region,
# give the region's events whatever --partition_items says, as a region is expanded
# over its reads and transcripts into one cluster
test-spliceevents-region-partition: mkdirs ${BASIC_ALIGN_BAM} ${SPLICEEVENTS_MANIFEST}
	flair spliceevents -m $(SPLICEEVENTS_MANIFEST) -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) --region=chr12:25200001-25235000 --region='chr12 60000000 60010000' --partition_items 1 -o $O/$@
	diff <(grep -F 'es|chr12:25215437-25215560(-)' $E/test-spliceevents.spliceevents.counts.tsv) <(sed 1d $O/$@.spliceevents.counts.tsv)
	diff <(grep -F 'es|chr12:25215437-25215560(-)' $E/test-spliceevents.spliceevents.psi_junc.tsv) <(sed 1d $O/$@.spliceevents.psi_junc.tsv)
	$(call cmp_expect,$@.spliceevents.bed)
	$(call cmp_expect,$@.spliceevents.psi_tot.tsv)
	$(call cmp_expect,$@.spliceevents.regions.bed)
	$(call cmp_expect,$@.spliceevents.totals.tsv)

# one partition per cluster of overlapping reads and transcripts, the most a
# chromosome can be cut into, gives the events of the default partitions: no read
# crosses a partition, so where the cuts fall changes nothing
test-spliceevents-smallparts: mkdirs ${BASIC_ALIGN_BAM} ${SPLICEEVENTS_MANIFEST}
	flair spliceevents -m $(SPLICEEVENTS_MANIFEST) -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) --partition_items 1 -o $O/$@
	for s in bed counts.tsv psi_junc.tsv psi_tot.tsv totals.tsv; do \
	    diff $E/test-spliceevents.spliceevents.$$s $O/$@.spliceevents.$$s || exit 1; done
	$(call cmp_expect,$@.spliceevents.regions.bed)
	$(call cmp_expect_named,$@.spliceevents.bed,test-spliceevents.spliceevents.bed)
	$(call cmp_expect_named,$@.spliceevents.counts.tsv,test-spliceevents.spliceevents.counts.tsv)
	$(call cmp_expect_named,$@.spliceevents.psi_junc.tsv,test-spliceevents.spliceevents.psi_junc.tsv)
	$(call cmp_expect_named,$@.spliceevents.psi_tot.tsv,test-spliceevents.spliceevents.psi_tot.tsv)
	$(call cmp_expect_named,$@.spliceevents.totals.tsv,test-spliceevents.spliceevents.totals.tsv)

# whole chromosomes, not cut, give the events of the cut ones
test-spliceevents-nopartition: mkdirs ${BASIC_ALIGN_BAM} ${SPLICEEVENTS_MANIFEST}
	flair spliceevents -m $(SPLICEEVENTS_MANIFEST) -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) --partition_items 0 -o $O/$@
	for s in bed counts.tsv psi_junc.tsv psi_tot.tsv totals.tsv; do \
	    diff $E/test-spliceevents.spliceevents.$$s $O/$@.spliceevents.$$s || exit 1; done
	$(call cmp_expect_named,$@.spliceevents.bed,test-spliceevents.spliceevents.bed)
	$(call cmp_expect_named,$@.spliceevents.counts.tsv,test-spliceevents.spliceevents.counts.tsv)
	$(call cmp_expect_named,$@.spliceevents.psi_junc.tsv,test-spliceevents.spliceevents.psi_junc.tsv)
	$(call cmp_expect_named,$@.spliceevents.psi_tot.tsv,test-spliceevents.spliceevents.psi_tot.tsv)
	$(call cmp_expect_named,$@.spliceevents.regions.bed,test-spliceevents.spliceevents.regions.bed)
	$(call cmp_expect_named,$@.spliceevents.totals.tsv,test-spliceevents.spliceevents.totals.tsv)

# a rerun with --keep_intermediate corrects only the samples without their corrected
# reads, here the one whose file is removed, and gives the same tables
test-spliceevents-restart: mkdirs ${BASIC_ALIGN_BAM} ${SPLICEEVENTS_MANIFEST}
	rm -rf $O/$@*
	flair spliceevents -m $(SPLICEEVENTS_MANIFEST) -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) --keep_intermediate -o $O/$@
	mkdir -p $O/$@.first && cp $O/$@.spliceevents.*.tsv $O/$@.first/
	rm $$(ls $O/$@.intermediate/*/*.read_chains.tsv | head -1)
	flair --log-level DEBUG spliceevents -m $(SPLICEEVENTS_MANIFEST) -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) --keep_intermediate -o $O/$@ 2> $O/$@.log
	test $$(grep -c 'reads to the annotated transcripts' $O/$@.log) -eq 1
	for f in $O/$@.first/*.tsv; do diff $$f $O/$@.$$(basename $$f | sed 's/^$@\.//') || exit 1; done
	$(call cmp_expect_named,$@.spliceevents.bed,test-spliceevents.spliceevents.bed)
	$(call cmp_expect_named,$@.spliceevents.counts.tsv,test-spliceevents.spliceevents.counts.tsv)
	$(call cmp_expect_named,$@.spliceevents.psi_junc.tsv,test-spliceevents.spliceevents.psi_junc.tsv)
	$(call cmp_expect_named,$@.spliceevents.psi_tot.tsv,test-spliceevents.spliceevents.psi_tot.tsv)
	$(call cmp_expect_named,$@.spliceevents.regions.bed,test-spliceevents.spliceevents.regions.bed)
	$(call cmp_expect_named,$@.spliceevents.totals.tsv,test-spliceevents.spliceevents.totals.tsv)

# --temp_dir contains the run's own directory, which is removed at the end, and nothing
# is written beside the output; the events are those of a run using $TMPDIR
test-spliceevents-temp-dir: mkdirs ${BASIC_ALIGN_BAM} ${SPLICEEVENTS_MANIFEST}
	mkdir -p $O/$@.tmp
	flair spliceevents -m $(SPLICEEVENTS_MANIFEST) -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) --temp_dir $O/$@.tmp -o $O/$@
	diff $E/test-spliceevents.spliceevents.counts.tsv $O/$@.spliceevents.counts.tsv
	test -z "$$(ls $O/$@.tmp)"
	test ! -e $O/$@.intermediate
	$(call cmp_expect_named,$@.spliceevents.bed,test-spliceevents.spliceevents.bed)
	$(call cmp_expect_named,$@.spliceevents.psi_junc.tsv,test-spliceevents.spliceevents.psi_junc.tsv)
	$(call cmp_expect_named,$@.spliceevents.psi_tot.tsv,test-spliceevents.spliceevents.psi_tot.tsv)
	$(call cmp_expect_named,$@.spliceevents.regions.bed,test-spliceevents.spliceevents.regions.bed)
	$(call cmp_expect_named,$@.spliceevents.totals.tsv,test-spliceevents.spliceevents.totals.tsv)

# the files kept are beside the output, so --temp_dir, which is accepted with it,
# names a directory they are not in and the run says so
test-spliceevents-keep-intermediate-with-temp-dir: mkdirs ${BASIC_ALIGN_BAM} ${SPLICEEVENTS_MANIFEST}
	mkdir -p $O/$@.tmp
	flair spliceevents -m $(SPLICEEVENTS_MANIFEST) -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) --keep_intermediate --temp_dir $O/$@.tmp -o $O/$@ > $O/$@.out 2>&1
	grep -q -- 'so --temp_dir is not used' $O/$@.out
	test -d $O/$@.intermediate
	test -z "$$(ls $O/$@.tmp)" 
	$(call cmp_expect_named,$@.spliceevents.bed,test-spliceevents.spliceevents.bed)
	$(call cmp_expect_named,$@.spliceevents.counts.tsv,test-spliceevents.spliceevents.counts.tsv)
	$(call cmp_expect_named,$@.spliceevents.psi_junc.tsv,test-spliceevents.spliceevents.psi_junc.tsv)
	$(call cmp_expect_named,$@.spliceevents.psi_tot.tsv,test-spliceevents.spliceevents.psi_tot.tsv)
	$(call cmp_expect_named,$@.spliceevents.regions.bed,test-spliceevents.spliceevents.regions.bed)
	$(call cmp_expect_named,$@.spliceevents.totals.tsv,test-spliceevents.spliceevents.totals.tsv)

$O/spliceevents.s%.bam: ${BASIC_ALIGN_BAM}
	samtools view -b -s $*.6 -o $@ $<
	samtools index $@

$O/spliceevents.%.bam: $I/spliceevents.%.reads.txt ${BASIC_ALIGN_BAM}
	samtools view -b -N $< -o $@ ${BASIC_ALIGN_BAM}
	samtools index $@


##
# flair_partition command
##
partition-tests: test-partition-basic test-partition-adjacent test-partition-min-reads \
	test-partition-sam test-partition-region test-partition-region-opt

test-partition-basic: mkdirs
	flair_partition --threads=2 --log-level=info --bed=input/partition.test.bed $O/$@.ranges.bed >& $O/$@.out
	$(call cmp_expect,$@.ranges.bed)
	$(call cmp_expect_log,$@.out)

# merge adjacent partitions
test-partition-adjacent: mkdirs
	flair_partition --threads=2 --part_merge_dist=10000 --log-level=info --bed=input/partition.test.bed $O/$@.ranges.bed >& $O/$@.out
	$(call cmp_expect,$@.ranges.bed)
	$(call cmp_expect_log,$@.out)

# minimum reads per partition
test-partition-min-reads: mkdirs
	flair_partition --threads=2 --min_partition_items=3  --log-level=info --bed=input/partition.test.bed $O/$@.ranges.bed >& $O/$@.out
	$(call cmp_expect,$@.ranges.bed)
	$(call cmp_expect_log,$@.out)

# --region_bed grows each region over the items reaching it, so the partitions are
# those of the whole input that the regions fall in
test-partition-region: mkdirs
	flair_partition --threads=2 --log-level=info --bed=input/partition.test.bed --region_bed=input/partition.region.bed \
	    $O/$@.ranges.bed >& $O/$@.out
	$(call cmp_expect,$@.ranges.bed)
	$(call cmp_expect_log,$@.out)

# --region, in either coordinate form, gives the partitions of the same --region_bed
# regions, so the comparison is with that test's expected files
test-partition-region-opt: mkdirs
	flair_partition --threads=2 --log-level=info --bed=input/partition.test.bed \
	    --region=chr19:368001-368100 --region='chr6 74285000 74285100' \
	    $O/$@.ranges.bed >& $O/$@.out
	diff $E/test-partition-region.ranges.bed $O/$@.ranges.bed
	sed -i -r -e 's/^[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2} //' $O/$@.out
	diff $E/test-partition-region.out $O/$@.out

# with bed, sam, and gtf
test-partition-sam: mkdirs
	flair_partition --threads=2 --log-level=info --bed=input/partition.test.bed --bam=input/tiny.sam --gtf=input/seg1.gencode_v47.gtf \
	    $O/$@.ranges.bed >& $O/$@.out
	$(call cmp_expect,$@.ranges.bed)
	$(call cmp_expect_log,$@.out)

##
# library tests with pytest.  test_flair_partition.py reads the basic alignment
# BAM, so a parallel run after `make clean' must not start pytest before
# test-alignments has built it.
##
lib-tests: mkdirs ${BASIC_ALIGN_BAM}
	${PYTEST_CMD} .

##
# build test input files
##

##
# The generated inputs under ${TI}.  They do not vary by python version, so runs
# in separate output directories share them and must build them first with `make
# test-inputs': two runs creating the same genome at once would corrupt it.  A
# new ${TI} file belongs in this list.
##
TEST_INPUT_FILES = ${GENOME_FA} ${REFSEQ_GENOME_FA} ${ONECHROM_GENOME_FA} ${BGZIP_GENOME_FA} \
	${COMBINE_BEDS_GZ} \
	${TI}/basic.annotation.gtf.gz ${TI}/basic.counts_matrix.diffexp.tsv.gz \
	${TI}/basic.counts_matrix.diffsplice.tsv.gz ${TI}/basic.isoforms.diffsplice.bed.gz \
	${TI}/basic.shortread_junctions.tab.gz ${TI}/test-collapse-annot.isoforms.bed.gz \
	${TI}/test-quantify.isoforms.counts.tsv.gz

test-inputs: ${TEST_INPUT_FILES}

CHROM_URL = https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes

# The chromosomes are downloaded and joined in a directory of their own, and the
# genome and its index are moved into place complete, so that two runs that both
# find the genome missing cannot write over each other's partial files.
${GENOME_FA}:
	@mkdir -p $(dir $@)
	tmpdir=$$(mktemp -d $(dir $@)/genome.XXXXXX) \
	    && (cd $$tmpdir && wget -nv ${CHROM_URL}/chr12.fa.gz ${CHROM_URL}/chr17.fa.gz ${CHROM_URL}/chr20.fa.gz) \
	    && zcat $$tmpdir/chr12.fa.gz $$tmpdir/chr17.fa.gz $$tmpdir/chr20.fa.gz > $$tmpdir/genome.fa \
	    && samtools faidx $$tmpdir/genome.fa \
	    && mv -f $$tmpdir/genome.fa.fai $@.fai \
	    && mv -f $$tmpdir/genome.fa $@ \
	    && rm -rf $$tmpdir

# one sequence, re-wrapped at 50 bases so it can be indexed
${REFSEQ_GENOME_FA}: ${GENOME_FA}
	sed -e 's/^>chr12$$/>NC_000012.12/' -e 's/^>chr17$$/>NC_000017.11/' -e 's/^>chr20$$/>NC_000020.11/' ${GENOME_FA} > $@.$$$$.tmp && mv -f $@.$$$$.tmp $@
	samtools faidx $@

${ONECHROM_GENOME_FA}: ${GENOME_FA}
	awk 'BEGIN {print ">chr12_17_20"} /^>/ {next} {buf = buf $$0; while (length(buf) >= 50) {print substr(buf, 1, 50); buf = substr(buf, 51)}} END {if (buf != "") print buf}' ${GENOME_FA} > $@.$$$$.tmp && mv -f $@.$$$$.tmp $@
	samtools faidx $@

${TOME_OUT_FA}: | mkdirs
	cp $E/test-transcriptome.isoforms.fa $@

${TOME_OUT_FAI}: ${TOME_OUT_FA}
	samtools faidx --output $@ $<

$O/%.genome_aligned.bam: $I/%.fa ${GENOME_FA} | mkdirs
	minimap2 -ax splice --secondary=no ${GENOME_FA} $< | samtools view -hb - | samtools sort - > $@.tmp.bam
	mv -f $@.tmp.bam $@
	samtools index $@

##
# test hub (work in progress)
# 
#  make hub-check huburl=https://hgwdev.gi.ucsc.edu/~markd/flair/flair-mddev/test/hub/hub.txt
##
HG38_CHROMS = https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromAlias.bb

TEST_TRANSCRIPTOME_MUTANT_SPLICE_BB = $O/test-transcriptome-mutant-splice.isoforms.bb
TEST_TRANSCRIPTOME_SEG1_BB =  $O/test-transcriptome-seg1.isoforms.bb

test-hub:  ${KRAS_MUT_SPLICE_BB} ${KRAS_MUT_SPLICE_SJ_BB} ${TEST_TRANSCRIPTOME_MUTANT_SPLICE_BB} \
	${TEST_TRANSCRIPTOME_SEG1_BB}

$O/%.bb: $O/%.bed
	bedToBigBed -sort -sizesIsChromAliasBb $< ${HG38_CHROMS} $@

$O/%.sj.bb: $I/%.sj.bed
	@mkdir -p $O
	bedToBigBed -sort -type=bed6 -sizesIsChromAliasBb $< ${HG38_CHROMS} $@

$O/%.bb: $I/%.bed
	@mkdir -p $O
	bedToBigBed -sort -sizesIsChromAliasBb $< ${HG38_CHROMS} $@

hub-check: test-hub
	@[ -n "${huburl}" ] || { echo "Error: huburl= must be specified on the commad line"; exit 1; }
	hubCheck ${huburl}

mkdirs:
	@mkdir -p $O

clean:
	rm -rf $O __pycache__
	rm -rf [0-9a-f]*-[0-9a-f]*-[0-9a-f]*-[0-9a-f]*-[0-9a-f]*

real-clean : clean
	rm -rf ${TI}

##
# Compression: every subcommand run with -z, its outputs compared against the
# expected files of its plain test
##

test-transcriptome-compress: mkdirs ${BASIC_ALIGN_BAM} ${TI}/basic.annotation.gtf.gz ${TI}/basic.shortread_junctions.tab.gz
	flair -z transcriptome -b ${BASIC_ALIGN_BAM} -g ${GENOME_FA} --junction_star ${TI}/basic.shortread_junctions.tab.gz \
	   -t 4 -f ${TI}/basic.annotation.gtf.gz -o $O/$@ --sample_name test-transcriptome --generate_map
	$(call cmp_expect_gz,$@.isoforms.bed,test-transcriptome.isoforms.bed)
	$(call cmp_expect_gz,$@.isoforms.gtf,test-transcriptome.isoforms.gtf)
	$(call cmp_expect_gz,$@.isoforms.fa,test-transcriptome.isoforms.fa)
	$(call cmp_expect_gz,$@.isoforms.read.map.tsv,test-transcriptome.isoforms.read.map.tsv)
	$(call cmp_expect_gz,$@.regions.bed,test-transcriptome.regions.bed)
	$(call check_fasta_indexable,$@.isoforms.fa)

test-spliceevents-compress: mkdirs ${BASIC_ALIGN_BAM} test-spliceevents \
		$O/spliceevents.manifest.tsv.gz ${TI}/basic.annotation.gtf.gz
	flair spliceevents -m $O/spliceevents.manifest.tsv.gz -g ${GENOME_FA} -t 4 \
	   -f ${TI}/basic.annotation.gtf.gz -o $O/$@
	$(call cmp_expect_sort_as,$@.spliceevents.counts.tsv,test-spliceevents.spliceevents.counts.tsv)
	$(call cmp_expect_sort_as,$@.spliceevents.totals.tsv,test-spliceevents.spliceevents.totals.tsv)
	$(call cmp_expect_sort_as,$@.spliceevents.bed,test-spliceevents.spliceevents.bed)

