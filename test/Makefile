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
	@echo "   " make align-tests >&2
	@echo to installed FLAIR for individual tests use the use_installed_flair option, such as: >&2
	@echo "   " make align-tests use_installed_flair=yes >&2
	@echo "Run only library tests" >&2
	@echo "   " make lib-tests >&2
	@echo "Remove output files:" >&2
	@echo "   " make clean >&2
	@echo Remove output and downloaded genomes files: >&2
	@echo "   " make real-clean >&2
	@echo run flake8: >&2
	@echo "   " make flake8 >&2
	@echo "build test case hub" >&2
	@echo "   " make test-hub >&2
	@echo "   " make make hub-check huburl=https://public.gi.ucsc.edu/~you/flair/test/hub/hub.txt >&2

test:
	${MAKE} test-base test-with-R

test-base:
	${MAKE} align-tests quantify-tests bed-to-gtf-tests \
		combine-tests diff-iso-usage-tests \
		plot-usage-tests transcriptome-tests fusion-tests variants-tests \
		spliceevents-tests partition-tests lib-tests test-doc-cli

test-no-align:
	${MAKE} quantify-tests bed-to-gtf-tests \
		combine-tests diff-iso-usage-tests \
		plot-usage-tests transcriptome-tests fusion-tests variants-tests \
		spliceevents-tests partition-tests lib-tests

test-with-R:
	${MAKE} diffexp-tests diffsplice-tests

test-installed:
	${MAKE} test use_installed_flair=yes

test-base-installed:
	${MAKE} test-base use_installed_flair=yes

test-help:
	${MAKE} test-main-help test-align-help \
		test-quantify-help test-diffexp-help test-diffsplice-help \
		test-diff-iso-usage-help \
		test-plot-usage-help test-transcriptome-help test-fusion-help \
		test-doc-cli

# the committed option pages must match what the parsers currently define
test-doc-cli: mkdirs
	rm -rf $O/cli
	${root}/dev/bin/flair-gen-cli-doc $O/cli >/dev/null
	diff -r ${root}/docs/source/cli $O/cli

test-fast:
	${MAKE} test-combine-bed-only \
		test-diff-iso-usage test-plot-usage test-transcriptome \
		test-fusion

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

I = input
E = expected
O = output
TI = tmp-input

# reference data
GENOME_FA = ${TI}/genome.fa
# the test genome with all chromosomes concatenated into one, chr12_17_20
ONECHROM_GENOME_FA = ${TI}/genome.onechrom.fa

##
# basic: general tests that programs seem to work
##
BASIC_READS_FA = $I/basic.reads.fa
BASIC_ANNOTATION = $I/basic.annotation.gtf
BASIC_ANNOTATION_INCOMPLETE = $I/basic.annotation.incomplete.gtf
BASIC_SHORTREAD_JUNCTIONS = $I/basic.shortread_junctions.tab
BASIC_PROMOTER_REGIONS = $I/basic.promoter_regions.bed
BASIC_READS_MANIFEST = $I/basic.reads_manifest.txt

# FIXME: test are based on old test data, however flair quantify
# does not produce sufficient pvalues to work diffexpr
BASIC_COUNTS_MATRIX = $I/basic.counts_matrix_diffexp.tsv
BASIC_DIFFSPLICE_MATRIX = $I/basic.counts_matrix_diffsplice.tsv
# isoforms matching BASIC_DIFFSPLICE_MATRIX; test-collapse-annot.isoforms.bed
# is a different generation of the data and does not match the counts
BASIC_DIFFSPLICE_BED = $I/basic.isoforms_diffsplice.bed
COLLAPSE_ANNOT_ISOFORMS_BED = $I/test-collapse-annot.isoforms.bed

TOME_OUT_FA = $O/test-transcriptome.isoforms.fa
TOME_OUT_FAI = $O/test-transcriptome.isoforms.fa.fai

QUANTIFY_INPUT_READS = basic.reads.1 basic.reads.2 basic.reads.3 \
                       basic.reads.4 basic.reads.5 basic.reads.6
QUANTIFY_INPUT_BAMS = ${QUANTIFY_INPUT_READS:%=$O/%.genomealigned.bam}

##
# Files generate by tests that are used by other tests.
# These enfoce files being create in the right order.
## 
TEST_ALIGN_BAM = $O/test-align.bam
# basic reads with poly(A) tails removed, so single-exon reads have no poly(A)
# evidence for strand correction
BASIC_READS_NOPOLYA_FA = $O/basic.reads.nopolyA.fa
NOPOLYA_ALIGN_BAM = $O/test-align-nopolya.bam
# poly(A)-trimmed basic reads with a poly(T) head added to every read, so single-exon
# reads from minus strand genes look like plus strand reads with a poly(A) tail
BASIC_READS_POLYT_FA = $O/basic.reads.polyT.fa
POLYT_ALIGN_BAM = $O/test-align-polyt.bam
# basic annotation plus a synthetic gene that shares one DDX5 intron
SHARED_JUNC_GENE_GTF = $I/shared_junction_gene.gtf
SHARED_JUNC_ANNOTATION = $O/basic.annotation.sharedjunc.gtf
# test-align.bam with only one chr20 read, so the chr20 partition has reads but
# too few to support any isoform
ONECHR20_ALIGN_BAM = $O/test-align-onechr20.bam
# basic reads aligned to the one-chromosome genome
ONECHROM_ALIGN_BAM = $O/test-align-onechrom.bam
TEST_QUANTIFY_B1_BAM = $O/test-quantify.B1.B.flair.aligned.bam

###
# seg1: single exon tests
###
SEG1_READS_FQ = $I/seg1.cdna-ont.fastq
SEG1_ANNOTATION = $I/seg1.gencodeV47.gtf
SEG1_PROMOTER_REGIONS = $I/seg1.promoter-regions.bed

SEG1_ALIGN_BAM = $O/test-align-seg1.bam

##
# mutant splice site, BED is truth
##
DDX5_READTHROUGH_FQ = $I/DDX5_readthrough.badread40x.fastq
DDX5_READTHROUGH_BAM = $O/DDX5_readthrough_sim_badread40x.bam
DDX5_READTHROUGH_ANNOT = $I/DDX5_POLG2_gencode38.gtf
DDX5_READTHROUGH_SJ = $I/DDX5_readthrough_IPSJ.bed
# DDX5 readthrough reads with all but 3 reads of the two supersets of the
# chr17:64500697-64506289 subset isoform dropped, so that subset outnumbers them
DDX5_BYSUPPORT_DROP_READS = $I/DDX5_readthrough_bysupport_drop_reads.txt
DDX5_BYSUPPORT_BAM = $O/DDX5_readthrough_bysupport.bam
# DDX5/POLG2 annotation plus the readthrough transcript as its own gene
DDX5_READTHROUGH_GENE_GTF = $I/DDX5_POLG2_readthrough_gene.gtf
DDX5_READTHROUGH_GENE_ANNOT = $O/DDX5_POLG2_gencode38.readthrough_gene.gtf

KRAS_MUT_SPLICE_FQ = $I/KRAS_sim_genomic_mut_driven_splicing.badreadsimPB.fastq
KRAS_MUT_SPLICE_BAM = $O/KRAS_sim_genomic_mut_driven_splicing.badreadsimPB.bam
KRAS_MUT_SPLICE_SJ_BED = $I/KRAS_sim_genomic_mut_driven_splicing.badreadsimPB.sj.bed
KRAS_MUT_SPLICE_SJ_BB = $O/KRAS_sim_genomic_mut_driven_splicing.badreadsimPB.sj.bb
KRAS_MUT_SPLICE_BED = $I/KRAS_sim_genomic_mut_driven_splicing.bed
KRAS_MUT_SPLICE_BB = $O/KRAS_sim_genomic_mut_driven_splicing.bb

COMBINE_MANIFEST = $I/combine_manifest_bed-only.txt
QUANTIFY_MANIFEST = $I/basic.reads_manifest_bam.txt
SPLICEEVENTS_MANIFEST = $I/spliceevents.manifest.tsv

TOME_ISOFORM_BED = $E/test-transcriptome.isoforms.bed
TOME_READ_MAP = $E/test-transcriptome.isoform.read.map.txt

# The help test results vary on screen width.  This forces them to be consistent.
# It might seem weird to see if help works, but some it did break in the past.
export COLUMNS = 80

##
# compare output, expected with result
##
define cmp_expect
	diff $E/$(1) $O/$(1)
endef

define cmp_expect_bed_simple
	diff <(cut -f 1,2,3,6,10,11,12 $E/$(1)) <(cut -f 1,2,3,6,10,11,12 $O/$(1))
endef

##
# Expected and output are sorted in cases where the order is not stable.  In
# this case, we first test for the file as sorted input subsitution failure
# isn't propigate by shell.
##
define cmp_expect_sort
	@test -e $O/${1} || (echo Error: output file missing: $O/${1}; exit 1)
	diff <(sort $E/${1}) <(sort $O/${1})
endef

###
# MAIN COMMAND
###

test-main-help: mkdirs
	flair --help | fgrep -v Version> $O/$@.out
	$(call cmp_expect,$@.out)

###
# FLAIR ALIGN
###

align-tests: test-align test-align-seg1 test-align-nopolya test-align-help

# just checks bed4 to make sure we are getting alignments
test-align : mkdirs ${GENOME_FA}
	flair align -r $(BASIC_READS_FA) --genome ${GENOME_FA} -t 4 -o $O/$@
	bedtools bamtobed -i $O/$@.bam > $O/$@.bed
	$(call cmp_expect,$@.bed)

${TEST_ALIGN_BAM}: test-align

${BASIC_READS_NOPOLYA_FA}: ${BASIC_READS_FA} bin/trim_polya_reads
	@mkdir -p $O
	bin/trim_polya_reads ${BASIC_READS_FA} $@.tmp
	mv -f $@.tmp $@

${SEG1_ALIGN_BAM}: test-align-seg1
${NOPOLYA_ALIGN_BAM}: test-align-nopolya

${BASIC_READS_POLYT_FA}: ${BASIC_READS_NOPOLYA_FA}
	awk '/^>/ {print; next} {print "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" $$0}' ${BASIC_READS_NOPOLYA_FA} > $@.tmp
	mv -f $@.tmp $@

${POLYT_ALIGN_BAM}: ${BASIC_READS_POLYT_FA} ${GENOME_FA}
	flair align -r ${BASIC_READS_POLYT_FA} --genome ${GENOME_FA} -t 4 -o $(basename $@)

${ONECHR20_ALIGN_BAM}: ${TEST_ALIGN_BAM}
	(samtools view -h ${TEST_ALIGN_BAM} chr12 chr17 && samtools view -F 0x900 ${TEST_ALIGN_BAM} chr20 | head -n 1) | samtools view -b -o $@.tmp.bam -
	mv -f $@.tmp.bam $@
	samtools index $@

${SHARED_JUNC_ANNOTATION}: ${BASIC_ANNOTATION} ${SHARED_JUNC_GENE_GTF}
	@mkdir -p $O
	cat ${BASIC_ANNOTATION} ${SHARED_JUNC_GENE_GTF} > $@.tmp
	mv -f $@.tmp $@

test-align-seg1 : mkdirs ${GENOME_FA}
	flair align -r ${SEG1_READS_FQ} --genome ${GENOME_FA} -t 4 -o $O/$@
	bedtools bamtobed -i $O/$@.bam > $O/$@.bed
	$(call cmp_expect,$@.bed)

# basic reads with poly(A) tails trimmed
test-align-nopolya : mkdirs ${GENOME_FA} ${BASIC_READS_NOPOLYA_FA}
	flair align -r ${BASIC_READS_NOPOLYA_FA} --genome ${GENOME_FA} -t 4 -o $O/$@
	bedtools bamtobed -i $O/$@.bam > $O/$@.bed
	$(call cmp_expect,$@.bed)

test-align-help: mkdirs
	flair align --help > $O/$@.out
	$(call cmp_expect,$@.out)

# macros to do comparions with expected results
# FIXME: make sort version go away by making results stable
define diff_bed_gtf_sort
	$(call cmp_expect_sort,$@.isoforms.bed)
	$(call cmp_expect_sort,$@.isoforms.gtf)
endef
define diff_bed_gtf
	$(call cmp_expect,$@.isoforms.bed)
	$(call cmp_expect,$@.isoforms.gtf)
endef

define collapse_promoter_diff
	$(call cmp_expect_sort,$@.promoter_supported.bed)
endef

# resulting read map differs between run (see #443)
# just check if it exists for now
# should do: $(call cmp_expect_sort,$@.isoform.read.map.txt)
define diff_bed_gtf_sort_read_map
	@test -e $O/${1} || (echo Error: output file missing: $O/${1}; exit 1)
endef

###
# FLAIR COMBINE
###

define combine_diff
	$(call cmp_expect,$@.combined.isoforms.bed)
	$(call cmp_expect,$@.combined.isoform.counts.tsv)
endef


combine-tests: test-combine test-combine-files test-combine-nose test-combine-10pct test-combine-multends test-combine-multends-10pct

update-combine-tests:
	for tname in "test-combine" "test-combine-files" "test-combine-nose" "test-combine-10pct" "test-combine-multends" "test-combine-multends-10pct"; \
	do \
	cp $O/$$tname.combined.isoforms.bed $E/$$tname.combined.isoforms.bed; \
	cp $O/$$tname.combined.isoform.counts.tsv $E/$$tname.combined.isoform.counts.tsv; \
	done

test-combine: mkdirs
	flair combine --manifest $(COMBINE_MANIFEST) --genome ${GENOME_FA} -o $O/$@
# 	$(call cmp_expect_bed_simple,$@.combined.isoforms.bed)
	$(combine_diff)

test-combine-files: mkdirs
	flair combine --genome ${GENOME_FA} -o $O/$@ --isoform_beds expected/test-transcriptome.isoforms.bed,expected/test-transcriptome-mult-ends.isoforms.bed,expected/test-transcriptome-mutant-splice.isoforms.bed,expected/test-transcriptome-readthrough.isoforms.bed,expected/test-fusion.fusions.isoforms.bed,expected/test-transcriptome-seg1-badannot.isoforms.bed,expected/test-transcriptome-seg1-badannot-multends.isoforms.bed,input/test-transcriptome-seg1-badannot-multends-chromshifted.isoforms.bed
# 	$(call cmp_expect_bed_simple,$@.combined.isoforms.bed)
	$(combine_diff)

test-combine-nose: mkdirs
	flair combine --manifest $(COMBINE_MANIFEST) --genome ${GENOME_FA} -o $O/$@ --remove_single_exon
# 	$(call cmp_expect_bed_simple,$@.combined.isoforms.bed)
	$(combine_diff)

test-combine-10pct: mkdirs
	flair combine --manifest $(COMBINE_MANIFEST) --genome ${GENOME_FA} -o $O/$@ --min_percent_usage 10
# 	$(call cmp_expect_bed_simple,$@.combined.isoforms.bed)
	$(combine_diff)

test-combine-multends: mkdirs
	flair combine --manifest $(COMBINE_MANIFEST) --genome ${GENOME_FA} -o $O/$@ --max_ends 10
# 	$(call cmp_expect_bed_simple,$@.combined.isoforms.bed)
	$(combine_diff)

test-combine-multends-10pct: mkdirs
	flair combine --manifest $(COMBINE_MANIFEST) --genome ${GENOME_FA} -o $O/$@ --max_ends 10 --min_percent_usage 10
# 	$(call cmp_expect_bed_simple,$@.combined.isoforms.bed)
	$(combine_diff)

###
# FLAIR QUANTIFY
###
quantify-tests: test-quantify test-quantify-normends test-quantify-help

test-quantify: mkdirs ${TOME_OUT_FAI} ${QUANTIFY_INPUT_BAMS} ${GENOME_FA}
	rm -rf $O/$@.intermediate
	flair quantify --manifest $(QUANTIFY_MANIFEST) --genome ${GENOME_FA} --isoform_bed $(TOME_ISOFORM_BED) --generate_map --tpm --sample_id_only -o $O/$@ --with_gene
	$(call cmp_expect,$@.counts.tsv)
	$(call cmp_expect,$@.tpm.tsv)

test-quantify-normends: mkdirs ${TOME_OUT_FAI} ${QUANTIFY_INPUT_BAMS} ${GENOME_FA}
	rm -rf $O/$@.intermediate
	flair quantify --manifest $(QUANTIFY_MANIFEST) --genome ${GENOME_FA} --isoform_bed $(TOME_ISOFORM_BED) --generate_map --tpm --sample_id_only -o $O/$@ --with_gene --normalize_ends
	$(call cmp_expect,$@.counts.tsv)
	$(call cmp_expect,$@.tpm.tsv)

test-quantify-help: mkdirs
	flair quantify --help > $O/$@.out
	$(call cmp_expect,$@.out)

###
# FLAIR DIFFEXP
###
diffexp-tests: test-diffexp test-diffexp-help

test-diffexp : mkdirs
	flair diffexp --counts_matrix $(BASIC_COUNTS_MATRIX) -o $O/$@ --min_expression 1 --overwrite_output
	$(call cmp_expect_sort,$@/genes_deseq2_MCF7_v_A549.tsv)
	$(call cmp_expect_sort,$@/isoforms_deseq2_MCF7_v_A549.tsv)
	$(call cmp_expect_sort,$@/isoforms_drimseq_MCF7_v_A549.tsv)

test-diffexp-help: mkdirs
	flair diffexp --help > $O/$@.out
	$(call cmp_expect,$@.out)

###
# FLAIR DIFFSPLICE
###
diffsplice-tests: test-diffsplice test-diffsplice-help

# various output files not long produce due to changes in test cases
test-diffsplice: mkdirs
	flair diffsplice --isoform_bed $(BASIC_DIFFSPLICE_BED) --counts_matrix $(BASIC_DIFFSPLICE_MATRIX) --test -o $O/$@ --overwrite_output
	$(call cmp_expect_sort,$@/diffsplice.alt3.events.quant.tsv)
	$(call cmp_expect_sort,$@/diffsplice.alt5.events.quant.tsv)
	$(call cmp_expect,$@/diffsplice.es.events.quant.tsv)
	$(call cmp_expect_sort,$@/diffsplice.ir.events.quant.tsv)
	$(call cmp_expect,$@/drimseq_alt3_A_v_B.tsv)
	$(call cmp_expect,$@/drimseq_es_A_v_B.tsv)

test-diffsplice-cond: mkdirs
	flair diffsplice --isoform_bed $(BASIC_DIFFSPLICE_BED) --counts_matrix $(BASIC_DIFFSPLICE_MATRIX) --test -o $O/$@ --overwrite_output --condition_a A --condition_b B
	$(call cmp_expect_sort,test-diffsplice/diffsplice.alt3.events.quant.tsv)
	$(call cmp_expect_sort,test-diffsplice/diffsplice.alt5.events.quant.tsv)
	$(call cmp_expect,test-diffsplice/diffsplice.es.events.quant.tsv)
	$(call cmp_expect_sort,test-diffsplice/diffsplice.ir.events.quant.tsv)
	$(call cmp_expect,test-diffsplice/drimseq_alt3_A_v_B.tsv)
	$(call cmp_expect,test-diffsplice/drimseq_es_A_v_B.tsv)

test-diffsplice-help: mkdirs
	flair diffsplice --help > $O/$@.out
	$(call cmp_expect,$@.out)


bed-to-gtf-tests: test-bed-to-gtf test-bed-to-gtf-extracols

test-bed-to-gtf : mkdirs
	bed_to_gtf $(COLLAPSE_ANNOT_ISOFORMS_BED) > $O/$@.gtf
	$(call cmp_expect,$@.gtf)

test-bed-to-gtf-extracols : 
	bed_to_gtf --is_flair_bed $(TOME_ISOFORM_BED) > $O/$@.gtf
	$(call cmp_expect,$@.gtf)


##
# DIFF ISO USAGE
##
diff-iso-usage-tests: test-diff-iso-usage test-diff-iso-usage-help

test-diff-iso-usage: mkdirs
	diff_iso_usage $(BASIC_COUNTS_MATRIX) MCF7cdnar1r3_MCF7_batch1 MCF7cdnastrandr3r2_MCF7_batch1 $O/$@.tsv
	$(call cmp_expect,$@.tsv)

test-diff-iso-usage-help: mkdirs
	diff_iso_usage --help > $O/$@.out
	$(call cmp_expect,$@.out)

# ##
# # plot isoform usage
# ##
plot-usage-tests: test-plot-usage test-plot-usage-help

# only works with quantify output run with with-gene
test-plot-usage: mkdirs
	plot_isoform_usage $(TOME_ISOFORM_BED) $E/test-quantify.counts.tsv FLG00000002 --min_reads 1 -o $O/$@
	file $O/test-plot-usage_FLG00000002_usage.pdf | grep 'PDF document' >/dev/null
	file $O/test-plot-usage_FLG00000002_isoforms.pdf | grep 'PDF document' >/dev/null

test-plot-usage-help: mkdirs
	plot_isoform_usage --help > $O/$@.out
	$(call cmp_expect,$@.out)

###
# flair transcriptome
###

transcriptome-tests: test-transcriptome test-transcriptome-nochecksplice test-transcriptome-nostringent test-transcriptome-nostringentsplice test-transcriptome-directrna test-transcriptome-partition test-transcriptome-lowannot test-transcriptome-noannot test-transcriptome-normends test-transcriptome-comprehensive test-transcriptome-comprehensive-basic \
                    test-transcriptome-seg1 test-transcriptome-mutant-splice test-transcriptome-readthrough test-transcriptome-readthrough-gene test-transcriptome-bysupport test-transcriptome-nopolya test-transcriptome-nopolya-trust-strand test-transcriptome-polyt test-transcriptome-shared-junc test-transcriptome-no-firstpass test-transcriptome-onechrom test-transcriptome-mult-ends \
					test-transcriptome-trust-ends test-transcriptome-mult-ends-trust-ends test-transcriptome-seg1-badannot test-transcriptome-seg1-badannot-multends test-transcriptome-seg1-badannot-comprehensive test-transcriptome-help

update-transcriptome-tests:
	for tname in "test-transcriptome" "test-transcriptome-nochecksplice" "test-transcriptome-nostringent" "test-transcriptome-nostringentsplice" "test-transcriptome-directrna" "test-transcriptome-partition" "test-transcriptome-lowannot" "test-transcriptome-noannot" "test-transcriptome-normends" "test-transcriptome-comprehensive" "test-transcriptome-comprehensive-basic" \
	"test-transcriptome-seg1" "test-transcriptome-mutant-splice" "test-transcriptome-readthrough" "test-transcriptome-readthrough-gene" "test-transcriptome-bysupport" "test-transcriptome-nopolya" "test-transcriptome-nopolya-trust-strand" "test-transcriptome-polyt" "test-transcriptome-shared-junc" "test-transcriptome-no-firstpass" "test-transcriptome-onechrom" "test-transcriptome-mult-ends" \
	"test-transcriptome-trust-ends" "test-transcriptome-mult-ends-trust-ends" "test-transcriptome-seg1-badannot" "test-transcriptome-seg1-badannot-multends" "test-transcriptome-seg1-badannot-comprehensive"; \
	do \
	cp $O/$$tname.isoforms.bed $E/$$tname.isoforms.bed; \
	cp $O/$$tname.isoforms.gtf $E/$$tname.isoforms.gtf; \
	done


test-transcriptome: mkdirs ${TEST_ALIGN_BAM}
	flair transcriptome -b ${TEST_ALIGN_BAM} -g ${GENOME_FA} --junction_tab ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --keep_intermediate --generate_map
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)
	$(call cmp_expect,$@.isoforms.fa)
	$(call cmp_expect,$@.isoform.read.map.txt)

test-transcriptome-nochecksplice: mkdirs ${TEST_ALIGN_BAM}
	flair transcriptome -b ${TEST_ALIGN_BAM} -g ${GENOME_FA} --junction_tab ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --no_check_splice
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

test-transcriptome-nostringent: mkdirs ${TEST_ALIGN_BAM}
	flair transcriptome -b ${TEST_ALIGN_BAM} -g ${GENOME_FA} --junction_tab ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --no_stringent
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

test-transcriptome-nostringentsplice: mkdirs ${TEST_ALIGN_BAM}
	flair transcriptome -b ${TEST_ALIGN_BAM} -g ${GENOME_FA} --junction_tab ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --no_stringent --no_check_splice
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

# --directRNA: trusts read strand and does not tolerate large indels in terminal exons
test-transcriptome-directrna: mkdirs ${TEST_ALIGN_BAM}
	flair transcriptome -b ${TEST_ALIGN_BAM} -g ${GENOME_FA} --junction_tab ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --directRNA
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

test-transcriptome-normends: mkdirs ${TEST_ALIGN_BAM}
	flair transcriptome -b ${TEST_ALIGN_BAM} -g ${GENOME_FA} --junction_tab ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --keep_intermediate --generate_map --normalize_ends
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

test-transcriptome-trust-ends: mkdirs ${TEST_ALIGN_BAM}
	flair transcriptome -b ${TEST_ALIGN_BAM} -g ${GENOME_FA} --junction_tab ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --trust_ends
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

# test of --filter comprehensive, which keeps the subset isoforms that the
# default nosubset filter drops from test-transcriptome-mutant-splice
test-transcriptome-comprehensive: mkdirs ${KRAS_MUT_SPLICE_BAM}
	flair transcriptome -b ${KRAS_MUT_SPLICE_BAM} -g ${GENOME_FA} \
	   --junction_bed=${KRAS_MUT_SPLICE_SJ_BED} \
	   -t 4 -f ${BASIC_ANNOTATION} -o $O/$@ --sample_name $@ --sjc_support 2 --filter comprehensive
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

# --filter comprehensive on the basic data, which has single-exon isoforms
# contained in spliced exons
test-transcriptome-comprehensive-basic: mkdirs ${TEST_ALIGN_BAM}
	flair transcriptome -b ${TEST_ALIGN_BAM} -g ${GENOME_FA} --junction_tab ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --filter comprehensive
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

test-transcriptome-mult-ends: mkdirs ${TEST_ALIGN_BAM}
	flair transcriptome  -b ${TEST_ALIGN_BAM} -g ${GENOME_FA} --junction_tab ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --max_ends 5 --keep_intermediate
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

test-transcriptome-mult-ends-trust-ends: mkdirs ${TEST_ALIGN_BAM}
	flair transcriptome  -b ${TEST_ALIGN_BAM} -g ${GENOME_FA} --junction_tab ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --max_ends 5 --trust_ends
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

test-transcriptome-partition: mkdirs ${TEST_ALIGN_BAM}
	flair transcriptome  -b ${TEST_ALIGN_BAM} -g ${GENOME_FA} --junction_tab ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --parallel_mode byregion
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

test-transcriptome-lowannot: mkdirs ${TEST_ALIGN_BAM}
	intron-prospector --genome-fasta=${GENOME_FA} --intron-bed6=$O/$@.IPjunctions.bed -C 0.0 ${TEST_ALIGN_BAM}
	flair transcriptome  -b ${TEST_ALIGN_BAM} -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) --junction_bed $O/$@.IPjunctions.bed --junction_support 2 -o $O/$@ --sample_name $@ --no_align_to_annot --sjc_support 3
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

test-transcriptome-noannot: mkdirs ${TEST_ALIGN_BAM}
	intron-prospector --genome-fasta=${GENOME_FA} --intron-bed6=$O/$@.IPjunctions.bed -C 0.0 ${TEST_ALIGN_BAM}
	flair transcriptome  -b ${TEST_ALIGN_BAM} -g ${GENOME_FA} -t 4 -o $O/$@ --sample_name $@ --junction_bed $O/$@.IPjunctions.bed --junction_support 2 --sjc_support 3
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

test-transcriptome-seg1: mkdirs  ${SEG1_ALIGN_BAM}
	flair transcriptome -b ${SEG1_ALIGN_BAM} -g ${GENOME_FA} -t 4 -f $(SEG1_ANNOTATION) -o $O/$@ --sample_name $@
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

test-transcriptome-seg1-badannot: mkdirs  ${SEG1_ALIGN_BAM}
	flair transcriptome -b ${SEG1_ALIGN_BAM} -g ${GENOME_FA} -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

test-transcriptome-seg1-badannot-multends: mkdirs  ${SEG1_ALIGN_BAM}
	flair transcriptome -b ${SEG1_ALIGN_BAM} -g ${GENOME_FA} -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --max_ends 4 --single_exon_support 2 --filter ginormous
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

# as seg1-badannot-multends, but comprehensive compares the overlapping single-exon
# isoforms by expression instead of keeping them all
test-transcriptome-seg1-badannot-comprehensive: mkdirs  ${SEG1_ALIGN_BAM}
	flair transcriptome -b ${SEG1_ALIGN_BAM} -g ${GENOME_FA} -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --max_ends 4 --single_exon_support 2 --filter comprehensive
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

# test of an artifical mutated splice site.
test-transcriptome-mutant-splice: mkdirs ${KRAS_MUT_SPLICE_BAM}
	flair transcriptome -b ${KRAS_MUT_SPLICE_BAM} -g ${GENOME_FA} \
	   --junction_bed=${KRAS_MUT_SPLICE_SJ_BED} \
	   -t 4 -f ${BASIC_ANNOTATION} -o $O/$@ --sample_name $@ --sjc_support 2
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

#test of an artificial read through transcript
test-transcriptome-readthrough: mkdirs ${DDX5_READTHROUGH_BAM}
	flair transcriptome -b ${DDX5_READTHROUGH_BAM} -g ${GENOME_FA} -f ${DDX5_READTHROUGH_ANNOT} -o $O/$@ --sample_name $@ --junction_bed=${DDX5_READTHROUGH_SJ} --sjc_support 2
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)


# the readthrough transcript annotated as its own gene: the readthrough isoform
# matches it, instead of being a novel isoform assigned to both DDX5 and POLG2
test-transcriptome-readthrough-gene: mkdirs ${DDX5_READTHROUGH_BAM} ${DDX5_READTHROUGH_GENE_ANNOT}
	flair transcriptome -b ${DDX5_READTHROUGH_BAM} -g ${GENOME_FA} -f ${DDX5_READTHROUGH_GENE_ANNOT} -o $O/$@ --sample_name $@ --junction_bed=${DDX5_READTHROUGH_SJ} --sjc_support 2
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

${DDX5_READTHROUGH_GENE_ANNOT}: ${DDX5_READTHROUGH_ANNOT} ${DDX5_READTHROUGH_GENE_GTF}
	@mkdir -p $O
	cat ${DDX5_READTHROUGH_ANNOT} ${DDX5_READTHROUGH_GENE_GTF} > $@.tmp
	mv -f $@.tmp $@

# test of --filter bysupport, which keeps a subset isoform with more than 1.2x
# the reads of its supersets; nosubset gives 2 isoforms, comprehensive 14
test-transcriptome-bysupport: mkdirs ${DDX5_BYSUPPORT_BAM}
	flair transcriptome -b ${DDX5_BYSUPPORT_BAM} -g ${GENOME_FA} -f ${DDX5_READTHROUGH_ANNOT} -o $O/$@ --sample_name $@ --junction_bed=${DDX5_READTHROUGH_SJ} --sjc_support 2 \
	   --filter bysupport
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

${DDX5_BYSUPPORT_BAM}: ${DDX5_READTHROUGH_BAM} ${DDX5_BYSUPPORT_DROP_READS}
	samtools view -h ${DDX5_READTHROUGH_BAM} | grep -v -F -w -f ${DDX5_BYSUPPORT_DROP_READS} | samtools view -b -o $@.tmp.bam -
	mv -f $@.tmp.bam $@
	samtools index $@

# basic reads with poly(A) trimmed, keeping single-exon isoforms.  Without poly(A),
# single-exon reads can only be stranded with --trust_strand, so compare these two
test-transcriptome-nopolya: mkdirs ${NOPOLYA_ALIGN_BAM}
	flair transcriptome -b ${NOPOLYA_ALIGN_BAM} -g ${GENOME_FA} --junction_tab ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ \
	   --filter ginormous --single_exon_support 2
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

test-transcriptome-nopolya-trust-strand: mkdirs ${NOPOLYA_ALIGN_BAM}
	flair transcriptome -b ${NOPOLYA_ALIGN_BAM} -g ${GENOME_FA} --junction_tab ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ \
	   --filter ginormous --single_exon_support 2 --trust_strand
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

# poly(T) heads on all reads, no annotation: the single-exon KRAS isoform gets re-assigned to the + strand
# all spliced reads strand as normal due to splice junction strand correction, only single exon reads get re-stranded
test-transcriptome-polyt: mkdirs ${POLYT_ALIGN_BAM}
	intron-prospector --genome-fasta=${GENOME_FA} --intron-bed6=$O/$@.IPjunctions.bed -C 0.0 ${POLYT_ALIGN_BAM}
	flair transcriptome -b ${POLYT_ALIGN_BAM} -g ${GENOME_FA} -t 4 -o $O/$@ --sample_name $@ --junction_bed $O/$@.IPjunctions.bed --junction_support 2 --sjc_support 3 \
	   --filter ginormous --single_exon_support 2
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

# novel DDX5 isoforms use the intron shared with the synthetic SHAREDJ1 gene; SHAREDJ1
# has fewer of their junctions, so gene assignment drops it and the results match
# test-transcriptome
test-transcriptome-shared-junc: mkdirs ${TEST_ALIGN_BAM} ${SHARED_JUNC_ANNOTATION}
	flair transcriptome -b ${TEST_ALIGN_BAM} -g ${GENOME_FA} --junction_tab ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f ${SHARED_JUNC_ANNOTATION} -o $O/$@ --sample_name $@
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

# the chr20 partition has one read, below --sjc_support 2, so it has no firstpass
# isoforms; results should be the chr12 and chr17 isoforms of test-transcriptome
test-transcriptome-no-firstpass: mkdirs ${ONECHR20_ALIGN_BAM}
	flair transcriptome -b ${ONECHR20_ALIGN_BAM} -g ${GENOME_FA} --junction_tab ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --sample_name $@ --generate_map --sjc_support 2
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)
	$(call cmp_expect,$@.isoform.read.map.txt)

# test-transcriptome-noannot with all loci on one chromosome, so KRAS and DDX5 (both
# minus strand) are separate novel genes on the same strand of one partition
test-transcriptome-onechrom: mkdirs ${ONECHROM_ALIGN_BAM}
	intron-prospector --genome-fasta=${ONECHROM_GENOME_FA} --intron-bed6=$O/$@.IPjunctions.bed -C 0.0 ${ONECHROM_ALIGN_BAM}
	flair transcriptome -b ${ONECHROM_ALIGN_BAM} -g ${ONECHROM_GENOME_FA} -t 4 -o $O/$@ --sample_name $@ --junction_bed $O/$@.IPjunctions.bed --junction_support 2 --sjc_support 3
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

${ONECHROM_ALIGN_BAM}: ${BASIC_READS_FA} ${ONECHROM_GENOME_FA}
	@mkdir -p $O
	flair align -r ${BASIC_READS_FA} --genome ${ONECHROM_GENOME_FA} -t 4 -o $(basename $@)

${KRAS_MUT_SPLICE_BAM}: ${KRAS_MUT_SPLICE_FQ}  ${GENOME_FA}
	@mkdir -p $O
	flair align -r ${KRAS_MUT_SPLICE_FQ} \
	   --genome ${GENOME_FA} -t 4 -o $(basename ${KRAS_MUT_SPLICE_BAM} .bam) \
	   --min_fragment_size 40 --max_intron_len 350k

#DDX5_readthrough.badread40x.fastq
${DDX5_READTHROUGH_BAM}: ${DDX5_READTHROUGH_FQ}  ${GENOME_FA}
	@mkdir -p $O
	flair align -r ${DDX5_READTHROUGH_FQ} \
	   --genome ${GENOME_FA} -t 4 -o $(basename ${DDX5_READTHROUGH_BAM} .bam) \
	   --min_fragment_size 40 --max_intron_len 350k

test-transcriptome-help: mkdirs
	flair transcriptome --help > $O/$@.out
	$(call cmp_expect,$@.out)


###
# flair variants
# FIXME need to add better checks for AAseq results, however the inconsistent sorting of the variants leads to inconsistent AAseq IDs which are hard to check
###

B1_GENOME_ALIGNED_BAM = output/B1_1.genomealigned.bam
B1_LONGSHOT = output/B1_1.genomealigned.longshot.vcf
ALIGN_LONGSHOT = output/test-align.longshot.vcf

${B1_GENOME_ALIGNED_BAM}: ${GENOME_FA} ${TEST_QUANTIFY_B1_BAM}
	@mkdir -p $O
	minimap2 -ax map-ont --secondary=no ${GENOME_FA} input/basic.reads.4.fa | samtools view -hb - | samtools sort - > ${B1_GENOME_ALIGNED_BAM}
	samtools index ${B1_GENOME_ALIGNED_BAM}

${B1_LONGSHOT}: ${GENOME_FA} ${TEST_QUANTIFY_B1_BAM} ${B1_GENOME_ALIGNED_BAM}
	@mkdir -p $O
	longshot --force_overwrite --bam ${B1_GENOME_ALIGNED_BAM} --ref ${GENOME_FA} --out ${B1_LONGSHOT} --min_cov 3 --min_alt_count 3 --strand_bias_pvalue_cutoff 0.000001

${ALIGN_LONGSHOT}: ${GENOME_FA} ${TEST_ALIGN_BAM}
	@mkdir -p $O
	longshot --force_overwrite --bam ${TEST_ALIGN_BAM} --ref ${GENOME_FA} --out ${ALIGN_LONGSHOT} --min_cov 3 --min_alt_count 3 --strand_bias_pvalue_cutoff 0.000001


# longshot is not currently include in the docker container, so tests needing
# it are skipped with a warning rather than failing
have_longshot := $(shell which longshot >/dev/null 2>&1 && echo yes || echo no)
longshot_tests = test-alleles test-alleles-withbed test-alleles-bedonly test-isoalleles test-variantquant

ifeq (${have_longshot},yes)
variants-tests: ${longshot_tests}
else
variants-tests:
	@echo "Note: skipping ${longshot_tests}: longshot is not in PATH, install it to run these tests" >&2
endif
variants-tests: test-identifyvars-snvs test-identifyvars-indels test-identifyvars-snvsindels

test-alleles: mkdirs ${GENOME_FA} ${TEST_ALIGN_BAM} ${ALIGN_LONGSHOT}
	flair alleles --tumor_bam ${TEST_ALIGN_BAM} --vcf ${ALIGN_LONGSHOT} -g ${GENOME_FA} -f $(BASIC_ANNOTATION) -o $O/$@ --generate_map --annotate_bams
	$(call cmp_expect,$@.allelegroups.counts.tsv)

test-alleles-withbed: mkdirs ${GENOME_FA} ${TEST_ALIGN_BAM} ${ALIGN_LONGSHOT}
	flair alleles --tumor_bam ${TEST_ALIGN_BAM} --vcf ${ALIGN_LONGSHOT} -g ${GENOME_FA} --isoform_bed $(TOME_ISOFORM_BED) -f $(BASIC_ANNOTATION) -o $O/$@ --generate_map --annotate_bams
	$(call cmp_expect,$@.allelegroups.counts.tsv)
	$(call cmp_expect,$@.allelegroups.vcf)
	$(call cmp_expect,$@.allelegroups.tumor.read.map.txt)

test-alleles-bedonly: mkdirs ${GENOME_FA} ${TEST_ALIGN_BAM} ${ALIGN_LONGSHOT}
	flair alleles --tumor_bam ${TEST_ALIGN_BAM} --vcf ${ALIGN_LONGSHOT} -g ${GENOME_FA} --isoform_bed $(TOME_ISOFORM_BED) -o $O/$@ --generate_map --annotate_bams
	$(call cmp_expect,$@.allelegroups.counts.tsv)

# test-isoalleles: mkdirs ${GENOME_FA} ${TEST_ALIGN_BAM} ${ALIGN_LONGSHOT}
# 	flair isoalleles --tumor_bam ${TEST_ALIGN_BAM} --vcf ${ALIGN_LONGSHOT} -g ${GENOME_FA}  --isoform_bed $(TOME_ISOFORM_BED) --read_map $(TOME_READ_MAP) -o $O/$@ --generate_map --annotate_bams
# 	$(call cmp_expect,$@.isoforms.allelegroups.bed)

test-isoalleles: mkdirs ${GENOME_FA} ${TEST_ALIGN_BAM} ${ALIGN_LONGSHOT}
	flair isoalleles --allele_vcf $E/test-alleles-withbed.allelegroups.vcf --allele_read_map $E/test-alleles-withbed.allelegroups.tumor.read.map.txt --genome ${GENOME_FA} \
	                 --isoform_bed $(TOME_ISOFORM_BED) --iso_read_map $(TOME_READ_MAP) -o $O/$@ --tumor_bam ${TEST_ALIGN_BAM} --generate_map
	$(call cmp_expect,$@.isoalleles.bed)
	$(call cmp_expect,$@.isoallele.counts.tsv)


test-variantquant: mkdirs ${GENOME_FA} ${TEST_ALIGN_BAM} ${ALIGN_LONGSHOT}
	flair variantquant --transcriptome_bam ${TEST_ALIGN_BAM} --vcf ${ALIGN_LONGSHOT} --isoform_bed $(TOME_ISOFORM_BED) -o $O/$@
	$(call cmp_expect_sort,$@.var.counts.tsv)

test-identifyvars-snvs: mkdirs ${GENOME_FA} ${TEST_ALIGN_BAM}
	python3 ../src/flair/identify_vars.py --genome ${GENOME_FA} --bam ${TEST_ALIGN_BAM} --identify_snvs yes --identify_indels no --output $O/$@
	$(call cmp_expect,$@.vcf)

test-identifyvars-indels: mkdirs ${GENOME_FA} ${TEST_ALIGN_BAM}
	python3 ../src/flair/identify_vars.py --genome ${GENOME_FA} --bam ${TEST_ALIGN_BAM} --identify_snvs no --identify_indels yes --output $O/$@
	$(call cmp_expect,$@.vcf)

test-identifyvars-snvsindels: mkdirs ${GENOME_FA} ${TEST_ALIGN_BAM}
	python3 ../src/flair/identify_vars.py --genome ${GENOME_FA} --bam ${TEST_ALIGN_BAM} --identify_snvs yes --identify_indels yes --output $O/$@
	$(call cmp_expect,$@.vcf)

###
# flair fusion
###
ERGIC3_KRAS_SIM_FA_GZ = $O/ERGIC3_KRAS_sim_fusion.badread20x.stranded.fastq.gz
FUSION_TEST_SHARED_OUT = $O/fusion_tests
FUSION_TEST_BAM = ${FUSION_TEST_SHARED_OUT}.bam

fusion-tests: test-fusion test-fusion-pary test-fusion-nongenic test-fusion-help

update-fusion-tests:
	for tname in "test-fusion" "test-fusion-pary" "test-fusion-nongenic"; \
	do \
	cp $O/$$tname.fusions.isoforms.bed $E/$$tname.fusions.isoforms.bed; \
	cp $O/$$tname.fusions.isoforms.fa $E/$$tname.fusions.isoforms.fa; \
	done

test-fusion: mkdirs ${GENOME_FA} ${ERGIC3_KRAS_SIM_FA_GZ} ${FUSION_TEST_BAM}
	flair fusion -g ${GENOME_FA} -f $(BASIC_ANNOTATION) -b ${FUSION_TEST_BAM} --max_loci 4 -o $O/$@ --sample_name $@ --min_support 2
	$(call cmp_expect_sort,$@.fusions.isoforms.bed)
	$(call cmp_expect_sort,$@.fusions.isoforms.fa)
# 	diff <(sort <(cut -f 1,2,3,5,6,7,8,9,10,11,12 $E/$@.fusions.isoforms.bed)) <(sort <(cut -f 1,2,3,5,6,7,8,9,10,11,12 $O/$@.fusions.isoforms.bed))
# 	diff <(sort <(grep -v ">" $E/$@.fusions.isoforms.fa)) <(sort <(grep -v ">" $O/$@.fusions.isoforms.fa))

test-fusion-pary: mkdirs ${GENOME_FA} ${ERGIC3_KRAS_SIM_FA_GZ} ${FUSION_TEST_BAM}
	flair fusion -g ${GENOME_FA} -f $I/parY.annotation.gtf -b ${FUSION_TEST_BAM} --max_loci 4 -o $O/$@ --sample_name $@ --min_support 2
	$(call cmp_expect_sort,$@.fusions.isoforms.bed)
	$(call cmp_expect_sort,$@.fusions.isoforms.fa)

test-fusion-nongenic: mkdirs ${GENOME_FA} ${ERGIC3_KRAS_SIM_FA_GZ} ${FUSION_TEST_BAM}
	flair fusion -g ${GENOME_FA} -f $I/basic.annotation.ERGIC3only.gtf -b ${FUSION_TEST_BAM} --max_loci 4 -o $O/$@ --sample_name $@ --min_support 2
	$(call cmp_expect_sort,$@.fusions.isoforms.bed)
	$(call cmp_expect_sort,$@.fusions.isoforms.fa)

${FUSION_TEST_BAM}: ${ERGIC3_KRAS_SIM_FA_GZ} ${GENOME_FA}
	@mkdir -p $O
	flair align -r $I/ERGIC3_KRAS_sim_fusion.badread20x.stranded.fastq \
	   --genome ${GENOME_FA} -t 4 -o ${FUSION_TEST_SHARED_OUT} \
	   --min_fragment_size 40 --max_intron_len 350k

test-fusion-help: mkdirs
	flair fusion  --help > $O/$@.out
	$(call cmp_expect,$@.out)

${ERGIC3_KRAS_SIM_FA_GZ}: $I/ERGIC3_KRAS_sim_fusion.badread20x.stranded.fastq
	@mkdir -p $(dir $@)
	gzip -c $< >$@


###
# flair spliceevents
###
spliceevents-tests: test-spliceevents test-spliceevents-extraannot

test-spliceevents: mkdirs ${TEST_ALIGN_BAM}
	python3 ../src/flair/flair_spliceevents.py -m $(SPLICEEVENTS_MANIFEST) -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@

test-spliceevents-extraannot: mkdirs ${TEST_ALIGN_BAM}
	python3 ../src/flair/flair_spliceevents.py -m $(SPLICEEVENTS_MANIFEST) -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) --annot_basic $(BASIC_ANNOTATION) -o $O/$@

# FIXME: orphan, this does not work and is not referenced; is this still needed
# the annotation file is there
test-refseq: mkdirs
	flair fusion -g genome.refSeq.fa -f $I/refSeq.test.annotation.gtf -r $I/ERGIC3_KRAS_sim_fusion.badread20x.stranded.fastq -b $O/$@_chimeric.bam --max_loci 4 -o $O/$@ --min_support 2

##
# flair_partition command
##
partition-tests: test-partition-basic test-partition-adjacent test-partition-min-reads \
	test-partition-sam

test-partition-basic: mkdirs
	flair_partition --threads=2 --log-level=info --bed=input/partition-test.bed output/$@.ranges.bed >& output/$@.out
	$(call cmp_expect,$@.ranges.bed)
	$(call cmp_expect,$@.out)

# merge adjacent partitions
test-partition-adjacent: mkdirs
	flair_partition --threads=2 --part_merge_dist=10000 --log-level=info --bed=input/partition-test.bed output/$@.ranges.bed >& output/$@.out
	$(call cmp_expect,$@.ranges.bed)
	$(call cmp_expect,$@.out)

# minimum reads per partition
test-partition-min-reads: mkdirs
	flair_partition --threads=2 --min_partition_items=3  --log-level=info --bed=input/partition-test.bed output/$@.ranges.bed >& output/$@.out
	$(call cmp_expect,$@.ranges.bed)
	$(call cmp_expect,$@.out)

# with bed, sam, and gtf
test-partition-sam: mkdirs
	flair_partition --threads=2 --log-level=info --bed=input/partition-test.bed --bam=input/tiny.sam --gtf=input/seg1.gencodeV47.gtf \
	    output/$@.ranges.bed >& output/$@.out
	$(call cmp_expect,$@.ranges.bed)
	$(call cmp_expect,$@.out)

##
# library tests with pytest
##
lib-tests:
	${PYTEST_CMD} .

##
# build test input files
##
${GENOME_FA}:
	mkdir -p $(dir $@)
	rm -f $@.fai
	wget -nv https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr12.fa.gz
	wget -nv https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr17.fa.gz
	wget -nv https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz
	zcat chr12.fa.gz chr17.fa.gz chr20.fa.gz > $@
	rm chr12.fa.gz chr17.fa.gz chr20.fa.gz
	samtools faidx $@

# one sequence, re-wrapped at 50 bases so it can be indexed
${ONECHROM_GENOME_FA}: ${GENOME_FA}
	awk 'BEGIN {print ">chr12_17_20"} /^>/ {next} {buf = buf $$0; while (length(buf) >= 50) {print substr(buf, 1, 50); buf = substr(buf, 51)}} END {if (buf != "") print buf}' ${GENOME_FA} > $@.tmp
	mv -f $@.tmp $@
	samtools faidx $@

${TOME_OUT_FA}: | mkdirs
	cp $E/test-transcriptome.isoforms.fa $@

${TOME_OUT_FAI}: ${TOME_OUT_FA}
	samtools faidx --output $@ $<

$O/%.genomealigned.bam: $I/%.fa ${GENOME_FA} | mkdirs
	minimap2 -ax splice --secondary=no ${GENOME_FA} $< | samtools view -hb - | samtools sort - > $@.tmp.bam
	mv -f $@.tmp.bam $@
	samtools index $@

##
# test hub (work in progress)
# 
#  make hub-check huburl=https://hgwdev.gi.ucsc.edu/~markd/flair/flair-mddev/test/hub/hub.txt
##
HUB_OUT_DIR = $O/hg38

HG38_CHROMS = https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromAlias.bb

TEST_TRANSCRIPOME_MUTANT_SPLICE_BB = output/test-transcripome-mutant-splice.isoforms.bb
TEST_TRANSCRIPTOME_SEG1_BB =  output/test-transcriptome-seg1.isoforms.bb

test-hub:  ${KRAS_MUT_SPLICE_BB} ${KRAS_MUT_SPLICE_SJ_BB} ${TEST_TRANSCRIPOME_MUTANT_SPLICE_BB} \
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
	rm -rf $O $D __pycache__
	rm -rf [0-9a-f]*-[0-9a-f]*-[0-9a-f]*-[0-9a-f]*-[0-9a-f]*

real-clean : clean
	rm -rf ${TI}
