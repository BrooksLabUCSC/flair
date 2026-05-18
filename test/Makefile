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
		combine-tests predict-productivity-tests diff-iso-usage-tests \
		plot-usage-tests transcriptome-tests fusion-tests variants-tests \
		spliceevents-tests partition-tests lib-tests

test-no-align:
	${MAKE} quantify-tests bed-to-gtf-tests \
		combine-tests predict-productivity-tests diff-iso-usage-tests \
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
		test-predict-productivity-help test-diff-iso-usage-help \
		test-plot-usage-help test-transcriptome-help test-variants-help test-fusion-help

test-fast:
	${MAKE} test-combine-bed-only \
		test-predict-productivity test-diff-iso-usage test-plot-usage test-transcriptome \
		test-fusion test-variants



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
BASIC_DIFFSPLICE_BED = $I/test-collapse-annot.isoforms.bed

# FIXME: recreate these files with transcriptome, and others below 
BASIC_ISOFORMS_FA = $E/test-collapse-annot.isoforms.fa
BASIC_ISOFORMS_BED = $E/test-collapse-annot.isoforms.bed

# copy so *.fai does end up in expect directory
BASIC_ISOFORMS_OUT_FA = $O/test-collapse-annot.isoforms.fa
BASIC_ISOFORMS_OUT_FAI = $O/test-collapse-annot.isoforms.fa.fai


##
# Files generate by tests that are used by other tests.
# These enfoce files being create in the right order.
## 
TEST_ALIGN_BAM = $O/test-align.bam
TEST_QUANTIFY_B1_BAM = $O/test-quantify.B1.B.flair.aligned.bam

###
# seg1: single exon tests
###
SEG1_READS_FQ = input/seg1.cdna-ont.fastq
SEG1_ANNOTATION = input/seg1.gencodeV47.gtf
SEG1_PROMOTER_REGIONS = input/seg1.promoter-regions.bed

SEG1_ALIGN_BED = $E/test-align-seg1.bed
SEG1_ALIGN_BAM = $O/test-align-seg1.bam
SEG1_ISOFORMS_FA = expected/test-collapse-seg1.isoforms.fa
SEG1_ISOFORMS_BED = expected/test-collapse-seg1.isoforms.bed

##
# mutant splice site, BED is truth
##
DDX5_READTHROUGH_FQ = input/DDX5_readthrough.badread40x.fastq
DDX5_READTHROUGH_BAM = output/DDX5_readthrough_sim_badread40x.bam
DDX5_READTHROUGH_ANNOT = input/DDX5_POLG2_gencode38.gtf
DDX5_READTHROUGH_SJ = input/DDX5_readthrough_IPSJ.bed

KRAS_MUT_SPLICE_FQ = input/KRAS_sim_genomic_mut_driven_splicing.badreadsimPB.fastq
KRAS_MUT_SPLICE_BAM = output/KRAS_sim_genomic_mut_driven_splicing.badreadsimPB.bam
KRAS_MUT_SPLICE_SJ_BED = input/KRAS_sim_genomic_mut_driven_splicing.badreadsimPB.sj.bed
KRAS_MUT_SPLICE_SJ_BB = output/KRAS_sim_genomic_mut_driven_splicing.badreadsimPB.sj.bb
KRAS_MUT_SPLICE_BED = input/KRAS_sim_genomic_mut_driven_splicing.bed
KRAS_MUT_SPLICE_BB = output/KRAS_sim_genomic_mut_driven_splicing.bb

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

align-tests: test-align test-align-seg1 test-align-help

# just checks bed4 to make sure we are getting alignments
test-align : mkdirs ${GENOME_FA}
	flair align -r $(BASIC_READS_FA) --genome ${GENOME_FA} -t 4 -o $O/$@
	bedtools bamtobed -i $O/$@.bam > $O/$@.bed
	$(call cmp_expect,$@.bed)

${TEST_ALIGN_BAM}: test-align
${SEG1_ALIGN_BAM}: test-align-seg1

test-align-seg1 : mkdirs ${GENOME_FA}
	flair align -r ${SEG1_READS_FQ} --genome ${GENOME_FA} -t 4 -o $O/$@
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
	$(call cmp_expect,$@/combined.bed)
	$(call cmp_expect,$@/combined.counts.tsv)
	$(call cmp_expect,$@/combined.isoform.map.txt)
endef
define combine_fa_diff
	$(call combine_diff)
	$(call cmp_expect,$@/combined.fa)
endef

combine-tests: test-combine-bed-only test-combine-bed-fa test-combine-bed-fa-map

update-combine-tests:
	for tname in "test-combine-bed-only" "test-combine-bed-fa" "test-combine-bed-fa-map"; \
	do \
	cp $O/$$tname/combined.bed $E/$$tname/combined.bed; \
	cp $O/$$tname/combined.counts.tsv $E/$$tname/combined.counts.tsv; \
	cp $O/$$tname/combined.isoform.map.txt $E/$$tname/combined.isoform.map.txt; \
	done
	for tname in "test-combine-bed-fa" "test-combine-bed-fa-map"; \
	do \
	cp $O/$$tname/combined.fa $E/$$tname/combined.fa; \
	done

test-combine-bed-only: mkdirs
	rm -rf $O/$@
	mkdir -p $O/$@
	flair combine --manifest=input/bed-only.combine_manifest.txt -o $O/$@/combined
	$(call combine_diff)

test-combine-bed-fa: mkdirs
	rm -rf $O/$@
	mkdir -p $O/$@
	flair combine --manifest=input/bed-fa.combine_manifest.txt -o $O/$@/combined
	$(call combine_fa_diff)

test-combine-bed-fa-map: mkdirs
	rm -rf $O/$@
	mkdir -p $O/$@
	flair combine --manifest=input/bed-fa-map.combine_manifest.txt -o $O/$@/combined
	$(call combine_fa_diff)

###
# FLAIR QUANTIFY
###
quantify-tests: test-quantify test-quantify-nobam test-quantify-help

test-quantify: mkdirs ${BASIC_ISOFORMS_OUT_FAI}
	rm -rf $O/$@.tmp
	flair quantify -r $(BASIC_READS_MANIFEST) -i $(BASIC_ISOFORMS_OUT_FA) --isoform_bed $(BASIC_ISOFORMS_BED) --generate_map --temp_dir $O/$@.tmp --tpm --sample_id_only -o $O/$@ --output_bam
	samtools view $O/$@.A1.A.flair.aligned.bam | cut -f 1,3 > $O/$@.A1.A.flair.aligned.readtotranscript.txt
	$(call cmp_expect,$@.counts.tsv)
	$(call cmp_expect_sort,$@.A1.A.flair.aligned.readtotranscript.txt)
	$(call cmp_expect,$@.tpm.tsv)

${TEST_QUANTIFY_B1_BAM}: test-quantify

# without map or bam outputs
test-quantify-nobam : mkdirs ${BASIC_ISOFORMS_OUT_FAI}
	rm -rf $O/$@.tmp
	flair quantify -r $(BASIC_READS_MANIFEST) -i $(BASIC_ISOFORMS_OUT_FA) --isoform_bed $(BASIC_ISOFORMS_BED) --temp_dir $O/$@.tmp --tpm --sample_id_only -o $O/$@
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
	flair diffexp -q $(BASIC_COUNTS_MATRIX) -o $O/$@ -e 1 -of
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
	flair diffsplice -i $(BASIC_DIFFSPLICE_BED) -q $(BASIC_DIFFSPLICE_MATRIX) --test -o $O/$@ -of
	$(call cmp_expect_sort,$@/diffsplice.alt3.events.quant.tsv)
	$(call cmp_expect_sort,$@/diffsplice.alt5.events.quant.tsv)
	$(call cmp_expect,$@/diffsplice.es.events.quant.tsv)
	$(call cmp_expect_sort,$@/diffsplice.ir.events.quant.tsv)
	$(call cmp_expect,$@/drimseq_alt3_A_v_B.tsv)
	$(call cmp_expect,$@/drimseq_es_A_v_B.tsv)

test-diffsplice-cond: mkdirs
	flair diffsplice -i $(BASIC_ISOFORMS_BED) -q $(BASIC_DIFFSPLICE_MATRIX) --test -o $O/$@ -of --conditionA A --conditionB B
	$(call cmp_expect_sort,test-diffsplice/diffsplice.alt3.events.quant.tsv)
	$(call cmp_expect_sort,test-diffsplice/diffsplice.alt5.events.quant.tsv)
	$(call cmp_expect,test-diffsplice/diffsplice.es.events.quant.tsv)
	$(call cmp_expect_sort,test-diffsplice/diffsplice.ir.events.quant.tsv)
	$(call cmp_expect,test-diffsplice/drimseq_alt3_A_v_B.tsv)
	$(call cmp_expect,test-diffsplice/drimseq_es_A_v_B.tsv)

test-diffsplice-help: mkdirs
	flair diffsplice --help > $O/$@.out
	$(call cmp_expect,$@.out)


bed-to-gtf-tests: test-bed-to-gtf test-bed-to-gtf-cds test-bed-to-gtf-extracols

test-bed-to-gtf : mkdirs
	bed_to_gtf $(BASIC_ISOFORMS_BED) > $O/$@.gtf
	$(call cmp_expect,$@.gtf)

test-bed-to-gtf-cds : 
	bed_to_gtf $E/test-predict-productivity.bed > $O/$@.gtf
	$(call cmp_expect,$@.gtf)

test-bed-to-gtf-extracols : 
	bed_to_gtf --gene_col 0 --as_file $E/test-transcriptome.isoforms.as --extra_cols_to_use ref_transcript_id,gene_desc,productivity $E/test-transcriptome.isoforms.bed > $O/$@.gtf
	$(call cmp_expect,$@.gtf)

##
# PREDICT PRODUCTIVITY
##
predict-productivity-tests: test-predict-productivity test-predict-productivity-help

test-predict-productivity : mkdirs ${GENOME_FA}
	predictProductivity -i $(BASIC_ISOFORMS_BED) -g $(BASIC_ANNOTATION) -f ${GENOME_FA} -o $O/$@
	$(call cmp_expect,$@.bed)
	$(call cmp_expect,$@.info.tsv)

test-predict-productivity-help: mkdirs
	predictProductivity --help > $O/$@.out
	$(call cmp_expect,$@.out)

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

##
# plot isoform usage
##
plot-usage-tests: test-plot-usage test-plot-usage-help

test-plot-usage: mkdirs
	plot_isoform_usage expected/test-predict-productivity.bed expected/test-collapse-annot.annotated_transcripts.alignment.counts ENSG00000133703.12 -o $O/$@
	file $O/test-plot-usage_isoforms.png | grep 'PNG image data' >/dev/null
	file $O/test-plot-usage_usage.png | grep 'PNG image data' >/dev/null

test-plot-usage-help: mkdirs
	plot_isoform_usage --help > $O/$@.out
	$(call cmp_expect,$@.out)

###
# flair transcriptome
###

transcriptome-tests: test-transcriptome test-transcriptome-partition test-transcriptome-lowannot test-transcriptome-noannot \
                    test-transcriptome-seg1 test-transcripome-mutant-splice test-transcriptome-readthrough test-transcriptome-mult-ends \
					test-transcriptome-trust-ends test-transcriptome-mult-ends-trust-ends test-transcriptome-help

update-transcriptome-tests:
	for tname in "test-transcriptome" "test-transcriptome-partition" "test-transcriptome-lowannot" "test-transcriptome-noannot" \
	"test-transcriptome-seg1" "test-transcripome-mutant-splice" "test-transcriptome-readthrough" "test-transcriptome-mult-ends" \
	"test-transcriptome-trust-ends" "test-transcriptome-mult-ends-trust-ends"; \
	do \
	cp $O/$$tname.isoforms.bed $E/$$tname.isoforms.bed; \
	cp $O/$$tname.isoforms.gtf $E/$$tname.isoforms.gtf; \
	done
	cp $O/test-transcriptome.isoforms.as $E/test-transcriptome.isoforms.as


test-transcriptome: mkdirs ${TEST_ALIGN_BAM}
	flair transcriptome -b ${TEST_ALIGN_BAM} -g ${GENOME_FA} --junction_tab ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --keep_intermediate
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

test-transcriptome-trust-ends: mkdirs ${TEST_ALIGN_BAM}
	flair transcriptome -b ${TEST_ALIGN_BAM} -g ${GENOME_FA} --junction_tab ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --trust_ends
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
# 	$(diff_bed_gtf)

test-transcriptome-mult-ends: mkdirs ${TEST_ALIGN_BAM}
	flair transcriptome  -b ${TEST_ALIGN_BAM} -g ${GENOME_FA} --junction_tab ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --max_ends 5 --keep_intermediate
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

test-transcriptome-mult-ends-trust-ends: mkdirs ${TEST_ALIGN_BAM}
	flair transcriptome  -b ${TEST_ALIGN_BAM} -g ${GENOME_FA} --junction_tab ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --max_ends 5 --trust_ends
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
# 	$(diff_bed_gtf)

test-transcriptome-partition: mkdirs ${TEST_ALIGN_BAM}
	flair transcriptome  -b ${TEST_ALIGN_BAM} -g ${GENOME_FA} --junction_tab ${BASIC_SHORTREAD_JUNCTIONS} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --parallel_mode byregion
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

test-transcriptome-lowannot: mkdirs ${TEST_ALIGN_BAM}
	intronProspector --genome-fasta=${GENOME_FA} --intron-bed6=$O/$@.IPjunctions.bed -C 0.0 ${TEST_ALIGN_BAM}
	flair transcriptome  -b ${TEST_ALIGN_BAM} -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) --junction_bed $O/$@.IPjunctions.bed --junction_support 2 -o $O/$@ --no_align_to_annot --sjc_support 3
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

test-transcriptome-noannot: mkdirs ${TEST_ALIGN_BAM}
	intronProspector --genome-fasta=${GENOME_FA} --intron-bed6=$O/$@.IPjunctions.bed -C 0.0 ${TEST_ALIGN_BAM}
	flair transcriptome  -b ${TEST_ALIGN_BAM} -g ${GENOME_FA} -t 4 -o $O/$@ --junction_bed $O/$@.IPjunctions.bed --junction_support 2 --sjc_support 3
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

test-transcriptome-seg1: mkdirs  ${SEG1_ALIGN_BAM}
	flair transcriptome -b ${SEG1_ALIGN_BAM} -g ${GENOME_FA} -t 4 -f $(SEG1_ANNOTATION) -o $O/$@
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

# test of an artifical mutated splice site.
test-transcripome-mutant-splice: mkdirs ${KRAS_MUT_SPLICE_BAM}
	flair transcriptome -b ${KRAS_MUT_SPLICE_BAM} -g ${GENOME_FA} \
	   --junction_bed=${KRAS_MUT_SPLICE_SJ_BED} \
	   -t 4 -f ${BASIC_ANNOTATION} -o $O/$@ --sjc_support 2
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)

#test of an artificial read through transcript
test-transcriptome-readthrough: mkdirs ${DDX5_READTHROUGH_BAM}
	flair transcriptome -b ${DDX5_READTHROUGH_BAM} -g ${GENOME_FA} -f ${DDX5_READTHROUGH_ANNOT} -o $O/$@ --junction_bed=${DDX5_READTHROUGH_SJ} --sjc_support 2
# 	$(call cmp_expect_bed_simple,$@.isoforms.bed)
	$(diff_bed_gtf)


${KRAS_MUT_SPLICE_BAM}: ${KRAS_MUT_SPLICE_FQ}  ${GENOME_FA}
	@mkdir -p $O
	flair align -r ${KRAS_MUT_SPLICE_FQ} \
	   --genome ${GENOME_FA} -t 4 -o $(basename ${KRAS_MUT_SPLICE_BAM} .bam) \
	   --minfragmentsize 40 --maxintronlen 350k

#DDX5_readthrough.badread40x.fastq
${DDX5_READTHROUGH_BAM}: ${DDX5_READTHROUGH_FQ}  ${GENOME_FA}
	@mkdir -p $O
	flair align -r ${DDX5_READTHROUGH_FQ} \
	   --genome ${GENOME_FA} -t 4 -o $(basename ${DDX5_READTHROUGH_BAM} .bam) \
	   --minfragmentsize 40 --maxintronlen 350k

test-transcriptome-help: mkdirs
	flair transcriptome --help > $O/$@.out
	$(call cmp_expect,$@.out)


###
# flair variants
# FIXME need to add better checks for AAseq results, however the inconsistent sorting of the variants leads to inconsistent AAseq IDs which are hard to check
###
variants-tests:  test-variants test-variants-help

# longshot is not currently include in the docker container
have_longshot := $(shell which longshot >/dev/null 2>&1 && echo yes || echo no)

ifeq (${have_longshot},yes)
test-variants: mkdirs ${GENOME_FA} ${TEST_QUANTIFY_B1_BAM} ${BASIC_ISOFORMS_OUT_FAI}
	minimap2 -ax map-ont --secondary=no ${GENOME_FA} input/basic.reads.4.fa | samtools view -hb - | samtools sort - > $O/$@.genomealigned.bam
	samtools index $O/$@.genomealigned.bam
	longshot --force_overwrite --bam $O/$@.genomealigned.bam --ref ${GENOME_FA} --out $O/$@.genomealigned.longshot.vcf --min_cov 3 --min_alt_count 3 --strand_bias_pvalue_cutoff 0.000001
	echo 's1    ${TEST_QUANTIFY_B1_BAM}   $O/$@.genomealigned.longshot.vcf' > $O/$@.varmanifest.tsv
	flair variants -m $O/$@.varmanifest.tsv -i $(BASIC_ISOFORMS_OUT_FA) -b $(BASIC_ISOFORMS_BED) -g ${GENOME_FA} -f $(BASIC_ANNOTATION) -o $O/$@
	$(call cmp_expect,$@.isoswithvars.counts.tsv)
	$(call cmp_expect,$@.aaseq.counts.tsv)
else
test-variants:
	@echo 'WARNING: ***' >&2
	@echo 'WARNING: *** longshot command not found, test-variants skipped' >&2
	@echo 'WARNING: ***' >&2
endif

test-variants-help: mkdirs
	flair variants --help > $O/$@.out
	$(call cmp_expect,$@.out)


###
# flair fusion
###
ERGIC3_KRAS_SIM_FA_GZ = $O/ERGIC3_KRAS_sim_fusion.badread20x.stranded.fastq.gz
FUSION_TEST_SHARED_OUT = $O/fusion_tests
FUSION_TEST_BAM = ${FUSION_TEST_SHARED_OUT}.bam

fusion-tests: test-fusion test-fusion-pary test-fusion-help

test-fusion: mkdirs ${GENOME_FA} ${ERGIC3_KRAS_SIM_FA_GZ} ${FUSION_TEST_BAM}
	flair fusion -g ${GENOME_FA} -f $(BASIC_ANNOTATION) -b ${FUSION_TEST_BAM} --maxloci 4 -o $O/$@ -s 2
	diff <(sort <(cut -f 1,2,3,5,6,7,8,9,10,11,12 $E/$@.fusions.isoforms.bed)) <(sort <(cut -f 1,2,3,5,6,7,8,9,10,11,12 $O/$@.fusions.isoforms.bed))
	diff <(sort <(grep -v ">" $E/$@.fusions.isoforms.fa)) <(sort <(grep -v ">" $O/$@.fusions.isoforms.fa))

test-fusion-pary: mkdirs ${GENOME_FA} ${ERGIC3_KRAS_SIM_FA_GZ} ${FUSION_TEST_BAM}
	flair fusion -g ${GENOME_FA} -f $I/parY.annotation.gtf -b ${FUSION_TEST_BAM} --maxloci 4 -o $O/$@ -s 2
	$(call cmp_expect_sort,$@.fusions.isoforms.bed)
	$(call cmp_expect_sort,$@.fusions.isoforms.fa)

${FUSION_TEST_BAM}: ${ERGIC3_KRAS_SIM_FA_GZ} ${GENOME_FA}
	@mkdir -p $O
	flair align -r $I/ERGIC3_KRAS_sim_fusion.badread20x.stranded.fastq \
	   --genome ${GENOME_FA} -t 4 -o ${FUSION_TEST_SHARED_OUT} \
	   --minfragmentsize 40 --maxintronlen 350k

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
	python3 ../src/flair/flair_spliceevents.py -m $I/spliceevents.manifest.tsv -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@

test-spliceevents-extraannot: mkdirs ${TEST_ALIGN_BAM}
	python3 ../src/flair/flair_spliceevents.py -m $I/spliceevents.manifest.tsv -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) --annot_basic $(BASIC_ANNOTATION) -o $O/$@

# FIXME: orphan, this does not work and is not referenced; is this still needed
# the annotation file is there
test-refseq: mkdirs
	flair fusion -g genome.refSeq.fa -f $I/refSeq.test.annotation.gtf -r $I/ERGIC3_KRAS_sim_fusion.badread20x.stranded.fastq -b $O/$@_chimeric.bam --maxloci 4 -o $O/$@ -s 2

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
	flair_partition --threads=2 -part_merge_dist=10000 --log-level=info --bed=input/partition-test.bed output/$@.ranges.bed >& output/$@.out
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

${BASIC_ISOFORMS_OUT_FA}: ${BASIC_ISOFORMS_FA}
	mkdir -p $(dir $@)
	cp -f $< $@

${BASIC_ISOFORMS_OUT_FAI}: ${BASIC_ISOFORMS_OUT_FA}
	samtools faidx --output $@ $<

##
# test hub (work in progress)
##
HUB_OUT_DIR = $O/hg38

HG38_CHROMS = https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromAlias.bb

TEST_TRANSCRIPOME_MUTANT_SPLICE_BB = output/test-transcripome-mutant-splice.isoforms.bb

test-hub:  ${KRAS_MUT_SPLICE_BB} ${KRAS_MUT_SPLICE_SJ_BB} ${TEST_TRANSCRIPOME_MUTANT_SPLICE_BB}

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
