root = ..
include ${root}/defs.mk

# If you want to use the install flair rather than the one in this
# tree use:
#   make test use_installed_flair=yes
ifneq (${use_installed_flair},yes)
    export PATH := ${root}/bin:${PATH}
    export PYTHONPATH := ${root}/src:${PYTHONPATH}
endif

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
BASIC_COUNTS_MATRIX = $I/basic.counts_matrix_diffexp.tsv
BASIC_DIFFSPLICE_MATRIX = $I/basic.counts_matrix_diffsplice.tsv

BASIC_ALIGN_BED = $E/test-align.bed
BASIC_CORRECT_READS_BED = $E/test-correct_all_corrected.bed
BASIC_ISOFORMS_FA = $E/test-collapse-annot.isoforms.fa
BASIC_ISOFORMS_BED = $E/test-collapse-annot.isoforms.bed


###
# seg1: single exon tests
###
SEG1_READS_FQ = input/seg1.cdna-ont.fastq
SEG1_ANNOTATION = input/seg1.gencodeV47.gtf
SEG1_PROMOTER_REGIONS = input/seg1.promoter-regions.bed

SEG1_ALIGN_BED = $E/test-align-seg1.bed
SEG1_CORRECT_READS_BED = $E/test-correct-seg1_all_corrected.bed
SEG1_ISOFORMS_FA = expected/test-collapse-seg1.isoforms.fa
SEG1_ISOFORMS_BED = expected/test-collapse-seg1.isoforms.bed


# The help test results vary on screen width.  This forces them to be consistent.
# It might seem weird to see if help works, but some it did break in the past.
export COLUMNS = 80

##
# compare output, expected with result
##
define cmp_expect
	diff $E/$(1) $O/$(1)
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

default:
	@echo To run all tests, use: >&2
	@echo "   " make test >&2
	@echo To run expression diff tests that depend on R >&2
	@echo "   " make test-diffexpress >&2
	@echo To run all tests >&2
	@echo "   " make test-all >&2
	@echo otherwise run individual tests by specifying test target, such as: >&2
	@echo "   " make align-tests >&2
	@echo remove output files: >&2
	@echo "   " make clean >&2
	@echo remove output and downloaded genomes files: >&2
	@echo "   " make real-clean >&2
	@exit 1

test: main-tests align-tests correct-tests collapse-tests quantify-tests \
	combine-tests predict-productivity-tests diff-iso-usage-tests \
	plot-usage-tests other-progs-tests

test-diffexpress: diffexp-tests diffsplice-tests

test-all: test test-diffexpress

test-help: test-main-help test-align-help test-correct-help test-collapse-help \
           test-quantify-help test-diffexp-help test-diffsplice-help \
           test-predict-productivity-help test-diff-iso-usage-help

###
# MAIN COMMAND
###
main-tests: test-main-help

test-main-help: mkdirs
	flair --help 2>&1 | fgrep -v Version> $O/$@.out
	$(call cmp_expect,$@.out)

###
# FLAIR ALIGN
###

align-tests: test-align test-align-intpriming test-align-seg1 test-align-help

test-align : mkdirs ${GENOME_FA}
	flair align -r $(BASIC_READS_FA) --genome ${GENOME_FA} -t 4 -o $O/$@ --filtertype keepsup
	$(call cmp_expect,$@.bed)

test-align-intpriming : mkdirs ${GENOME_FA}
	flair align -r $(BASIC_READS_FA) --genome ${GENOME_FA} -t 4 -o $O/$@ --remove_internal_priming -f $(BASIC_ANNOTATION) --remove_singleexon --filtertype keepsup
	$(call cmp_expect,$@.bed)

test-align-seg1 : mkdirs ${GENOME_FA}
	flair align -r ${SEG1_READS_FQ} --genome ${GENOME_FA} -t 4 -o $O/$@ --filtertype keepsup
	$(call cmp_expect,$@.bed)

test-align-help: mkdirs
	flair align --help >&$O/$@.out
	$(call cmp_expect,$@.out)

###
# FLAIR CORRECT
###
correct-tests: test-correct test-correct-seg1 test-correct-help

test-correct: mkdirs ${GENOME_FA}
	flair correct -q ${BASIC_ALIGN_BED} -j $(BASIC_SHORTREAD_JUNCTIONS) -f $(BASIC_ANNOTATION_INCOMPLETE) -g ${GENOME_FA} -o $O/$@
	$(call cmp_expect_sort,$@_all_inconsistent.bed)
	$(call cmp_expect_sort,$@_all_corrected.bed)

test-correct-seg1: mkdirs ${GENOME_FA}
	flair correct -q ${SEG1_ALIGN_BED} -g ${GENOME_FA} -f ${SEG1_ANNOTATION} -o $O/$@
	$(call cmp_expect_sort,$@_all_inconsistent.bed)
	$(call cmp_expect_sort,$@_all_corrected.bed)

test-correct-help: mkdirs
	flair correct --help >&$O/$@.out
	$(call cmp_expect,$@.out)

###
# FLAIR COLLAPSE
###
collapse-tests: test-collapse test-collapse-annot \
	test-collapse-support test-collapse-seg1 \
	test-collapse-help

# macros to do comparions with expected results
define collapse_diff
	$(call cmp_expect_sort,$@.isoforms.bed)
	$(call cmp_expect_sort,$@.isoforms.fa)
	$(call cmp_expect_sort,$@.isoforms.gtf)
endef

define collapse_promoter_diff
	$(call cmp_expect_sort,$@.promoter_supported.bed)
endef

# resulting read map differs between run (see #443)
# just check if it exists for now
# should do: $(call cmp_expect_sort,$@.isoform.read.map.txt)
define collapse_diff_read_map
	@test -e $O/${1} || (echo Error: output file missing: $O/${1}; exit 1)
endef

# basic collapse
test-collapse: mkdirs ${GENOME_FA}
	rm -rf $O/$@.tmp
	flair collapse -r $(BASIC_READS_FA) -q $(BASIC_CORRECT_READS_BED) -g ${GENOME_FA} -t 4 --generate_map --temp_dir $O/$@.tmp --keep_intermediate -f $(BASIC_ANNOTATION) -o $O/$@ -p $(BASIC_PROMOTER_REGIONS)
	$(collapse_diff)
	$(collapse_diff_read_map)
	$(collapse_promoter_diff)

# collapse annotation-reliant
##these reads are low quality, many isoforms fail with the check_splice parameter
test-collapse-annot: mkdirs ${GENOME_FA}
	rm -rf $O/$@.tmp
	flair collapse -r $(BASIC_READS_FA) -q $(BASIC_CORRECT_READS_BED) -g ${GENOME_FA} -t 4 --generate_map --temp_dir $O/$@.tmp --keep_intermediate -f $(BASIC_ANNOTATION) -o $O/$@ -p $(BASIC_PROMOTER_REGIONS) --annotation_reliant generate --stringent
	$(collapse_diff)
	$(collapse_diff_read_map)
	$(collapse_promoter_diff)

test-collapse-seg1: mkdirs ${GENOME_FA}
	rm -rf $O/$@.tmp
	flair collapse -r ${SEG1_READS_FQ} -q ${SEG1_CORRECT_READS_BED} -g ${GENOME_FA} -t 4 --generate_map --temp_dir $O/$@.tmp --keep_intermediate -f $(SEG1_ANNOTATION) -o $O/$@ -p ${SEG1_PROMOTER_REGIONS} --stringent
	$(collapse_diff)
	$(collapse_diff_read_map)
	$(collapse_promoter_diff)

# collapse with support
test-collapse-support: mkdirs ${GENOME_FA}
	rm -rf $O/$@.tmp
	flair collapse -s 0.25 -r $(BASIC_READS_FA) -q $(BASIC_CORRECT_READS_BED) -g ${GENOME_FA} -t 4 --generate_map --temp_dir $O/$@.tmp --keep_intermediate -f $(BASIC_ANNOTATION) -o $O/$@ -p $(BASIC_PROMOTER_REGIONS)
	$(collapse_diff)
	$(collapse_diff_read_map)
	$(collapse_promoter_diff)

test-collapse-help: mkdirs
	flair collapse --help >&$O/$@.out
	$(call cmp_expect,$@.out)


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

test-quantify: mkdirs
	rm -rf $O/$@.tmp
	flair quantify -r $(BASIC_READS_MANIFEST) -i $(BASIC_ISOFORMS_FA) --isoform_bed $(BASIC_ISOFORMS_BED) --generate_map --temp_dir $O/$@.tmp --tpm --sample_id_only -o $O/$@ --output_bam
	samtools view $O/$@.A1.A.flair.aligned.bam | cut -f 1,3 > $O/$@.A1.A.flair.aligned.readtotranscript.txt
	$(call cmp_expect_sort,$@.A1.A.flair.aligned.readtotranscript.txt)
	$(call cmp_expect,$@.tpm.tsv)

# without map or bam output
test-quantify-nobam : mkdirs
	rm -rf $O/$@.tmp
	flair quantify -r $(BASIC_READS_MANIFEST) -i $(BASIC_ISOFORMS_FA) --isoform_bed $(BASIC_ISOFORMS_BED) --temp_dir $O/$@.tmp --tpm --sample_id_only -o $O/$@
	$(call cmp_expect,$@.tpm.tsv)

test-quantify-help: mkdirs
	flair quantify --help >&$O/$@.out
	$(call cmp_expect,$@.out)

###
# FLAIR DIFFEXP
###
diffexp-tests: test-diffexp test-diffexp-help

test-diffexp : mkdirs
	flair diffexp -q $(BASIC_COUNTS_MATRIX) -o $O/$@ -e 1 -of
	$(call cmp_expect,$@/genes_deseq2_MCF7_v_A549.tsv)
	$(call cmp_expect,$@/isoforms_deseq2_MCF7_v_A549.tsv)
	$(call cmp_expect_sort,$@/isoforms_drimseq_MCF7_v_A549.tsv)

test-diffexp-help: mkdirs
	flair diffexp --help >&$O/$@.out
	$(call cmp_expect,$@.out)

###
# FLAIR DIFFSPLICE
###
diffsplice-tests: test-diffsplice test-diffsplice-help

# NOTE: the order of the ids generated in diffsplice.es.events.quant.tsv is not stable
# (issue #417)
# various output files not long produce due to changes in test cases
test-diffsplice: mkdirs
	flair diffsplice -i $(BASIC_ISOFORMS_BED) -q $(BASIC_DIFFSPLICE_MATRIX) --test -o $O/$@ -of
	$(call cmp_expect_sort,$@/diffsplice.alt3.events.quant.tsv)
	$(call cmp_expect_sort,$@/diffsplice.alt5.events.quant.tsv)
#	$(call cmp_expect,$@/diffsplice.es.events.quant.tsv)
	$(call cmp_expect_sort,$@/diffsplice.ir.events.quant.tsv)
#	$(call cmp_expect,$@/drimseq_alt3_A_v_B.tsv)
#	$(call cmp_expect,$@/drimseq_es_A_v_B.tsv)

test-diffsplice-help: mkdirs
	flair diffsplice --help >&$O/$@.out
	$(call cmp_expect,$@.out)


##
# PREDICT PRODUCTIVITY
##
predict-productivity-tests: test-predict-productivity test-predict-productivity-help

test-predict-productivity : mkdirs ${GENOME_FA}
	predictProductivity -i $(BASIC_ISOFORMS_BED) -g $(BASIC_ANNOTATION) -f ${GENOME_FA} -o $O/$@ --longestORF
	$(call cmp_expect,$@.bed)
	$(call cmp_expect,$@.info.tsv)

test-predict-productivity-help: mkdirs
	predictProductivity --help >&$O/$@.out
	$(call cmp_expect,$@.out)

##
# DIFF ISO USAGE
##
diff-iso-usage-tests: test-diff-iso-usage test-diff-iso-usage-help

test-diff-iso-usage: mkdirs
	diff_iso_usage $(BASIC_COUNTS_MATRIX) MCF7cdnar1r3_MCF7_batch1 MCF7cdnastrandr3r2_MCF7_batch1 $O/$@.tsv
	$(call cmp_expect,$@.tsv)

test-diff-iso-usage-help: mkdirs
	diff_iso_usage --help >&$O/$@.out
	$(call cmp_expect,$@.out)

##
# plot isoform usage
##
plot-usage-tests: test-plot-usage

test-plot-usage: mkdirs
	plot_isoform_usage expected/test-predict-productivity.bed expected/test-collapse-annot.annotated_transcripts.alignment.counts ENSG00000133703.12 -o output/$@
	file output/test-plot-usage_isoforms.png | grep 'PNG image data' >/dev/null
	file output/test-plot-usage_usage.png | grep 'PNG image data' >/dev/null


other-progs-tests: test-bam2bed12 test-transcriptome test-fusion test-variants


###
# flair transcriptome
###


test-transcriptome: mkdirs
	flair transcriptome  -b output/test-align.bam -g ${GENOME_FA} -t 4 -f $(BASIC_ANNOTATION) -o $O/$@ --stringent
	$(call cmp_expect,$@.isoforms.bed)

###
# flair variants
# FIXME need to add better checks for AAseq results, however the inconsistent sorting of the variants leads to inconsistent AAseq IDs which are hard to check
###

test-variants: mkdirs
	minimap2 -ax map-ont --secondary=no ${GENOME_FA} input/basic.reads.4.fa | samtools view -hb - | samtools sort - > $O/$@.genomealigned.bam
	samtools index $O/$@.genomealigned.bam
	samtools faidx $(BASIC_ISOFORMS_FA)
	longshot --force_overwrite --bam $O/$@.genomealigned.bam --ref ${GENOME_FA} --out $O/$@.genomealigned.longshot.vcf --min_cov 3 --min_alt_count 3 --strand_bias_pvalue_cutoff 0.000001
	echo 's1    output/test-quantify.B1.B.flair.aligned.bam   $O/$@.genomealigned.longshot.vcf' > $O/$@.varmanifest.tsv
	flair variants -m $O/$@.varmanifest.tsv -i $(BASIC_ISOFORMS_FA) -b $(BASIC_ISOFORMS_BED) -g ${GENOME_FA} -f $(BASIC_ANNOTATION) -o $O/$@
	$(call cmp_expect,$@.isoswithvars.counts.tsv)
	$(call cmp_expect,$@.aaseq.counts.tsv)


###
# flair fusion
###
test-fusion: mkdirs
	flair align -r $(BASIC_READS_FA) --genome ${GENOME_FA} -t 4 -o $O/$@ --filtertype separate --minfragmentsize 40 --maxintronlen 350k
	flair fusion -g ${GENOME_FA} -f $(BASIC_ANNOTATION) -r $(BASIC_READS_FA) -b $O/$@_chimeric.bam --maxloci 4 -o $O/$@ --annotated_fa output/test-collapse-annot.annotated_transcripts.fa -s 2
	$(call cmp_expect,$@-syntheticBreakpointLoc.bed)

##
# Other programs
##

test-bam2bed12: mkdirs
	samtools view -b input/tiny.sam >output/$@.bam
	samtools index output/$@.bam
	bam2Bed12 -i output/$@.bam > output/$@.bed
	$(call cmp_expect,$@.bed)

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

mkdirs:
	@mkdir -p $O

clean:
	rm -rf $O $D

realclean : clean
	rm -rf ${TI}
