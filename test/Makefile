root = ..
include ${root}/defs.mk

# If you want to use the install flair rather than the one it this
# tree use:
#   make test use_installed_flair=yes
ifneq (${use_installed_flair},yes)
    export PATH := ${root}/src/flair:${root}/bin:${PATH}
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
SEG1_READS_MANIFEST = input/seg1.reads_manifest.txt

SEG1_ALIGN_BED = $E/test-align-seg1.bed
SEG1_CORRECT_READS_BED = $E/test-correct-seg1_all_corrected.bed
SEG1_ISOFORMS_FA = expected/test-collapse-seg1.isoforms.fa
SEG1_ISOFORMS_BED = expected/test-collapse-seg1.isoforms.bed


# The help test results vary on screen width.  This forces them to be consistent.
# It might seem weird to see if help works, but some it did break in the past.
export COLUMNS = 80


default:
	@echo To run all tests, use: >&2
	@echo "   " make test >&2
	@echo To run tests that do not depend on R >&2
	@echo "   " make base-tests >&2
	@echo otherwise run individual tests by specifying test target, such as: >&2
	@echo "   " make align-tests >&2
	@exit 1

base-tests: main-tests align-tests correct-tests collapse-tests quantify-tests \
	predict-productivity-tests diff-iso-usage-tests
test : base-tests diffexp-tests diffsplice-tests

###
# MAIN COMMAND
###
main-tests: test-main-help

test-main-help: mkdirs
	flair --help >& $O/$@.out
	diff $E/$@.out $O/$@.out

###
# FLAIR ALIGN
###

align-tests: test-align test-align-intpriming test-align-seg1 test-align-help

test-align : mkdirs ${GENOME_FA}
	flair align -r $(BASIC_READS_FA) --genome ${GENOME_FA} -t 4 -o $O/$@
	diff $E/$@.bed $O/$@.bed

test-align-intpriming : mkdirs ${GENOME_FA}
	flair align -r $(BASIC_READS_FA) --genome ${GENOME_FA} -t 4 -o $O/$@ --remove_internal_priming -f $(BASIC_ANNOTATION) --remove_singleexon
	diff $E/$@.bed $O/$@.bed

test-align-seg1 : mkdirs ${GENOME_FA}
	flair align -r ${SEG1_READS_FQ} --genome ${GENOME_FA} -t 4 -o $O/$@
	diff $E/$@.bed $O/$@.bed

test-align-help: mkdirs
	flair align --help >&$O/$@.out
	diff $E/$@.out $O/$@.out

###
# FLAIR CORRECT
###
correct-tests: test-correct test-correct-help

test-correct: mkdirs ${GENOME_FA}
	flair correct -q ${BASIC_ALIGN_BED} -j $(BASIC_SHORTREAD_JUNCTIONS) -f $(BASIC_ANNOTATION_INCOMPLETE) -g ${GENOME_FA} -o $O/$@
	diff <(sort $E/$@_all_inconsistent.bed) <(sort $O/$@_all_inconsistent.bed)
	diff <(sort $E/$@_all_corrected.bed) <(sort $O/$@_all_corrected.bed)

test-correct-seg1: mkdirs ${GENOME_FA}
	flair correct -q ${SEG1_ALIGN_BED} -g ${GENOME_FA} -f ${SEG1_ANNOTATION} -o $O/$@
	diff <(sort $E/$@_all_inconsistent.bed) <(sort $O/$@_all_inconsistent.bed)
	diff <(sort $E/$@_all_corrected.bed) <(sort $O/$@_all_corrected.bed)

test-correct-help: mkdirs
	flair correct --help >&$O/$@.out
	diff $E/$@.out $O/$@.out

###
# FLAIR COLLAPSE
###
collapse-tests: test-collapse test-collapse-annot test-collapse-faster test-collapse-faster-annot test-collapse-faster-annot-removeintprim \
	test-collapse-support test-collapse-seg1 \
	test-collapse-help test-collapse-faster-help

# basic collapse
test-collapse: mkdirs ${GENOME_FA}
	rm -rf $O/$@.tmp
	flair collapse -r $(BASIC_READS_FA) -q $(BASIC_CORRECT_READS_BED) -g ${GENOME_FA} -t 4 --generate_map --temp_dir $O/$@.tmp --keep_intermediate -f $(BASIC_ANNOTATION) -o $O/$@ -p $(BASIC_PROMOTER_REGIONS)
	diff $E/$@.isoforms.gtf $O/$@.isoforms.gtf

# collapse annotation-reliant
##these reads are low quality, many isoforms fail with the check_splice parameter
test-collapse-annot: mkdirs ${GENOME_FA}
	rm -rf $O/$@.tmp
	flair collapse -r $(BASIC_READS_FA) -q $(BASIC_CORRECT_READS_BED) -g ${GENOME_FA} -t 4 --generate_map --temp_dir $O/$@.tmp --keep_intermediate -f $(BASIC_ANNOTATION) -o $O/$@ -p $(BASIC_PROMOTER_REGIONS) --annotation_reliant generate --stringent
	diff $E/$@.isoforms.bed $O/$@.isoforms.bed
	diff $E/$@.isoforms.gtf $O/$@.isoforms.gtf
	diff $E/$@.isoforms.fa $O/$@.isoforms.fa

test-collapse-seg1: mkdirs ${GENOME_FA}
	rm -rf $O/$@.tmp
	@echo "WARNING: test-collapse-seg1 not currently working" >&2
	#flair collapse -r ${SEG1_READS_FQ} -q ${SEG1_CORRECT_READS_BED} -g ${GENOME_FA} -t 4 --generate_map --temp_dir $O/$@.tmp --keep_intermediate -f $(SEG1_ANNOTATION) -o $O/$@ -p ${SEG1_PROMOTER_REGIONS} --annotation_reliant generate --stringent

test-collapse-help: mkdirs
	flair collapse --help >&$O/$@.out
	diff $E/$@.out $O/$@.out


# speedy collapse
###this works badly on the test dataset because the test reads are error-prone and fail the built-in correct step
test-collapse-faster : mkdirs ${GENOME_FA}
	rm -rf $O/$@.tmp
	${PYTHON} ../src/flair/fastercollapse.py -r $(BASIC_READS_FA) -g ${GENOME_FA} -t 4 --generate_map --temp_dir $O/$@.tmp --keep_intermediate -f $(BASIC_ANNOTATION) -o $O/$@ -p $(BASIC_PROMOTER_REGIONS)
	diff $E/$@.isoforms.gtf $O/$@.isoforms.gtf

test-collapse-faster-help: mkdirs
	${PYTHON} ../src/flair/fastercollapse.py --help >&$O/$@.out
	diff $E/$@.out $O/$@.out

# speedy collapse annotation-reliant
test-collapse-faster-annot : mkdirs ${GENOME_FA}
	rm -rf $O/$@.tmo
	${PYTHON} ../src/flair/fastercollapse.py -r $(BASIC_READS_FA) --transcriptfasta generate --stringent -g ${GENOME_FA} -t 4 --generate_map --temp_dir $O/$@.tmp --keep_intermediate -f $(BASIC_ANNOTATION) -o $O/$@ -p $(BASIC_PROMOTER_REGIONS)
	diff $E/$@.isoforms.gtf $O/$@.isoforms.gtf

test-collapse-faster-annot-removeintprim : mkdirs ${GENOME_FA}
	rm -rf $O/$@.tmp
	${PYTHON} ../src/flair/fastercollapse.py -r $(BASIC_READS_FA) --transcriptfasta generate --stringent -g ${GENOME_FA} -t 4 --generate_map --temp_dir $O/$@.tmp --keep_intermediate -f $(BASIC_ANNOTATION) -o $O/$@ -p $(BASIC_PROMOTER_REGIONS) --remove_internal_priming
	diff $E/$@.isoforms.gtf $O/$@.isoforms.gtf

# collapse with support
# Note: output order not stable
test-collapse-support : mkdirs ${GENOME_FA}
	rm -rf $O/$@.tmp
	flair collapse -s 0.25 -r $(BASIC_READS_FA) -q $(BASIC_CORRECT_READS_BED) -g ${GENOME_FA} -t 4 --generate_map --temp_dir $O/$@.tmp --keep_intermediate -f $(BASIC_ANNOTATION) -o $O/$@ -p $(BASIC_PROMOTER_REGIONS)
	diff <(sort $E/$@.isoforms.gtf) <(sort $O/$@.isoforms.gtf)

###
# FLAIR QUANTIFY
###
quantify-tests: test-quantify test-quantify-nobam test-quantify-seg1 test-quantify-help

test-quantify: mkdirs
	rm -rf $O/$@.tmp
	flair quantify -r $(BASIC_READS_MANIFEST) -i $(BASIC_ISOFORMS_FA) --isoform_bed $(BASIC_ISOFORMS_BED) --generate_map --temp_dir $O/$@.tmp --tpm --sample_id_only -o $O/$@ --output_bam
	samtools view $O/$@.A1.A.flair.aligned.bam | cut -f 1,3 > $O/$@.A1.A.flair.aligned.readtotranscript.txt
	diff $E/$@.tpm.tsv $O/$@.tpm.tsv
	diff $E/$@.A1.A.flair.aligned.readtotranscript.txt $O/$@.A1.A.flair.aligned.readtotranscript.txt

# without map or bam output
test-quantify-nobam : mkdirs
	rm -rf $O/$@.tmp
	flair quantify -r $(BASIC_READS_MANIFEST) -i $(BASIC_ISOFORMS_FA) --isoform_bed $(BASIC_ISOFORMS_BED) --temp_dir $O/$@.tmp --tpm --sample_id_only -o $O/$@
	diff $E/$@.tpm.tsv $O/$@.tpm.tsv

# currently produces no output, (#437)
test-quantify-seg1: mkdirs
	rm -rf $O/$@.tmp
	@echo "WARNING: test-quantify-seg1 not currently working" >&2
	#flair quantify -r ${SEG1_READS_MANIFEST} -i ${SEG1_ISOFORMS_FA} --isoform_bed ${SEG1_ISOFORMS_BED} --generate_map --temp_dir $O/$@.tmp --tpm --sample_id_only -o $O/$@ --output_bam

test-quantify-help: mkdirs
	flair quantify --help >&$O/$@.out
	diff $E/$@.out $O/$@.out

###
# FLAIR DIFFEXP
###
diffexp-tests: test-diffexp test-diffexp-help

test-diffexp : mkdirs
	flair_diffExp -q $(BASIC_COUNTS_MATRIX) -o $O/$@ -e 1 -of
	diff $E/$@/genes_deseq2_MCF7_v_A549.tsv $O/$@/genes_deseq2_MCF7_v_A549.tsv
	diff $E/$@/isoforms_deseq2_MCF7_v_A549.tsv $O/$@/isoforms_deseq2_MCF7_v_A549.tsv
	diff $E/$@/isoforms_drimseq_MCF7_v_A549.tsv $O/$@/isoforms_drimseq_MCF7_v_A549.tsv

test-diffexp-help: mkdirs
	flair_diffExp --help >&$O/$@.out
	diff $E/$@.out $O/$@.out

###
# FLAIR DIFFSPLICE
###
diffsplice-tests: test-diffsplice test-diffsplice-help

# NOTE: the order of the ids generated in diffsplice.es.events.quant.tsv is not stable
# (issue #417)
test-diffsplice : mkdirs
	flair_diffSplice -i $(BASIC_ISOFORMS_BED) -q $(BASIC_DIFFSPLICE_MATRIX) --test -o $O/$@ -of
	diff $E/$@/diffsplice.alt3.events.quant.tsv $O/$@/diffsplice.alt3.events.quant.tsv
	diff $E/$@/diffsplice.alt5.events.quant.tsv $O/$@/diffsplice.alt5.events.quant.tsv
	@echo "WARNING: diffsplice.es.events.quant.tsv check disabled due to instability in id generation" >&2
	#diff $E/$@/diffsplice.es.events.quant.tsv $O/$@/diffsplice.es.events.quant.tsv
	diff $E/$@/diffsplice.ir.events.quant.tsv $O/$@/diffsplice.ir.events.quant.tsv
	diff $E/$@/drimseq_alt3_A_v_B.tsv $O/$@/drimseq_alt3_A_v_B.tsv
	diff $E/$@/drimseq_es_A_v_B.tsv $O/$@/drimseq_es_A_v_B.tsv

test-diffsplice-help: mkdirs
	flair_diffSplice --help >&$O/$@.out
	diff $E/$@.out $O/$@.out


##
# PREDICT PRODUCTIVITY
##
predict-productivity-tests: test-predict-productivity test-predict-productivity-help

test-predict-productivity : mkdirs ${GENOME_FA}
	predictProductivity -i $(BASIC_ISOFORMS_BED) -g $(BASIC_ANNOTATION) -f ${GENOME_FA} -o $O/$@ --longestORF
	diff $E/$@.bed $O/$@.bed
	diff $E/$@.info.tsv $O/$@.info.tsv

test-predict-productivity-help: mkdirs
	predictProductivity --help >&$O/$@.out
	diff $E/$@.out $O/$@.out

##
# DIFF ISO USAGE
##
diff-iso-usage-tests: test-diff-iso-usage test-diff-iso-usage-help

test-diff-iso-usage: mkdirs
	diff_iso_usage $(BASIC_COUNTS_MATRIX) MCF7cdnar1r3_MCF7_batch1 MCF7cdnastrandr3r2_MCF7_batch1 $O/$@.tsv
	diff $E/$@.tsv $O/$@.tsv

test-diff-iso-usage-help: mkdirs
	diff_iso_usage --help >&$O/$@.out
	diff $E/$@.out $O/$@.out

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
