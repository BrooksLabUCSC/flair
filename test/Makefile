root = ..
include ${root}/defs.mk

# If you want to use the install flair rather than the one it this
# tree use:
#   make test use_installed_flair=yes
ifneq (${use_installed_flair},yes)
    export PATH := ${root}/src/flair:${root}/bin:${PATH}
endif

I = input
E = expected
O = output
D = diffs
IT = input.tmp
READS_FA = $I/reads.fa
ANNOTATION = $I/annotation.gtf
ANNOTATION_INCOMPLETE = $I/annotation.incomplete.gtf
SHORTREAD_JUNCTIONS = $I/shortread_junctions.tab
UNCORRECTED_READS_BED = $I/test.align.bed
READS_BED = $I/test.corrected.bed
# ISOFORMS_FA = $I/test.collapse.isoforms.fa
# ISOFORMS_BED = $I/test.collapse.isoforms.bed
ISOFORMS_FA = $I/test-collapse-annot.isoforms.fa
ISOFORMS_BED = $I/test-collapse-annot.isoforms.bed
PROMOTER_REGIONS = $I/promoter_regions.bed
READS_MANIFEST = $I/reads_manifest.txt
COUNTS_MATRIX = $I/counts_matrix_diffexp.tsv
DIFFSPLICE_MATRIX = $I/counts_matrix_diffsplice.tsv
GENOME_FA = ${IT}/genome.fa

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

align-tests: test-align test-align-intpriming test-align-help

test-align : mkdirs ${GENOME_FA}
	flair align -r $(READS_FA) --genome ${GENOME_FA} -t 4 -o $O/$@
	diff <(sort $E/$@.bed) <(sort $O/$@.bed) > $D/$@.diff

test-align-intpriming : mkdirs ${GENOME_FA}
	flair align -r $(READS_FA) --genome ${GENOME_FA} -t 4 -o $O/$@ --remove_internal_priming -f $(ANNOTATION) --remove_singleexon
	diff <(sort $E/$@.bed) <(sort $O/$@.bed) > $D/$@.diff

test-align-help: mkdirs
	flair align --help >&$O/$@.out
	diff $E/$@.out $O/$@.out

###
# FLAIR CORRECT
###
correct-tests: test-correct test-correct-help

test-correct : mkdirs ${GENOME_FA}
	flair correct -q $(UNCORRECTED_READS_BED) -j $(SHORTREAD_JUNCTIONS) -f $(ANNOTATION_INCOMPLETE) -g ${GENOME_FA} -o $O/$@
	diff <(sort $E/$@_all_inconsistent.bed) <(sort $O/$@_all_inconsistent.bed) > $D/$@.inconsistent.diff
	diff <(sort $E/$@_all_corrected.bed) <(sort $O/$@_all_corrected.bed) > $D/$@.corrected.diff

test-correct-help: mkdirs
	flair correct --help >&$O/$@.out
	diff $E/$@.out $O/$@.out

###
# FLAIR COLLAPSE
###
collapse-tests: test-collapse test-collapse-annot test-collapse-faster test-collapse-faster-annot test-collapse-faster-annot-removeintprim \
	test-collapse-support test-collapse-help test-collapse-faster-help

# basic collapse
test-collapse : mkdirs ${GENOME_FA}
	rm -rf $O/$@.tmp
	flair collapse -r $(READS_FA) -q $(READS_BED) -g ${GENOME_FA} -t 4 --generate_map --temp_dir $O/$@.tmp --keep_intermediate -f $(ANNOTATION) -o $O/$@ -p $(PROMOTER_REGIONS)
	diff <(sort $E/$@.isoforms.gtf) <(sort $O/$@.isoforms.gtf) > $D/$@.isoforms.gtf.diff

# collapse annotation-reliant
##these reads are low quality, many isoforms fail with the check_splice parameter
test-collapse-annot : mkdirs ${GENOME_FA}
	rm -rf output/$@.tmo
	flair collapse -r $(READS_FA) -q $(READS_BED) -g ${GENOME_FA} -t 4 --generate_map --temp_dir $O/$@.tmp --keep_intermediate -f $(ANNOTATION) -o $O/$@ -p $(PROMOTER_REGIONS) --annotation_reliant generate --stringent

test-collapse-help: mkdirs
	flair collapse --help >&$O/$@.out
	diff $E/$@.out $O/$@.out
	diff <(sort $E/$@.isoforms.gtf) <(sort $O/$@.isoforms.gtf) > $D/$@.isoforms.gtf.diff


# speedy collapse
###this works badly on the test dataset because the test reads are error-prone and fail the built-in correct step
test-collapse-faster : mkdirs ${GENOME_FA}
	rm -rf $O/$@.tmp
	${PYTHON} ../src/flair/fastercollapse.py -r $(READS_FA) -g ${GENOME_FA} -t 4 --generate_map --temp_dir $O/$@.tmp --keep_intermediate -f $(ANNOTATION) -o $O/$@ -p $(PROMOTER_REGIONS)
	diff <(sort $E/$@.isoforms.gtf) <(sort $O/$@.isoforms.gtf) > $D/$@.isoforms.gtf.diff

test-collapse-faster-help: mkdirs
	${PYTHON} ../src/flair/fastercollapse.py --help >&$O/$@.out
	diff $E/$@.out $O/$@.out

# speedy collapse annotation-reliant
test-collapse-faster-annot : mkdirs ${GENOME_FA}
	rm -rf $O/$@.tmo
	${PYTHON} ../src/flair/fastercollapse.py -r $(READS_FA) --transcriptfasta generate --stringent -g ${GENOME_FA} -t 4 --generate_map --temp_dir $O/$@.tmp --keep_intermediate -f $(ANNOTATION) -o $O/$@ -p $(PROMOTER_REGIONS)
	diff <(sort $E/$@.isoforms.gtf) <(sort $O/$@.isoforms.gtf) > $D/$@.isoforms.gtf.diff

test-collapse-faster-annot-removeintprim : mkdirs ${GENOME_FA}
	rm -rf $O/$@.tmp
	${PYTHON} ../src/flair/fastercollapse.py -r $(READS_FA) --transcriptfasta generate --stringent -g ${GENOME_FA} -t 4 --generate_map --temp_dir $O/$@.tmp --keep_intermediate -f $(ANNOTATION) -o $O/$@ -p $(PROMOTER_REGIONS) --remove_internal_priming
	diff <(sort $E/$@.isoforms.gtf) <(sort $O/$@.isoforms.gtf) > $D/$@.isoforms.gtf.diff

# collapse with support
test-collapse-support : mkdirs ${GENOME_FA}
	rm -rf $O/$@.tmp
	flair collapse -s 0.25 -r $(READS_FA) -q $(READS_BED) -g ${GENOME_FA} -t 4 --generate_map --temp_dir $O/$@.tmp --keep_intermediate -f $(ANNOTATION) -o $O/$@ -p $(PROMOTER_REGIONS)
	diff <(sort $E/$@.isoforms.gtf) <(sort $O/$@.isoforms.gtf) > $D/$@.gtf.diff

###
# FLAIR QUANTIFY
###
quantify-tests: test-quantify test-quantify-help

test-quantify : mkdirs
	rm -rf $O/$@.tmp
	flair quantify -r $(READS_MANIFEST) -i $(ISOFORMS_FA) --isoform_bed $(ISOFORMS_BED) --generate_map --temp_dir $O/$@.tmp --tpm --sample_id_only -o $O/$@ --output_bam
	samtools view $O/$@.A1.A.flair.aligned.bam | cut -f 1,3 > $O/$@.A1.A.flair.aligned.readtotranscript.txt
	diff <(sort $E/$@.tpm.tsv) <(sort $O/$@.tpm.tsv) > $D/$@.diff
	diff <(sort $E/$@.A1.A.flair.aligned.readtotranscript.txt) <(sort $O/$@.A1.A.flair.aligned.readtotranscript.txt) > $D/$@.outputbam.diff

test-quantify-help: mkdirs
	flair quantify --help >&$O/$@.out
	diff $E/$@.out $O/$@.out

###
# FLAIR DIFFEXP
###
diffexp-tests: test-diffexp test-diffexp-help

test-diffexp : mkdirs
	flair_diffExp -q $(COUNTS_MATRIX) -o $O/$@ -e 1 -of
	diff <(sort $E/$@.genes_deseq2_MCF7_v_A549.tsv) <(sort $O/$@/genes_deseq2_MCF7_v_A549.tsv) > $D/$@.dge.diff
	diff <(sort $E/$@.isoforms_deseq2_MCF7_v_A549.tsv) <(sort $O/$@/isoforms_deseq2_MCF7_v_A549.tsv) > $D/$@.die.diff
	diff <(sort $E/$@.isoforms_drimseq_MCF7_v_A549.tsv) <(sort $O/$@/isoforms_drimseq_MCF7_v_A549.tsv) > $D/$@.drimseq.diff

test-diffexp-help: mkdirs
	flair_diffExp --help >&$O/$@.out
	diff $E/$@.out $O/$@.out

###
# FLAIR DIFFSPLICE
###
diffsplice-tests: test-diffsplice test-diffsplice-help

test-diffsplice : mkdirs
	flair_diffSplice -i $(ISOFORMS_BED) -q $(DIFFSPLICE_MATRIX) --test -o $O/$@ -of
	diff <(sort $O/$@/drimseq_alt3_A_v_B.tsv) <(sort $E/$@.drimseq_alt3_A_v_B.tsv) > $D/$@.alt3.diff
	diff <(sort $O/$@/diffsplice.alt3.events.quant.tsv) <(sort $E/$@.alt3.events.quant.tsv) > $D/$@.alt3.quant.diff

test-diffsplice-help: mkdirs
	flair_diffSplice --help >&$O/$@.out
	diff $E/$@.out $O/$@.out


##
# PREDICT PRODUCTIVITY
##
predict-productivity-tests: test-predict-productivity test-predict-productivity-help

test-predict-productivity : mkdirs ${GENOME_FA}
	predictProductivity -i $(ISOFORMS_BED) -g $(ANNOTATION) -f ${GENOME_FA} -o $O/$@ --longestORF
	diff <(sort $E/$@.bed) <(sort $O/$@.bed) > $D/$@.diff
	diff <(sort $E/$@.info.tsv) <(sort $O/$@.info.tsv) > $D/$@.info.diff

test-predict-productivity-help: mkdirs
	predictProductivity --help >&$O/$@.out
	diff $E/$@.out $O/$@.out

##
# DIFF ISO USAGE
##
diff-iso-usage-tests: test-diff-iso-usage test-diff-iso-usage-help

test-diff-iso-usage: mkdirs
	diff_iso_usage $(COUNTS_MATRIX) MCF7cdnar1r3_MCF7_batch1 MCF7cdnastrandr3r2_MCF7_batch1 $O/$@.tsv
	diff <(sort $E/$@.tsv) <(sort $O/$@.tsv) > $D/$@.diff

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

mkdirs :
	@mkdir -p $O $D

clean :
	rm -rf $O $D

realclean : clean
	rm -rf ${IT}
