SHELL := /bin/bash

# These tests assume that flair is in your $PATH
# Flair should be in your path if you conda installed it

I = test_input
E = test_expected
O = test_output
D = test_diffs
READS_FA = $I/reads.fa
ANNOTATION = $I/annotation.gtf
ANNOTATION_INCOMPLETE = $I/annotation.incomplete.gtf
SHORTREAD_JUNCTIONS = $I/shortread_junctions.tab
UNCORRECTED_READS_BED = $I/test.align.bed
READS_BED = $I/test.corrected.bed
ISOFORMS_FA = $I/test.collapse.isoforms.fa
ISOFORMS_BED = $I/test.collapse.isoforms.bed
PROMOTER_REGIONS = $I/promoter_regions.bed
READS_MANIFEST = $I/reads_manifest.txt
COUNTS_MATRIX = $I/counts_matrix_diffexp.tsv
DIFFSPLICE_MATRIX = $I/counts_matrix_diffsplice.tsv

default:
	@echo To run all tests, use: >&2
	@echo "   " make test >&2
	@echo otherwise run individual tests with: >&2
	@echo "   " make test-align >&2
	@echo "   " make test-correct >&2
	@echo "   " make test-collapse >&2
	@echo "   " make test-quantify >&2
	@echo "   " make test-diffexp >&2
	@echo "   " make test-diffsplice >&2
	@echo "   " make test-predict-productivity >&2
	@echo "   " make test-diff-iso-usage >&2
	@exit 1


# FLAIR ALIGN
test-align : mkdirs test_input/genome.fa
	flair align -r $(READS_FA) --genome test_input/genome.fa -t 4 -o $O/test.align
	diff <(sort $E/test.align.bed) <(sort $O/test.align.bed) > $D/test.align.diff

test_input/genome.fa :
	rm -f test_input/genome.fa.fai
	wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr12.fa.gz
	wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr17.fa.gz
	wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz
	zcat chr12.fa.gz chr17.fa.gz chr20.fa.gz > test_input/genome.fa
	rm chr12.fa.gz chr17.fa.gz chr20.fa.gz

mkdirs :
	mkdir -p $O
	mkdir -p $D

# FLAIR CORRECT
test-correct : mkdirs test_input/genome.fa
	flair correct -q $(UNCORRECTED_READS_BED) -j $(SHORTREAD_JUNCTIONS) -f $(ANNOTATION_INCOMPLETE) -g test_input/genome.fa -o $O/test.correct
	diff <(sort $E/test.correct_all_inconsistent.bed) <(sort $O/test.correct_all_inconsistent.bed) > $D/test.correct.inconsistent.diff
	diff <(sort $E/test.correct_all_corrected.bed) <(sort $O/test.correct_all_corrected.bed) > $D/test.correct.corrected.diff

# FLAIR COLLAPSE
test-collapse : mkdirs test_input/genome.fa
	flair collapse -r $(READS_FA) -q $(READS_BED) -g test_input/genome.fa -t 4 --generate_map --temp_dir $O/temp_collapse --keep_intermediate -f $(ANNOTATION) -o $O/test.collapse -p $(PROMOTER_REGIONS)
	diff <(sort $E/test.collapse.isoforms.gtf) <(sort $O/test.collapse.isoforms.gtf) > $D/test.collapse.isoforms.gtf.diff

# FLAIR QUANTIFY
test-quantify : mkdirs
	flair quantify -r $(READS_MANIFEST) -i $(ISOFORMS_FA) --isoform_bed $(ISOFORMS_BED) --generate_map --temp_dir $O/test.quant.temp --tpm --sample_id_only --check_splice --stringent -o $O/test.quantify --output_bam
	samtools view $O/test.quantify.A1.A.flair.aligned.bam | cut -f 1,3 > $O/test.quantify.A1.A.flair.aligned.readtotranscript.txt
	diff <(sort $E/test.quantify.tpm.tsv) <(sort $O/test.quantify.tpm.tsv) > $D/test.quantify.diff
	diff <(sort $E/test.quantify.A1.A.flair.aligned.readtotranscript.txt) <(sort $O/test.quantify.A1.A.flair.aligned.readtotranscript.txt) > $D/test.quantify.outputbam.diff

# FLAIR DIFFEXP
test-diffexp : mkdirs
	flair_diffExp -q $(COUNTS_MATRIX) -o $O/test.diffexp -e 1 -of
	diff <(sort $E/test.diffexp.genes_deseq2_MCF7_v_A549.tsv) <(sort $O/test.diffexp/genes_deseq2_MCF7_v_A549.tsv) > $D/test.diffexp.dge.diff
	diff <(sort $E/test.diffexp.isoforms_deseq2_MCF7_v_A549.tsv) <(sort $O/test.diffexp/isoforms_deseq2_MCF7_v_A549.tsv) > $D/test.diffexp.die.diff
	diff <(sort $E/test.diffexp.isoforms_drimseq_MCF7_v_A549.tsv) <(sort $O/test.diffexp/isoforms_drimseq_MCF7_v_A549.tsv) > $D/test.diffexp.drimseq.diff

# FLAIR DIFFSPLICE
test-diffsplice : mkdirs
	flair_diffSplice -i $(ISOFORMS_BED) -q $(DIFFSPLICE_MATRIX) --test -o $O/test.diffsplice -of
	diff <(sort $O/test.diffsplice/drimseq_alt3_A_v_B.tsv) <(sort $E/test.diffsplice.drimseq_alt3_A_v_B.tsv) > $D/test.diffsplice.alt3.diff
	diff <(sort $O/test.diffsplice/diffsplice.alt3.events.quant.tsv) <(sort $E/test.diffsplice.alt3.events.quant.tsv) > $D/test.diffsplice.alt3.quant.diff

# FLAIR SCRIPTS
test-predict-productivity : mkdirs test_input/genome.fa
	$(DCMD) predictProductivity -i $(ISOFORMS_BED) -g $(ANNOTATION) -f test_input/genome.fa > $O/test.predict_productivity.bed --longestORF
	diff <(sort $E/test.predict_productivity.bed) <(sort $O/test.predict_productivity.bed) > $D/test.predict_productivity.diff

test-diff-iso-usage : mkdirs
	$(DCMD) diff_iso_usage $(COUNTS_MATRIX) MCF7cdnar1r3_MCF7_batch1 MCF7cdnastrandr3r2_MCF7_batch1 $O/test.diff_iso_usage.tsv
	diff <(sort $E/test.diff_iso_usage.tsv) <(sort $O/test.diff_iso_usage.tsv) > $D/test.diff.iso.usage.diff

test : test-align test-correct test-collapse test-quantify test-diffexp test-diffsplice test-predict-productivity test-diff-iso-usage

clean :
	rm -rf $O/* $D/*
