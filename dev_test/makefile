SHELL := /bin/bash

FLAIR=flair

R = test_input
E = test_expected
O = test_output
D = test_diffs
READS_FQ = $R/reads.fq
GENOME = $R/genome.fa 
ANNOTATION = $R/annotation.gtf
SHORTREAD_JUNCTIONS = $R/sr_junctions.bed
UNCORRECTED_READS = $R/uncorrected_reads.bed
READS_BED = $R/corrected_reads.bed # output from correct
PROMOTERS = $R/gc.starts
GENE_ENDS = $R/gc.ends
#ISOFORMS_FA = $R/collapse.isoforms.fa
#ISOFORMS_BED = $R/collapse.isoforms.bed
#READS_MANIFEST = $R/reads_manifest.txt
#COUNTS_MATRIX = $R/counts_matrix_diffexp.tsv
#DIFFSPLICE_MATRIX = $R/counts_matrix_diffsplice.tsv

test_input/genome.fa :
	rm -f test_input/genome.fa.fai
	wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr7.fa.gz
	wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr13.fa.gz
	zcat chr7.fa.gz chr13.fa.gz > test_input/genome.fa
	rm chr7.fa.gz chr13.fa.gz 

test_input/reads.fq :
	gunzip test_input/reads.fq.gz

mkdirs :
	mkdir -p $O
	mkdir -p $D

# FLAIR ALIGN TODO: make tiny input 
test-align : mkdirs test_input/genome.fa test_input/reads.fq
	$(FLAIR) align --reads $(READS_FQ) --genome $(GENOME) --threads 8 --output $O/align
	diff <(sort $E/align.bed) <(sort $O/align.bed) > $D/align.diff

# junction_bed improves output dramatically
test-align-jncts : mkdirs test_input/genome.fa
	$(FLAIR) align --junction_bed $(SHORTREAD_JUNCTIONS) --reads $(READS_FQ) --genome $(GENOME) --threads 8 --output $O/align_jncts
	diff <(sort $E/align_jncts.bed) <(sort $O/align_jncts.bed) > $D/align_jncts.diff

# FLAIR CORRECT
test-correct : mkdirs 
	$(FLAIR) correct --query $(UNCORRECTED_READS) --gtf $(ANNOTATION) --genome $(GENOME) --shortread $(SHORTREAD_JUNCTIONS) --threads 8 --output $O/correct
	diff <(sort $E/correct_all_corrected.bed) <(sort $O/correct_all_corrected.bed) > $D/correct.diff
	diff <(sort $E/correct_all_inconsistent.bed) <(sort $O/correct_all_inconsistent.bed) > $D/correct.diff

# FLAIR COMBI
#flair 123 \
#            -g $ref/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna \
#            -t 1024 \
#            -r $input_dir/$sample-gridion-final_pass.fastq.gz \
#            -f $ref/GCA_000001405.15_GRCh38_full_analysis_set.refseq_annotation.gtf \
#            -m \
#            -nvrna \
#            --generate_map \
#            --annotated_bed \
#            -o flair4/$sample-flair.output \
#            --temp_dir flair4/temp_dir \
#            --keep_intermediate
test-flair123 : mkdirs
	$(FLAIR) 123 --genome $(GENOME) --reads $(READS_FQ) --gtf $(ANNOTATION) -m -nvrna --generate_map --annotated_bed --output $O/flair123 --keep-intermediate

# FLAIR COLLAPSE
# TODO: Flair can be run without a gtf input but the result is quite bad. Maybe require?
test-collapse-basic : mkdirs
	$(FLAIR) collapse --gtf $(ANNOTATION) --reads $(READS_FQ) --query $(READS_BED) --genome $(GENOME) --threads 4 --output $O/collapse 
	diff <(sort $E/collapse.isoforms.bed) <(sort $O/collapse.isoforms.bed) > $D/collapse.diff

#test-collapse-annotation-reliant-infile : mkdirs
#	$(FLAIR) collapse --gtf $(ANNOTATION) --reads $(READS_FQ) --query $(READS_BED) --genome $(GENOME) --threads 4 --output $O/collapse_annotation_reliant_infile --annotation_reliant <file>

# TODO: annotation_reliant generate always outputs isoform.read.map.txt, combined.isoform.read.map.txt, and annotated_transcripts.read.map.txt. Necessary?
test-collapse-annotation-reliant-generate : mkdirs
	$(FLAIR) collapse --gtf $(ANNOTATION) --reads $(READS_FQ) --query $(READS_BED) --genome $(GENOME) --threads 4 --output $O/collapse_annotation_reliant_generate --annotation_reliant generate --keep_intermediate
	diff <(sort $E/collapse_annotation_reliant_generate.isoforms.bed) <(sort $O/collapse_annotation_reliant_generate.isoforms.bed) > $D/collapse_annotation_reliant_generate.diff

# this has the exact same output as without check splice, TODO figure this out
#test-collapse-annotation-reliant-generate-checksplice : mkdirs
#	$(FLAIR) collapse --reads $(READS_FQ) --query $(READS_BED) --genome $(GENOME) --threads 4 --output $O/collapse_annotation_reliant_generate_checksplice --annotation_reliant generate --check_splice --gtf $(ANNOTATION)

test-collapse-support : mkdirs
	$(FLAIR) collapse --gtf $(ANNOTATION) --reads $(READS_FQ) --query $(READS_BED) --genome $(GENOME) --threads 4 --output $O/collapse_support --support 5 
	diff <(sort $E/collapse_support.isoforms.bed) <(sort $O/collapse_support.isoforms.bed) > $D/collapse_support.diff

# fract requires input gtf (TODO add to README)
test-collapse-support-fract : mkdirs
	$(FLAIR) collapse --gtf $(ANNOTATION) --reads $(READS_FQ) --query $(READS_BED) --genome $(GENOME) --threads 4 --output $O/collapse_support_fract --support 0.8
	diff <(sort $E/collapse_support_fract.isoforms.bed) <(sort $O/collapse_support_fract.isoforms.bed) > $D/collapse_support_fract.diff

test-collapse-stringent : mkdirs
	$(FLAIR) collapse --gtf $(ANNOTATION) --reads $(READS_FQ) --query $(READS_BED) --genome $(GENOME) --threads 4 --output $O/collapse_stringent --stringent 
	diff <(sort $E/collapse_stringent.isoforms.bed) <(sort $O/collapse_stringent.isoforms.bed) > $D/collapse_stringent.diff

test-collapse-trust-ends : mkdirs
	$(FLAIR) collapse --gtf $(ANNOTATION) --reads $(READS_FQ) --query $(READS_BED) --genome $(GENOME) --threads 4 --output $O/collapse_trust_ends --trust_ends
	diff <(sort $E/collapse_trust_ends.isoforms.bed) <(sort $O/collapse_trust_ends.isoforms.bed) > $D/collapse_trust_ends.diff

test-collapse-quality : mkdirs
	$(FLAIR) collapse --gtf $(ANNOTATION) --reads $(READS_FQ) --query $(READS_BED) --genome $(GENOME) --threads 4 --output $O/collapse_quality --quality 2
	diff <(sort $E/collapse_quality.isoforms.bed) <(sort $O/collapse_quality.isoforms.bed) > $D/collapse_quality.diff

test-collapse-end-window : mkdirs
	$(FLAIR) collapse --gtf $(ANNOTATION) --reads $(READS_FQ) --query $(READS_BED) --genome $(GENOME) --threads 4 --output $O/collapse_end_window --end_window 50
	diff <(sort $E/collapse_end_window.isoforms.bed) <(sort $O/collapse_end_window.isoforms.bed) > $D/collapse_end_window.diff

test-collapse-promoters : mkdirs
	$(FLAIR) collapse --gtf $(ANNOTATION) --reads $(READS_FQ) --query $(READS_BED) --genome $(GENOME) --threads 4 --output $O/collapse_promoters --promoters $(PROMOTERS) 
	diff <(sort $E/collapse_promoters.isoforms.bed) <(sort $O/collapse_promoters.isoforms.bed) > $D/collapse_promoters.diff

test-collapse-3prime-regions : mkdirs
	$(FLAIR) collapse --gtf $(ANNOTATION) --reads $(READS_FQ) --query $(READS_BED) --genome $(GENOME) --threads 4 --output $O/collapse_3prime_regions --3prime_regions $(GENE_ENDS)
	diff <(sort $E/collapse_3prime_regions.isoforms.bed) <(sort $O/collapse_3prime_regions.isoforms.bed) > $D/collapse_3prime_regions.diff

test-collapse-bothsides : mkdirs
	$(FLAIR) collapse --generate_map --gtf $(ANNOTATION) --reads $(READS_FQ) --query $(READS_BED) --genome $(GENOME) --threads 4 --output $O/collapse_bothends --promoters $(PROMOTERS) --3prime_regions $(GENE_ENDS)
	diff <(sort $E/collapse_bothends.isoforms.bed) <(sort $O/collapse_bothends.isoforms.bed) > $D/collapse_bothends.diff

test-collapse-no-redundant-longest : mkdirs
	$(FLAIR) collapse --gtf $(ANNOTATION) --reads $(READS_FQ) --query $(READS_BED) --genome $(GENOME) --threads 4 --output $O/collapse_no_redundant_longest --no_redundant longest
	diff <(sort $E/collapse_no_redundant_longest.isoforms.bed) <(sort $O/collapse_no_redundant_longest.isoforms.bed) > $D/collapse_no_redundant_longest.diff

# TODO: if you say best-only it doesn't warn, FIX
test-collapse-no-redundant-bestonly : mkdirs
	$(FLAIR) collapse --gtf $(ANNOTATION) --reads $(READS_FQ) --query $(READS_BED) --genome $(GENOME) --threads 4 --output $O/collapse_no_redundant_bestonly --no_redundant best_only
	diff <(sort $E/collapse_no_redundant_bestonly.isoforms.bed) <(sort $O/collapse_no_redundant_bestonly.isoforms.bed) > $D/collapse_no_redundant_bestonly.diff

test-collapse-isoformtss : mkdirs
	$(FLAIR) collapse --gtf $(ANNOTATION) --reads $(READS_FQ) --query $(READS_BED) --genome $(GENOME) --threads 4 --output $O/collapse_isoformtss --isoformtss
	diff <(sort $E/collapse_isoformtss.isoforms.bed) <(sort $O/collapse_isoformtss.isoforms.bed) > $D/collapse_isoformtss.diff

test-collapse-no-gtf-end-adjustment : mkdirs
	$(FLAIR) collapse --gtf $(ANNOTATION) --reads $(READS_FQ) --query $(READS_BED) --genome $(GENOME) --threads 4 --output $O/collapse_no_gtf_end_adjustment --no_gtf_end_adjustment
	diff <(sort $E/collapse_no_gtf_end_adjustment.isoforms.bed) <(sort $O/collapse_no_gtf_end_adjustment.isoforms.bed) > $D/collapse_no_gtf_end_adjustment.diff

test-collapse-max-ends : mkdirs
	$(FLAIR) collapse --gtf $(ANNOTATION) --reads $(READS_FQ) --query $(READS_BED) --genome $(GENOME) --threads 4 --output $O/collapse_max_ends --max_ends 1
	diff <(sort $E/collapse_max_ends.isoforms.bed) <(sort $O/collapse_max_ends.isoforms.bed) > $D/collapse_max_ends.diff

test-collapse-filter-nosubset : mkdirs
	$(FLAIR) collapse --gtf $(ANNOTATION) --reads $(READS_FQ) --query $(READS_BED) --genome $(GENOME) --threads 4 --output $O/collapse_filter_nosubset --filter nosubset
	diff <(sort $E/collapse_filter_nosubset.isoforms.bed) <(sort $O/collapse_filter_nosubset.isoforms.bed) > $D/collapse_filter_nosubset.diff

test-collapse-filter-comprehensive : mkdirs
	$(FLAIR) collapse --gtf $(ANNOTATION) --reads $(READS_FQ) --query $(READS_BED) --genome $(GENOME) --threads 4 --output $O/collapse_filter_comprehensive --filter comprehensive
	diff <(sort $E/collapse_filter_comprehensive.isoforms.bed) <(sort $O/collapse_filter_comprehensive.isoforms.bed) > $D/collapse_filter_comprehensive.diff

test-collapse-filter-ginormous : mkdirs
	$(FLAIR) collapse --gtf $(ANNOTATION) --reads $(READS_FQ) --query $(READS_BED) --genome $(GENOME) --threads 4 --output $O/collapse_filter_ginormous --filter ginormous
	diff <(sort $E/collapse_filter_ginormous.isoforms.bed) <(sort $O/collapse_filter_ginormous.isoforms.bed) > $D/collapse_filter_ginormous.diff


#test-align-jncts  test-correct-jncts test-collapse-annotation-reliant-generate-checksplice
test-collapse : test-collapse-basic test-collapse-annotation-reliant-generate  test-collapse-support test-collapse-support-fract  test-collapse-stringent test-collapse-trust-ends test-collapse-quality test-collapse-end-window test-collapse-promoters test-collapse-3prime-regions test-collapse-bothsides test-collapse-no-redundant-longest test-collapse-no-redundant-bestonly test-collapse-isoformtss test-collapse-no-gtf-end-adjustment test-collapse-max-ends test-collapse-filter-nosubset test-collapse-filter-comprehensive test-collapse-filter-ginormous

clean :
	rm -rf $O/* $D/*
