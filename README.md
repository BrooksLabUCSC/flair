# flair
FLAIR (Full-Length Alternative Isoform analysis of RNA) for the correction, isoform definition, and alternative splicing analysis of noisy reads. FLAIR has primarily been used for nanopore cDNA, native RNA, and PacBio sequencing reads. 

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [FLAIR modules](#modules)
	- [align](#align)
	- [correct](#correct)
		- [short-read junctions](#short)
	- [collapse](#collapse)
	- [quantify](#quantify)
	- [diffExp](#diffExp)
	- [diffSplice](#diffSplice)
- [Scripts](#scripts)
- [Docker](#docker)
- [Conda environment](#condaenv)
- [Example Files](#exfiles)
- [Cite FLAIR](#cite)

## <a name="overview"></a>Overview
FLAIR can be run optionally with short-read data to help increase splice site accuracy of the long read splice junctions. FLAIR uses multiple alignment steps and splice site filters to increase confidence in the set of isoforms defined from noisy data. FLAIR was designed to be able to sense subtle splicing changes in nanopore data from [Tang et al. (2018)](https://www.biorxiv.org/content/early/2018/09/06/410183). Please read for more description of the methods.

<img src='https://people.ucsc.edu/~atang14/flair/flair_workflow_compartmentalized.png' alt='flair workflow' width='680'/>

It is recommended to combine all samples together prior to running flair-collapse for isoform assembly by concatenating corrected read `psl` files together. Following the creation of an isoform reference from flair-collapse, consequent steps will assign reads from each sample individually to isoforms of the combined assembly for downstream analyses.

It is also good to note that `bed12` and `PSL` can be converted using [kentUtils](https://github.com/ENCODE-DCC/kentUtils/tree/master/src/hg/utils) bedToPsl or pslToBed, or using `bin/bed_to_psl.py` and `bin/psl_to_bed.py`.

## <a name="requirements"></a>Requirements

1. python v2.7+ and python modules: intervaltree, kerneltree, tqdm, pybedtools, pysam v0.8.4+
2. bedtools, samtools
3. [minimap2](https://github.com/lh3/minimap2)

## <a name="modules"></a>FLAIR modules 
`flair.py` is a wrapper script with modules for running various processing scripts located in `bin/`. Modules are assumed to be run in order (align, correct, collapse), but the user can forgo the wrapper if a more custom build is desired. 

### <a name="align"></a>flair align
Aligns reads to the genome using minimap2, and converts the aligned minimap2 `sam` output to [BED12](https://genome.ucsc.edu/FAQ/FAQformat.html#format14) and optionally [PSL](https://genome.ucsc.edu/FAQ/FAQformat.html#format2). Aligned reads in `psl` format can be visualized in IGV or the UCSC Genome browser. As for which human reference genome to use, Heng Li has written a [blog post](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use) on this.

Alternatively, the user can align the reads themselves with their aligner of choice and convert sorted `bam` output to `bed12` using `bin/bam2Bed12.py` to supply for flair-correct. This step smooths gaps in the alignment.

**Usage:**
```sh
python flair.py align -g genome.fa -r <reads.fq>|<reads.fa> [options]
```
run with `--help` for a description of optional arguments. Outputs (1) `sam` of raw aligned reads and (2) smoothed `bed12` file of aligned reads to be supplied to flair-correct.

### <a name="correct"></a>flair correct
Corrects misaligned splice sites using genome annotations and/or short-read splice junctions. Based on common user issues we have encountered, for flair-correct to run properly, please ensure/note that (1) the genome annotation and genome sequences are compatible, (2) `gtf` is preferred over `gff` for annotation and annotations that do not split single exons into multiple entries are ideal, (3) Bedtools is in your $PATH, and (4) kerneltree is properly installed (you may need to install Cython first). You may also want to refer to the [installation requirements](#requirements) and/or use the [conda environment](#condaenv) for flair.

**Usage:**
```sh
python flair.py correct -q query.bed12 -g genome.fa [options]
```
run with `--help` for description of optional arguments.
Outputs (1) `bed12` of corrected reads, (2) `bed12` of reads that weren't able to be corrected, (3) `psl` of corrected reads if the -c chromsizes file is provided. Either (1) or (3) can be supplied to flair-collapse as the query.

#### <a name="short"></a>Short-read junctions
To use short-read splice sites to aid with correction, one option is `bin/junctions_from_sam.py` to extract splice junctions from short-read alignments. The `-s` option accepts either `sam` or `bam` files, and if there are multiple sams/bams they can be provided in a comma-separated list.

**Usage:**
```sh
python junctions_from_sam.py -s <shortreads.sam>|<shortreads.bam> -n outname
```
The file that can be supplied to flair-correct with `-j` is in the output file `outname_junctions.bed`. It is recommended that the user remove infrequently used junctions i.e. junctions with few supporting junction reads, which are in the 5th column of the junction bed file. For example, if you wanted to do the filter out junctions with fewer than 3 short reads, you could use `awk '{ if ($5 >= 3) { print } }' outname.sj_junctions.bed > outname.filtered.bed`.

Alternatively, the `-j` argument for flair-correct can also be generated using STAR. STAR 2-pass alignment of short reads produces a compatible splice junction file (`SJ.out.tab`). We recommend filtering out junctions with few uniquely mapping reads (column 7).


### <a name="collapse"></a>flair collapse
Defines high-confidence isoforms from corrected reads. As FLAIR does not use annotations to collapse isoforms, FLAIR will pick the name of a read that shares the same splice junction chain as the isoform to be the isoform name. It is recommended to still provide an annotation with `-f`, which is used to rename FLAIR isoforms that match isoforms in existing annotation according to the transcript_id field in the gtf. Intermediate files generated by this step are removed by default, but can be retained for debugging purposes by supplying the argument `--keep_intermediate` and optionally supplying a directory to keep those files with `--temp_dir`.

If there are multiple samples to be compared, the flair-corrected read `psl` files should be concatenated prior to running flair-collapse. In addition, all raw read fastq/fasta files should either be specified after `-r` with space/comma separators or concatenated into a single file.

**Usage:**
```sh
python flair.py collapse -g genome.fa -r <reads.fq>|<reads.fa> -q <query.psl>|<query.bed> [options]
```
run with `--help` for description of optional arguments.

Outputs the high-confidence isoforms in several formats: (1) `*isoforms.psl` or `.bed`, (2) `*isoforms.gtf`, as well as (3) an `*isoforms.fa` file of isoform sequences. If an annotation file is provided, the isoforms ID format will contain the transcript id, underscore, and then the gene id, like so: `EN*T*_EN*G*`. If multiple TSSs/TESs are allowed (toggle with `--max_ends` or `--no_redundant`), then a `-1` or higher will be appended to the end of the isoform name for the isoforms that have identical splice junction chains and differ only by their TSS/TES. For the gene field, the gene that is assigned to the isoform is based on whichever annotated gene has the greatest number of splice junctions shared with the isoform. If there are no genes in the annotation which can be assigned to the isoform, a genomic coordinate is used (e.g. `chr*:100000`.

### <a name="quantify"></a>flair quantify
Convenience function to quantifying FLAIR isoform usage across samples using minimap2. If isoform quantification in TPM is desired, please use the `--tpm` option. If the user prefers [salmon](https://combine-lab.github.io/salmon/getting_started/) to quantify transcripts using their nanopore reads, please specify a path to salmon using `--salmon`. For all options run flair-quantify with `--help`.

**Usage:**
```sh
python flair.py quantify -r reads_manifest.tsv -i isoforms.fasta [options]
```

**Inputs:**</br>
(1) `reads_manifest.tsv` is a tab-delimited file containing the sample name, condition, batch\*, and path to reads.fq/fa.
For exmaple:
```tsv
sample1	conditionA	batch1	./sample1_reads.fq
sample2	conditionA	batch1	./sample2_reads.fq
sample3	conditionA	batch2	./sample3_reads.fq
sample4	conditionB	batch1	./sample4_reads.fq
sample5	conditionB	batch1	./sample5_reads.fq
sample6	conditionB	batch2	./sample6_reads.fq
```
\* The batch descriptor is used in the downstream flair-diffExp analysis to model unintended variability due to secondary factors such as batch or sequencing replicate. If unsure about this option, leave this column defined as `batch1` for all samples.

(2) `isoforms.fasta` contains FLAIR collapsed isoforms produced by the [`flair-collapse`](#collapse) module.

**Outputs:**</br>
(1) `counts_matrix.tsv` which is a tab-delimited file containing isoform counts for each sample. In the output, the values in the manifest file are concatenated with underscores so please do not use underscores in the manifest file. For example:

```tsv
ids	samp1_conditionA_batch1	samp2_conditionA_batch1 samp3_conditionA_batch2	...
0042c9e7-b993_ENSG00000131368.3	237.0	156.0	165.0	150.0	...
0042d216-6b08_ENSG00000101940.13	32.0	14.0 	25.0	...
```


### <a name="diffExp"></a>flair diffExp
Performs differential isoform expression, differential gene expression, and differential isoform usage analyses between multiple samples with 3 or more replicates. For differential isoform usage analysis between samples without replicates, please use the [diff_iso_usage.py](#diffisoscript) standalone script. This module requires additional python modules and R packages which are described below: 

#### Additional Requirements
1. python v2.7+ and python modules: pandas, numpy, rpy2
2. [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html)
3. [ggplot2](https://ggplot2.tidyverse.org)
4. [qqman](https://cran.r-project.org/web/packages/qqman/index.html)
5. [DRIMSeq](http://bioconductor.org/packages/release/bioc/html/DRIMSeq.html)
6. [stageR](http://bioconductor.org/packages/release/bioc/html/stageR.html)

**Usage:**
```sh
python flair.py diffExp -q counts_matrix.tsv -o output_directory [options]
```

**Inputs:**</br>
(1) `counts_matrix.tsv` is a tab-delimited file generated by the [`flair-quantify`](#quantify) module.

**Outputs:**</br>
(1) Files contained in the `output_directory` are tables and plots generated from the various R-packages used in this analysis, including raw DESeq2/DRIMSeq output tables with foldChange, isoform frequency and adjusted pvalues. 


### <a name="diffsplice"></a>flair diffSplice
Calls alternative splicing events from isoforms. Currently we support the following AS events: intron retention, alternative 3' splicing, alternative 5' splicing, and cassette exons. 

**Usage:**
```sh
python flair.py diffSplice -i <isoforms.bed>|<isoforms.psl> -q counts_matrix.tsv [options]
```
If there are 3 or more samples per condition, then the user can run with `--test` and DRIMSeq will be used to calculate differential usage of the alternative splicing events between two conditions. Run with `--help` to see more DRIMSeq-specific arguments. If conditions were sequenced without replicates, then the flair-diffsplice output files can be input to the [diffsplice_fishers_exact.py](#diffsplice_fishers) script for statistical testing instead.


**Inputs:**</br>
(1) `-i` is a tab-delimited isoforms.bed/psl file generated by the [`flair-collapse`](#collapse) module.<br>
(2) `-q` is a tab-delimited counts_matrix.tsv file generated by the [`flair-quantify`](#quantify) module.

**Outputs:**</br>
(1-4) 4 tab-delimited files for each AS event type. If DRIMSeq was run, and a file with PSIs for each sample and the corresponding p-values for each event type (5-8).

For a complex splicing example, please note the 2 alternative 3' SS, 3 intron retention, and 4 exon skipping events in the following set of isoforms that flair-diffSplice would call and the isoforms that are considered to include or exclude the each event:</br>
<img src='https://people.ucsc.edu/~atang14/flair/toy.isoforms.coord.png' alt='isoforms' width='780'/></br>
```tsv
a3ss_feature_id		coordinate			sample1	sample2	...	isoform_ids
inclusion_chr1:80	chr1:80-400_chr1:80-450		75.0	35.0	...	a,e
exclusion_chr1:80	chr1:80-400_chr1:80-450		3.0	13.0	...	c
inclusion_chr1:500	chr1:500-650_chr1:500-700	4.0	18.0	...	d
exclusion_chr1:500	chr1:500-650_chr1:500-700	70.0	17.0	...	e
```
```tsv
ir_feature_id		coordinate	sample1	sample2	...	isoform_ids
inclusion_chr1:500-650	chr1:500-650	46.0	13.0	...	g
exclusion_chr1:500-650	chr1:500-650	4.0	18.0	...	d
inclusion_chr1:500-700	chr1:500-700	46.0	13.0	...	g
exclusion_chr1:500-700	chr1:500-700	70.0	17.0	...	e
inclusion_chr1:250-450	chr1:250-450	50.0	31.0	...	d,g
exclusion_chr1:250-450	chr1:250-450	80.0	17.0	...	b
```
```tsv
es_feature_id		coordinate	sample1	sample2	...	isoform_ids
inclusion_chr1:450-500	chr1:450-500	83.0	30.0	...	b,c
exclusion_chr1:450-500	chr1:450-500	56.0	15.0	...	f
inclusion_chr1:200-250	chr1:200-250	80.0	17.0	...	b
exclusion_chr1:200-250	chr1:200-250	3.0	13.0	...	c
inclusion_chr1:200-500	chr1:200-500	4.0	18.0	...	d
exclusion_chr1:200-500	chr1:200-500	22.0	15.0	...	h
inclusion_chr1:400-500	chr1:400-500	75.0	35.0	...	e,a
exclusion_chr1:400-500	chr1:400-500	56.0	15.0	...	f
```


## <a name="scripts"></a>Scripts

We have also provided standalone scripts for splicing and productivity analysis of quantified isoforms from flair-collapse output.

### predictProductivity.py

Annotated start codons from the annotation are used to identify the longest ORF for each isoform for predicting isoform productivity. Requires three arguments to classify isoforms according to productivity: (1) isoforms in `psl` or `bed` format, (2) `gtf` genome annotation, (3) `fasta` genome sequences. Bedtools must be in your $PATH for predictProductivity.py to run properly.

**Usage:**
```sh
python predictProductivity.py -i <isoforms.bed>|<isoforms.psl> -g annotation.gtf -f genome.fa --longestORF > productivity.bed
```
Outputs a bed file with either the values `PRO` (productive), `PTC` (premature termination codon, i.e. unproductive), `NGO` (no start codon), or `NST` (has start codon but no stop codon) appended to the end of the isoform name. When isoforms are visualized in the UCSC genome browser or IGV, the isoforms will be colored accordingly and have thicker exons to denote the coding region.

### mark_intron_retention.py

Requires three positional arguments to identify intron retentions in isoforms: (1) a `psl` of isoforms, (2) `psl` output filename, (3) `txt` output filename for coordinates of introns found.

**Usage:**
```sh
python mark_intron_retention.py <isoforms.psl>|<isoforms.bed> out_isoforms.psl out_coords.txt
```
Outputs (1) an extended `psl` with an additional column containing either values 0 or 1 classifying the isoform as either spliced or intron-retaining, respectively; (2) `txt` file of intron retentions with format `isoform name` `chromosome` `intron 5' coordinate` `intron 3' coordinate`. Note: A psl or bed file with more additional columns will not be displayed in the genome browser, but can be displayed in IGV.

### <a name="diffisoscript"></a>diff_iso_usage.py
Requires four positional arguments to identify and calculate significance of alternative isoform usage between two samples using Fisher's exact tests: (1) counts_matrix.tsv from flair-quantify, (2) the name of the column of the first sample, (3) the name of the column of the second sample, (4) `txt` output filename containing the p-value associated with differential isoform usage for each isoform. The more differentially used the isoforms are between the first and second condition, the lower the p-value.

**Usage:**
```sh
python diff_iso_usage.py counts_matrix.tsv colname1 colname2 diff_isos.txt
```
Output file format columns are as follows: 
`gene name` `isoform name` `p-value` `sample1 isoform count` `sample2 isoform count` `sample1 alternative isoforms for gene count` `sample2 alternative isoforms for gene count` 

### plot_isoform_usage.py
Visualization script for FLAIR isoform structures and the percent usage of each isoform in each sample for a given gene. If you supply the isoforms.bed file from running `predictProductivity.py`, then isoforms will be filled according to the predicted productivity (solid for `PRO`, hatched for `PTC`, faded for `NGO` or `NST`). The gene name supplied should correspond to a gene name in your isoform file and counts file.

**Usage:**
```sh
python plot_isoform_usage.py <isoforms.psl>|<isoforms.bed> counts_matrix.tsv gene_name 
```

Outputs (1) gene_name_isoforms.png of isoform structures and (2) gene_name_usage.png of isoform usage by sample. 

For example:

<img src='https://people.ucsc.edu/~atang14/flair/toy_diu_isoforms.png' alt='isoforms' width='480'/></br>
<img src='https://people.ucsc.edu/~atang14/flair/toy_diu_usage.png' alt='usage' width='480'/>

### <a name="diffsplice_fishers"></a>diffsplice_fishers_exact.py

Identifies and calculates the significance of alternative splicing events between two samples without replicates using Fisher's exact tests. Requires four positional arguments: (1) flair-diffSplice `tsv` of alternative splicing calls for a splicing event type, (2) the name of the column of the first sample, (3) the name of the column of the second sample, and (4) `tsv` output filename containing the p-values from Fisher's exact tests of each event. 

**Usage:**
```sh
python diffsplice_fishers_exact.py events.quant.tsv colname1 colname2 out.fishers.tsv 
```

The output file contains the original columns with an additional column containing the p-values appended.


## Docker <a name="docker"></a>
If the user wishes to run FLAIR using [docker](https://docs.docker.com/) instead of cloning this repository, the following commands can be used:
`docker pull quay.io/brookslab/flair`
`docker run -w /usr/data -v [your_path_to_data]:/usr/data  -t -d [image_id]`
`docker exec [container_id] python3 /usr/local/flair/flair.py [your_command]`

## Conda Environment <a name="condaenv"></a>
Users can run FLAIR within the conda environment provided in `misc/flair_conda_env.yaml`. FLAIR should run smoothly in this environment. Refer to the conda docs for how to [create an environment from an environment.yml file](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file).

## Example Files <a name="exfiles"></a>
We have provided the following example files [here](https://people.ucsc.edu/~atang14/flair/example_files/):  
- `star.firstpass.gm12878.junctions.3.tab`, a file of splice junctions observed from short read sequencing of GM18278 that can be used in the correction step with `-j`. Junctions with fewer than 3 uniquely mapping reads have been filtered out
- `promoter.gencode.v27.20.bed`, promoter regions determined from [ENCODE promoter chromatin states for GM12878](http://hgdownload.cse.ucsc.edu/goldenPath/hg18/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmGm12878HMM.bed.gz) and 20 bp around annotated TSS in GENCODE v27. Can be supplied to flair-collapse with `-p` to build the initial firstpass set with only reads with start positions falling within these regions

Other downloads:
- [Native RNA Pass reads](https://github.com/nanopore-wgs-consortium/NA12878/blob/master/RNA.md) Running these 10 million nanopore reads from `fastq` through flair align, correct, and collapse modules to assembled isoforms with 8 threads requires ~3.5 hours (includes ~2.5 hours of minimap2 alignment)
- [NanoSim_Wrapper.py](https://github.com/BrooksLabUCSC/labtools/blob/master/NanoSim_Wrapper.py), a wrapper script written for simulating nanopore transcriptome data using [Nanosim](https://github.com/bcgsc/NanoSim)

## Cite FLAIR <a name="cite"></a>
If you use or discuss FLAIR, please cite the following [paper](https://www.nature.com/articles/s41467-020-15171-6):
>Tang, A.D., Soulette, C.M., van Baren, M.J. et al. Full-length transcript characterization of SF3B1 mutation in chronic lymphocytic leukemia reveals downregulation of retained introns. Nat Commun 11, 1438 (2020). https://doi.org/10.1038/s41467-020-15171-6

