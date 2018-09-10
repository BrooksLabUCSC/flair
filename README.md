# flair
FLAIR (Full-Length Alternative Isoform analysis of RNA) for the alignment, correction, isoform definition, and alternative splicing analysis of noisy reads. 

## Table of Contents

- [Overview](#overview)
- [Requirements](#requirements)
- [FlAIR modules](#modules)
	- [align](#align)
	- [correct](#correct)
	- [collapse](#collapse)
- [Scripts](#scripts)

## <a name="overview"></a>Overview
FLAIR can be run optionally with short-read data to help splice site accuracy of the long read splice junctions. FLAIR uses multiple alignment steps and splice site filters to increase confidence in the set of isoforms defined from noisy data. FLAIR was designed to be able to sense subtle splicing changes in nanopore data from [Tang et al. (2018)](https://www.biorxiv.org/content/early/2018/09/06/410183). Please read for more description of some methods.
![flair_workflow](misc/flair_workflow.png) 
It is recommended to combine all samples together prior to running FLAIR modules for isoform assembly, followed by individual sample read assignment to isoforms of the combined assembly.

## <a name="requirements"></a>Requirements

1. python v2.7+
2. [mRNAtoGene](https://github.com/ENCODE-DCC/kentUtils/tree/master/src/hg/mrnaToGene) in $PATH
3. [minimap2](https://github.com/lh3/minimap2)

## <a name="modules"></a>flair modules 
flair.py is a wrapper script with modules for running various processing scripts located in bin/. Modules are assumed to be run in order (align, correct, collapse), but the user can forgo the wrapper if a more custom build is desired. 

### <a name="align"></a>flair align
Aligns reads to the genome using minimap2, and converts `.sam` output to [.psl](https://genome.ucsc.edu/FAQ/FAQformat.html#format2), the predominant format used in consequent steps. Aligned reads in PSL format can be visualized in IGV; alternatively, the UCSC Genome browser can also be used if a chromosome sizes tab-separated file is provided with `-c`.

Usage:
```sh
python flair.py align -r <reads.fq>/<reads.fa> -g genome.fa [options]
```
run with `--help` for a description of optional arguments.

### flair correct

Smooths gaps and corrects misaligned splice sites using genome annotations. To use short-read splice sites to aid with correction, use [junctionsFromSam.py](https://github.com/BrooksLabUCSC/labtools/blob/master/junctionsFromSam.py) to extract splice junctions.

Usage:
```sh
python flair.py correct -a annotated.gp -g genome.fa -q query.psl [options]
```
run with `--help` for description of optional arguments.
Outputs (1) PSL of raw reads with strand and gene inferred and (2) PSL of corrected reads within directory specified by `-o`

### flair collapse
Defines isoforms from correct reads. If a GTF is provided with `-f`, isoforms that match isoforms in existing annotation will be named using the Ensembl ID in existing annotation. By default, redundant isoforms (those that are proper subsets of another isoform in the set) are filtered out, an option that can be toggled with `-e`. Isoforms in PSL format can be visualized again in IGV, or the UCSC genome browser if columns 22, number of supporting reads, is removed. 

Usage:
```sh
python flair.py collapse -r <reads.fq>/<reads.fa> -q <query.psl>/<query.bed12> -g genome.fa [options]
```
run with `--help` for description of optional arguments.
Outputs (1) extended PSL containing the data-specific isoforms and (2) fasta file of isoform sequences.

## Scripts

### find_alt3prime_5prime_ss.py

Requires four positional arguments to identify and calculate significance of alternative 3' and 5' splicing between two samples using Fisher's exact tests. First, an extended PSL of isoforms containing two extra columns for read counts of each isoform per sample type; second, the column number of the two extra columns (assumed to be last two); third, txt file output name for alternative 3' SS; fourth, txt file output name for alternative 5' SS.
Usage: 
```sh
python find_alt3prime_5prime_ss.py isoforms.psl colnum alt3.txt alt5.txt 
```
Output file format:
`chrom` `intron 3'` `intron 5'` `p-value` `strand` `sample1 intron count` `sample2 intron count` `sample1 alternative introns counts` `sample2 alternative introns counts` `gene name` `SS distance from predominant alternative SS` `predominant alternative SS`

### diff_iso_usage.py
Requires three positional arguments to identify and calculate significance of alternative 3' and 5' splicing between two samples using Fisher's exact tests. First, an extended PSL of isoforms containing two extra columns for read counts of each isoform per sample type; second, the column number of the two extra columns (assumed to be last two); third, txt file output name for differentially used isoforms.
```sh
python diff_iso_usage.py isoforms.psl colnum diff_isos.txt
```
Output file format: 
`gene name` `isoform name` `p-value` `sample1 isoform count` `sample2 isoform count` `sample1 alternative isoforms for gene count` `sample2 alternative isoforms for gene count` 


