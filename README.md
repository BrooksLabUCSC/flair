# flair
FLAIR (Full-Length Alternative Isoform analysis of RNA) for the correction, isoform definition, and alternative splicing analysis of noisy reads. FLAIR has primarily been used for nanopore cDNA, native RNA, and PacBio sequencing reads.

The complete Flair manual is available via [readthedocs](https://flair.readthedocs.io/en/latest/)

## FLAIR Announce Mailing list
**If you are using FLAIR, please subscribe to the FLAIR mailing list:**

<https://groups.google.com/a/ucsc.edu/g/flair-announce-group>

This is a read-only, low volume list that will only have announcement of new
FLAIR releases, publications and other FLAIR-related user information.


# FLAIR3 usage 
variant + fusion-aware isoform detection + AA-seq prediction

In misc there is a file for running this whole pipeline if your files are named consistently with sequential numbers

## For each individual sample:

- FLAIR align with options: -f separate --minfragmentsize 40
- FLAIR correct standard
- FLAIR collapse with options: --generate_map --annotation_reliant generate --stringent --check_splice --quality 0 --isoformtss -n longest
- flair_detectfusions.py with options: -b sample.flair.align_chimeric.bam --annotated_fa sample.flair.collapse.annotated_transcripts.fa
    - uses files from FLAIR align and collapse as input

## For combining all samples and doing futher analysis:

- Make combine manifest (tab separated) with two lines per sample (one line for normal isoforms, one for fusion isoforms) as follows:
    - sampleid, [isoform OR fusionisoform], sample.isoforms.bed, sample.isoforms.fa, sample.combined.isoform.read.map.txt
- FLAIR combine with default parameters
- Make quantify manifest (tab separated) with one line per sample as follows:
    - sampleid, condition, batch, sample.reads.fa
- FLAIR quantify with options: -i experiment.combined.isoforms.fa --output_bam --quality 0
- Make variant calling manifest file (tab separated) with one line per sample as follows:
    - sampleid, sample.flair.quantify.aligned.bam
- flair_variants.py with normal inputs
    - this requires longshot to be installed and runnable in your path

## Cite FLAIR <a name="cite"></a>
If you use or discuss FLAIR, please cite the following [paper](https://www.nature.com/articles/s41467-020-15171-6):
>Tang, A.D., Soulette, C.M., van Baren, M.J. et al. Full-length transcript characterization of SF3B1 mutation in chronic lymphocytic leukemia reveals downregulation of retained introns. Nat Commun 11, 1438 (2020). https://doi.org/10.1038/s41467-020-15171-6

