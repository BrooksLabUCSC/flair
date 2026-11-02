# Major user-visible changes

## [unreleased]
* General
  * supported Python versions are 3.12 through 3.14 in both the pip and Conda
    installs.  pysam is now 0.24, the first release with wheels above 3.12
  * the `--help` tests no longer fail on Python 3.12.  Their expected files are
    in the 3.13 argparse format, and the comparison is skipped, naming the
    python version, on anything older
  * reverse complement takes every IUPAC nucleotide code, so a genome with
    degenerate bases no longer stops FLAIR fusion with `KeyError: 'R'` (#715).
    A character that is not a nucleotide now reports itself
  * dropped the unused setuptools dependency, which was the last one carrying
    an open advisory (GHSA-h35f-9h28-mq5c)
  * updated certifi, pillow, urllib3 and httpx2 in the lock file
  * the conda environments carry `zlib`.  mappy and ncls publish no linux wheel
    above python 3.12, so pip builds them from source, and mappy needs the zlib
    headers, which the `libzlib` that came in by itself does not have
* FLAIR fusion
  * keep the 500 base gene padding inside the sequence.  A gene starting within
    500 bases of a sequence start failed with `ValueError: start out of range`,
    and the far end truncated silently, leaving coordinates that claimed a
    longer span than the sequence held (#674)
  * error out naming the sequence and both files when the annotation names a
    sequence the genome does not have, rather than failing later
  * an exon starting exactly on a synthetic locus boundary is no longer lost,
    and an isoform whose exons cannot all be placed in one locus is skipped
    with a message rather than written with the wrong blocks
* Documentation
  * corrected the v3.0.0 release date in this file to 2026-01-30
  * installing.rst no longer shows a release version for a `+master` tree

## [v3.0.1] 2026-08-10
* General
  * fixed problems with some diffsplice_fishers_exact, and other
    auxiliary console script installs
  * dropped the `mark_productivity` console script and the `bin/bam2Bed12`
    and `bin/mark_productivity` links, which were never functional
  * supported Python versions widened from 3.12 only to 3.12 through 3.14,
    in both the pip and Conda installs.  The expected files for the `--help`
    tests are now generated with Python 3.13, so those tests fail on 3.12
* FLAIR collapse
  * worked around `--remove_internal_priming` failing when applied during the
    transcriptome alignment phase (#629)
* FLAIR combine
  * error out with a clear message when the manifest contains no samples
  * added debug logging of sample, intron-chain and final isoform counts
* Documentation
  * new FLAIR pipeline figure on the main page; the compartmentalized workflow
    figure moved to the modules page
  * fixed figures mangled by the conversion to SVG and a bad margin
  * removed the obsolete FLAIR2 functions page
  * removed the stale note that the protocol applies only to a GitHub master
    checkout

## [v3.0.0] 2026-01-30
* General
  * Bug fixes since v3.0.0b1

## [v3.0.0b1] 2025-11-27 (beta 1)
* General
  * Removed dependencies on pandas and rpy2.
  * Added proper logging throughout
* FLAIR align
  * Cleaned up output files - now outputs unfiltered BAM and filtered BED
* FLAIR correct
  * Changed orthogonal junction inputs - now specify whether your file is 
  a bed (--junction_bed) or whether it comes from STAR short-read RNA alignment (--junction_tab)
  * Now support detecting junctions directly from long read data: 
  run intronProspector to generate junction bed, input that via --junction_bed, 
  specify desired junction read support with --junction_support
  * Removed -g genome option, increases speed
  * Now produces sorted output files
* FLAIR collapse
  * cleaned up output files, now only gives: isoforms.bed, isoforms.gtf, isoforms.fa, isoform.read.map.txt
  * (experimental) added filter for removing reads with internal priming. 
  Options: --remove_internal_priming, --intprimingthreshold, --intprimingfracAs
  * allow CDS prediction in collapse with --predictCDS
  * removed genomic range option (too fragile). For parallelization, run FLAIR transcriptome
  * made gtf (annotation) parsing more robust, especially for unsorted files 
  * improved fractional support filtering (with --support < 1)
  * improved isoform haplotyping through longshot - this is being deprecated though, please use FLAIR variants
* FLAIR quantify
  * improved fraction of reads assigned to isoforms and recovered (better handling of ambiguous alignments)
  * improved processing of multiple samples so running requires less available memory and is faster
* FLAIR diffexp and diffsplice
  * recoded directly in R to improve performance
* New modules
  * FLAIR transcriptome
    * Combines the functions of correct and collapse
    * Runs directly from an aligned BAM file
    * Performs more effective parallelization (specify in --parallelmode)
  * FLAIR fusion
    * Detects gene fusions and fusion isoforms
    * Fusion detection accuracy is comparable to JAFFAL and CTAT-lr-fusion
    * Fusion isoform detection is much more accurate
  * FLAIR combine
    * Allows combining transcriptomes generated from different samples
    * Fusion isoforms can be combined with collapsed isoforms to form a full transcriptome
    * Manipulate filtering of single exon isoforms with --include_se
  * FLAIR variants
    * Allows detection of variant-aware transcripts through 
    identification of read clusters with shared variants
    * Uses variants detected either from WGS or from lr-RNA-seq. 
    We recommend Longshot.
    * Good for identifying splicing of specific variants
    * Can be used in conjunction with FLAIR diffexp or diff_iso_usage 
    to identify changes in variant-aware transcripts between samples/groups
* Added FLAIR protocol to improve guidance on running flair


## [v2.2.1] 2025-04-06
* Fixed case where --quality=0 default didn't work, as it tested for greater
  than the cutoff.
* Converted all Python code to use 4-space indentation.  Previously some code
  used tabs, other used 4-spaces.
* Added missing dependencies for plot_isoform_usage.


## [v2.2.0] 2025-05-06
* Returned diffexp and diffsplice as standard modules.  The BioConda
  environment does not include the dependencies for these modules
  and required software does not run on Apple Silicon (ARM64) systems.
* The flair combine functionality is now to a module.  The
  `flair_combine` program is run with `flair combine`.
* Changed default MAPQ minimum quality score to 0. This allows more reads to
  be used in identifying isoforms, which tends to improve the overall models
  with out adversely affecting the accuracy.
* GitHub releases include the Conda YAML files for building FLAIR
  environments.  Useful if the BioConda release has not been manually
  reviewed.
* The FLAIR Docker now includes all dependencies to run diffexp and diffsplice.
* Reorganized the installation documentation.
* Fixed flair_quantify --output_bam crash
* Other bug fixes

## [v2.1.2] 2025-04-17
* Address issue getting BioConda to work
* Bug fixes for collapse command line parsing


## [v2.1.1] 2025-04-10
* converted all programs to use console scripts to allow BioConda to work


## [v2.1.0] 2025-03-27
* Numerous bug fixes.
* Removed support for PSL format.
* Remove `flair 123` to run multiple modules at once.
* Compatibility with Python 3.12 
* Compatibility with Apple ARM64 systems.
* Deprecated `diffExp` and `diffSplice`, they will be removed in a future release.
  Lets us know if you use this functionality.  Their dependencies are no longer
  part of the conda package, they can be added the conda environment with
  `misc/flair_diffexp_conda_env.yaml`.
