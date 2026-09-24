# Major user-visible changes

## [v3.x.xx]
* General
  * fixed problems with some diffsplice_fishers_exact, and other
    auxiliary console script installs
  * flair now uses real subcommands, so `flair <subcommand> --help` and the
    subcommand list in `flair --help` describe what each one accepts.
  * Command line documentation is generated from the programs themselves, so it
    can no longer drift from what they accept.  flair variantquant, alleles and
    isoalleles are documented for the first time.
  * All subcommands now honor the logging options; previously only flair align did.
* Incompatibles
  * Removed flair correct and collapse modules, the functionality is replaced
    by flair transcriptome.
  * Remove flair align options that are no longer need without flair correct.
  * Logging options (--log-level, --log-stderr, --log-conf, --log-debug) and
    --version must now precede the subcommand: `flair --log-debug align ...`.
  * Options renamed for one name per concept and one meaning per short option.
    Old names are not accepted.
    * everywhere: --out_dir, --output_prefix to --output; --bed_paths to
      --isoform_beds; --bedisoforms and diffsplice --isoforms to --isoform_bed;
      --out_dir_force to --overwrite_output; --norm_ends to --normalize_ends
    * align: --filtertype to --filter_type; --minfragmentsize to
      --min_fragment_size; --maxintronlen to --max_intron_len; --nvrna to
      --native_rna
    * transcriptome: --se_support to --single_exon_support; --keep_sup to
      --keep_supplementary
    * combine: --endwindow to --end_window; --minpercentusage to
      --min_percent_usage; --remove_se to --remove_single_exon
    * fusion: --support to --min_support; --maxloci to --max_loci;
      --max_dist_to_TSS to --max_dist_to_tss; --min_dist_between_bp to
      --min_dist_between_breakpoints
    * diffexp: --exp_thresh to --min_expression
    * diffsplice: --drim1 to --drim4 to --min_samps_gene_expr,
      --min_samps_feature_expr, --min_gene_expr, --min_feature_expr, after the
      DRIMSeq dmFilter arguments they are passed to; --conditionA and
      --conditionB to --condition_a and --condition_b
    * variantquant: --input_bam to --transcriptome_bam; --threshold to
      --min_coverage
    * alleles and isoalleles: --bam to --tumor_bam; --norm_bam to --normal_bam;
      --norm_vcf to --normal_vcf; --allele_read_map_norm and --iso_read_map_norm
      to --allele_read_map_normal and --iso_read_map_normal
  * Short options dropped where a letter meant two different things: -w, -p, -q,
    -e, -s, -i, -k, -m, -v and -of.  -b is now always the aligned BAM, -r the
    reads, -t threads, -g the genome, -f the annotation and -o the output.
  * Removed options that were accepted and never used: flair align --quiet and
    flair quantify --quality.
  * flair align now requires exactly one of --genome and --mm_index, rather than
    failing later when neither was given.
  * junctions_from_sam, fasta_seq_lengths, mark_intron_retention,
    identify_annotated_gene and diffsplice_fishers_exact now use standard option
    parsing and accept --help; bed_to_gtf --noCDS is now --no_cds.
  * annotate_aaseq_with_uniprot is now installed as a command.
    
## [v3.0.0] 2025-11-31
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
