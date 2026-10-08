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
  * Bad input is reported as one line rather than a Python traceback: every
    installed command now reports errors this way, and the errors for mismatched
    annotation fasta and isoform names, alignments without AS or MD tags, and the
    option and input problems in junctions_from_sam and identify_vars say what to
    do about them.
  * `flair diffexp` gained `--condition_a` and `--condition_b`, matching
    `flair diffsplice`.  `condition_a` is the reference that fold changes are
    measured against.
  * The `diffexp` and `diffsplice` R dependencies are now part of
    `misc/flair_conda_env.yaml`, and `misc/flair_diffexp_conda_env.yaml` is gone.
    The BioConda package still does not carry them; `installing.rst` gives the
    `conda install` command to add them.
* flair transcriptome
  * About 3x faster: full WTC11 PacBio (3 million reads, GENCODE, 24 threads)
    went from 75 to 23 minutes and from 1,490 to 500 CPU minutes.  Setup before
    the partitions run is under 2 minutes (was 11); minimap2 threads stay within
    `--threads`, and the largest partitions run first.  Reads that match a
    first-pass isoform exactly are assigned without the final realignment, and
    reads no isoform can take are not realigned, which cuts the final
    realignment's minimap2 time 6-fold.  Reads are aligned to the annotation only
    when it can change their correction.
  * Isoform ends come from the reads assigned to each isoform in the final
    realignment, against transcripts padded so that reads running past an
    isoform's end are no longer rejected: more reads are assigned to isoforms.
    The densest cluster of read ends gives its mode.
  * Single-exon isoforms are reported at each clustered pair of read ends,
    whatever `--max_ends` is; single-exon isoforms that are 3' UTR fragments of a
    spliced transcript, and internally primed single-exon reads (an A-rich aligned
    end, or SQANTI3's 60% A after an untailed end), are removed.
  * A spliced isoform removed as a subset is reported after all when every
    isoform it is a subset of fails support and it has at least
    `--subset_backup_support` reads (default 10).
  * Isoforms are assigned to genes by junctions and splice sites, and novel genes
    are grouped by shared splicing.
  * Semi-canonical junction motifs (GC-AG, AT-AC) are weak strand evidence:
    without an annotation, the strand of a read whose junctions all have them
    comes from gene identification.  `--trust_junctions` trusts every input
    junction's strand instead.
  * `--trust_strand` applies to spliced reads too; `--directRNA` declares direct
    RNA input (it implies `--trust_strand` and allows no large UTR deletions);
    `--total_rna` keeps single-exon reads that end at genomic A runs and needs no
    poly(A) tails (it requires `--trust_strand`).
  * Temporary files go in `--temp_dir`, by default $TMPDIR.
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
  * `flair diffexp` no longer takes the two conditions from the first and last
    column of the counts matrix.  With `--condition_a` and `--condition_b` left
    out, the two conditions are used in sorted order.  For a counts matrix whose
    first and last columns are not the conditions in sorted order, the reference
    changes, which flips the sign of the reported fold changes and renames the
    output files `prefix_X_v_Y.tsv`.  Results are not comparable across this
    change; name the conditions to get a specific direction.
  * `flair diffexp` previously took the two conditions from the first and last
    column even when they were the same condition, as in a column order A,B,B,A,
    which filtered one condition twice and ignored the other.  Such runs gave
    wrong results and now give correct ones.
  * `diff_iso_usage` output columns are named after the two samples given rather
    than sample1 and sample2.
  * `flair quantify` writes `<output>.sample_info.tsv` beside the counts file,
    giving the condition and batch of each counts column, and the counts columns
    are now named after the sample alone rather than `sample_condition_batch`.
    `flair diffexp` and `flair diffsplice` read that file.  A counts matrix from
    an earlier FLAIR is still read by taking the fields out of the column names.
  * Because condition and batch are columns of their own, the id, condition and
    batch fields of the quantify manifest may now contain underscores.
  * `flair quantify --sample_id_only` is gone; the counts columns always name the
    sample alone, which is what it asked for.
  * `flair transcriptome --trust_ends` is removed: ends now come from the reads
    assigned to each isoform.
  * `flair transcriptome --filter bysupport` is replaced by `--filter <N>X`, which
    keeps a subset isoform with more than N times its supersets' reads;
    `bysupport` was `1.2X`.  The default is still `nosubset`.
  * `flair transcriptome --normalize_ends` now reports, for each spliced isoform,
    the ends it shares with the isoforms of its gene whose terminal splice sites
    are near its own.
  * `flair transcriptome --ss_window` defaults to 10 (was 15).
  * `count_sam_transcripts` no longer takes `--quality`, `--fusion_dist` or
    `--end_norm_dist`.
    
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
