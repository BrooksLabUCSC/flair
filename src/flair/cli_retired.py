"""The subcommand and option names earlier releases accepted, with what replaced
each one.

The command line is checked against these tables before it is parsed, rather than
the names being hidden options of the parsers.  A hidden option is only reached
once the rest of the command line parses, and a command line being converted from
an earlier release usually has several old names and not yet the arguments the new
ones require.  Hidden long names also make the abbreviations of the names that
remain ambiguous: with --out_dir in the parser, --out no longer selects --output.
"""
import sys
import itertools
import argparse
from flair import FlairError

# names dropped in every subcommand
RETIRED_EVERYWHERE = {
    '--out_dir': "is now --output, a file name base rather than a directory",
    '--output_prefix': "is now --output",
    '--out_dir_force': "has been removed: --output names a file name base, and its outputs are replaced",
    '--overwrite_output': "has been removed: --output names a file name base, and its outputs are replaced",
    '-of': "has been removed: --output names a file name base, and its outputs are replaced",
    '--bed_paths': "is now --isoform_beds",
    '--bedisoforms': "is now --isoform_bed",
    '--isoforms': "is now --isoform_bed",
    '--norm_ends': "is now --normalize_ends",
    '--outDir': "is now --output, a file name base rather than a directory",
    '--quiet': ("has been removed: the logging options --log-level, --log-stderr, --log-conf and "
                "--log-debug set what is logged, and come before the subcommand"),
}

_STRINGENT = ("has been removed: a read must be full length on an isoform to count toward it, "
              "which --stringent asked for")
_CHECK_SPLICE = ("has been removed: a read must align accurately at the splice sites to count "
                 "toward an isoform, which --check_splice asked for")
_NO_STRINGENT = ("has been removed, with {other}: they counted reads that do not support an "
                 "isoform toward it, the first by dropping the requirement that the read be full "
                 "length on it, the second by dropping the check that it aligns accurately at the "
                 "splice sites.  Drop them, or count the reads against an isoform set with flair "
                 "quantify")
_ENDS_FROM_READS = ("has been removed: an isoform's ends come from the ends of the reads assigned "
                    "to it, clustered within --end_window, at most --max_ends of them, or "
                    "--normalize_ends")
_ALIGN_TO_ANNOT = ("has been removed: reads are aligned to the transcripts of -f by default, with "
                   "no FASTA or BED of them to build or pass, and --no_align_to_annot turns that "
                   "alignment off")
_INTERNAL_PRIMING = ("has been removed: a single exon read whose 3' end is in A-rich genome is "
                     "always removed as oligo-dT mispriming, and --keep_internally_primed keeps it")
_PREDICT_CDS = ("has been removed: the CDS and productivity of an isoform are reported whenever "
                "the annotation given with -f has start codons")
_SUPPORT_SPLIT = ("has been split into --sjc_support, for spliced isoforms, "
                  "--single_exon_support, for single exon isoforms, and --frac_support, the share "
                  "of its gene's reads an isoform must have")
_PERCENT_USAGE = ("is now --frac_support, a fraction of its gene's reads where this was a "
                  "percentage of them: 5 becomes 0.05, the default of both")
_PARALOGS = ("has been removed: a read aligning equally well to several transcripts is assigned "
             "to one of them, the same one for every such read, rather than to all of them")

# names dropped in one subcommand or one program, by the name it was dropped from
RETIRED_BY_COMMAND = {
    'transcriptome': {
        '--genomealignedbam': "is now -b, --genome_aligned_bam",
        '--se_support': "is now --single_exon_support",
        '--keep_sup': "is now --keep_supplementary",
        '--noaligntoannot': "is now --no_align_to_annot",
        '--support': _SUPPORT_SPLIT,
        '-s': _SUPPORT_SPLIT,
        '-w': "is now --end_window",
        '--no_redundant': ("has been removed: the ends reported for a splice junction chain are set "
                           "by --max_ends and --end_window"),
        '-n': ("has been removed: the ends reported for a splice junction chain are set by "
               "--max_ends and --end_window"),
        '--stringent': _STRINGENT,
        '--check_splice': _CHECK_SPLICE,
        '--no_stringent': _NO_STRINGENT.format(other="--no_check_splice"),
        '--no_check_splice': _NO_STRINGENT.format(other="--no_stringent"),
        '--trust_ends': ("has been removed: each isoform's ends come from the ends of the reads "
                         "assigned to it"),
        '--allow_paralogs': _PARALOGS,
        '--parallel_mode': ("has been removed: the reads and annotated transcripts are partitioned "
                            "into clusters of overlapping features, and -t sets how many run at once"),
        '--parallelmode': ("has been removed: the reads and annotated transcripts are partitioned "
                           "into clusters of overlapping features, and -t sets how many run at once"),
        '--predictCDS': _PREDICT_CDS,
        '--longestORF': _PREDICT_CDS,
        '--annotation_reliant': _ALIGN_TO_ANNOT,
        '--annotated_bed': _ALIGN_TO_ANNOT,
        '--transcriptomefasta': _ALIGN_TO_ANNOT,
        '--remove_internal_priming': _INTERNAL_PRIMING,
        '--intprimingthreshold': _INTERNAL_PRIMING,
        '--intprimingfracAs': _INTERNAL_PRIMING,
        '--promoters': _ENDS_FROM_READS,
        '-p': _ENDS_FROM_READS,
        '--3prime_regions': _ENDS_FROM_READS,
        '--gene_tss': _ENDS_FROM_READS,
        '-i': _ENDS_FROM_READS,
        '--no_gtf_end_adjustment': _ENDS_FROM_READS,
        '--nvrna': ("has been removed: correcting a read's strand against the annotation is what "
                    "correction does, and --trust_strand and --directRNA keep the alignment's "
                    "strand instead"),
        '--mm2_args': "has been removed: the minimap2 runs are internal",
        '--longshot_bam': "has been removed: haplotypes come from flair isoalleles",
        '--longshot_vcf': "has been removed: haplotypes come from flair isoalleles",
    },
    'quantify': {
        '--reads_manifest': "is now --manifest",
        '-r': "is now --manifest",
        '-i': "is now --isoform_bed",
        '--isoformbed': "is now --isoform_bed",
        '--isoforms': ("has been removed: the isoforms are given as --isoform_bed with -g, "
                       "--genome, rather than as a FASTA"),
        '--sample_id_only': "has been removed: a counts column always names the sample alone",
        '--tpm': ("has been removed: <output>.isoforms.cpm.tsv is written by every run, and holds "
                  "the counts per million this option wrote as a tpm file"),
        '--with_gene': "has been removed: the counts matrix always names the gene, in its gene_id column",
        '--quality': "has been removed: it was accepted and never used",
        '--stringent': _STRINGENT,
        '--check_splice': _CHECK_SPLICE,
        '--output_bam': "has been removed: --generate_map reports the isoform each read was assigned to",
    },
    'combine': {
        '--endwindow': "is now --end_window",
        '-w': "is now --end_window",
        '--minpercentusage': _PERCENT_USAGE,
        '--min_percent_usage': _PERCENT_USAGE,
        '-p': _PERCENT_USAGE,
        '--remove_se': "is now --remove_single_exon",
        '--include_se': ("has been removed: single exon isoforms are combined unless "
                         "--remove_single_exon is given"),
        '-s': ("has been removed: single exon isoforms are combined unless --remove_single_exon "
               "is given"),
        '--convert_gtf': "has been removed: a GTF of the combined isoforms is always written",
        '-c': "has been removed: a GTF of the combined isoforms is always written",
        '--filter': "has been removed: --frac_support and --min_reads select the isoforms kept",
        '-m': "is now --manifest",
    },
    'fusion': {
        '--genomechimbam': "is now -b, --genome_aligned_bam",
        '--transcriptchimbam': ("has been removed: the reads are aligned to the annotated "
                                "transcripts by the run itself, which --no_transcriptome_alignment skips"),
        '--minfragmentsize': "is now --min_fragment_size",
        '--support': "is now --min_support",
        '-s': "is now --min_support",
        '--maxloci': "is now --max_loci",
        '--max_dist_to_TSS': "is now --max_dist_to_tss",
        '--min_dist_between_bp': "is now --min_dist_between_breakpoints",
        '--allow_paralogs': _PARALOGS,
        '--reads': "has been removed: fusion takes the aligned reads, -b, --genome_aligned_bam",
        '-r': "has been removed: fusion takes the aligned reads, -b, --genome_aligned_bam",
        '--stringent': _STRINGENT,
        '--check_splice': _CHECK_SPLICE,
    },
    'diffexp': {
        '--exp_thresh': "is now --min_expression",
        '-e': "is now --min_expression",
        '-q': "is now --counts_matrix",
    },
    'diffsplice': {
        '-q': "is now --counts_matrix",
        '-i': "is now --isoform_bed",
        '--drim1': "is now --min_samps_gene_expr",
        '--drim2': "is now --min_samps_feature_expr",
        '--drim3': "is now --min_gene_expr",
        '--drim4': "is now --min_feature_expr",
        '--conditionA': "is now --condition_a",
        '--conditionB': "is now --condition_b",
    },
    'spliceevents': {
        '--annot': "is now -f, --gtf",
        '--junc_support': "is now --junction_support",
        '--keep_sup': "is now --keep_supplementary",
        '--noaligntoannot': "is now --no_align_to_annot",
        '--check_outliers': ("has been removed: the spliceevents_outliers program calls outliers "
                             "from <output>.spliceevents.counts.tsv and totals.tsv, of either PSI"),
        '--output_read_ends': ("has been removed, with its <output>.spliceevents.readends.bed: the "
                               "corrected reads are passed to event calling as a line per junction "
                               "chain with its reads' ends"),
    },
    'isoalleles': {
        '--norm_bam': "is now --normal_bam",
        '--iso_read_map': "is now --read_map",
        '--iso_read_map_norm': "is now --normal_read_map",
        '--norm_vcf': ("has been removed: the normal's phased variants are --phased_vcf, with "
                       "--germline_vcf for the calls it does not phase"),
        '--allele_vcf': ("has been removed: the variants are read from --phased_vcf, --germline_vcf "
                         "and the somatic options"),
        '--allele_read_map': ("has been removed: allele groups are no longer called from the reads; "
                              "each isoform's reads are split by the haplotypes of --phased_vcf"),
        '--allele_read_map_norm': ("has been removed: allele groups are no longer called from the "
                                   "reads; each isoform's reads are split by the haplotypes of "
                                   "--phased_vcf"),
    },
    'bed_to_gtf': {
        '--noCDS': "is now --no_cds",
        '--dontuseCDS': "is now --no_cds",
        '--add_reference_transcript_id': ("has been removed: --is_flair_bed writes the extra "
                                          "columns of a flair BED as GTF attributes, the reference "
                                          "transcript id among them"),
    },
    'plot_isoform_usage': {
        '--vcf': "has been removed: the plot is of isoform usage alone",
    },
    'count_sam_transcripts': {
        '--quality': ("has been removed: reads are aligned to transcript sequences, where most align "
                      "about equally well to several isoforms and get a mapping quality of 0; "
                      "choosing between those isoforms is what this program does"),
        '--fusion_dist': "has been removed: it only turned on --trust_ends",
        '--end_norm_dist': ("has been removed: it relaxed the coverage check for single exon "
                            "transcripts, which end normalization never extends"),
    },
}

RETIRED_SUBCOMMANDS = {
    'align': ("has been removed: align the reads with minimap2, as the documentation on aligning "
              "reads shows, and give flair transcriptome the BAM with -b"),
    'correct': "has been removed: flair transcriptome builds transcript models from the aligned reads",
    'collapse': "has been removed: flair transcriptome builds transcript models from the aligned reads",
    'alleles': ("has been removed: flair isoalleles splits each isoform's reads by haplotype, from "
                "reads haplotagged with a phasing tool such as WhatsHap"),
    'variants': ("has been removed: flair isoalleles groups each isoform's reads by the variants "
                 "they contain"),
    'variantquant': ("has been removed: flair isoalleles groups each isoform's reads by the variants "
                     "they contain"),
}

def _option_names(argv):
    "the option names of a command line, with the value of --name=value dropped"
    args = itertools.takewhile(lambda arg: arg != '--', argv)
    return [arg.split('=', 1)[0] for arg in args if arg.startswith('-') and arg != '-']

def subcommand_parsers(parser):
    "the subparsers of a parser by subcommand name, or None if it has no subcommands"
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            return action.choices
    return None

def _report(parser, name, explanation):
    parser.error(f"{name} {explanation}")

def _check_options(parser, retired, argv):
    for name in _option_names(argv):
        if name in retired:
            _report(parser, name, retired[name])

def _check_subcommands(parser, argv):
    for name in argv:
        if name in RETIRED_SUBCOMMANDS:
            _report(parser, f"flair {name}", RETIRED_SUBCOMMANDS[name])

def check_retired(parser, argv):
    """Report a subcommand or option name an earlier release accepted, stating what
    replaced it.  A parser with subcommands is checked against the tables of the
    subcommand its command line names, one without them against the table of its
    program."""
    subcommands = subcommand_parsers(parser)
    if subcommands is None:
        _check_options(parser, RETIRED_BY_COMMAND.get(parser.prog, {}), argv)
    else:
        named = [arg for arg in argv if arg in subcommands]
        if len(named) == 0:
            _check_subcommands(parser, argv)
        else:
            _check_options(subcommands[named[0]],
                           RETIRED_EVERYWHERE | RETIRED_BY_COMMAND.get(named[0], {}), argv)

class RetiredCheckParser(argparse.ArgumentParser):
    """ArgumentParser that reports a name an earlier release accepted before parsing
    the command line.  The tables are keyed by subcommand, or by the prog of a
    parser without subcommands, so a program using this must set prog."""

    def parse_args(self, args=None, namespace=None):
        check_retired(self, sys.argv[1:] if args is None else list(args))
        return super().parse_args(args, namespace)

def check_tables(parser):
    """Report a subcommand with no table of its own, and any retired name that is
    also a live option, which would make the live option unusable.  Called by the
    tests, over the whole flair parser."""
    for subcommand, subparser in subcommand_parsers(parser).items():
        if subcommand not in RETIRED_BY_COMMAND:
            raise FlairError(f"flair {subcommand} has no table in cli_retired; add one, empty if "
                             f"the subcommand is new and has dropped no names")
        retired = RETIRED_EVERYWHERE | RETIRED_BY_COMMAND[subcommand]
        live = sorted(set(retired) & set(subparser._option_string_actions))
        if live:
            raise FlairError(f"flair {subcommand} has {', '.join(live)} as both a live and a "
                             f"retired option; remove it from cli_retired")
