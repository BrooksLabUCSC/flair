"""The subcommand and option names earlier releases accepted.  None of them is an
option of any parser: each is reported before the command line is parsed, with the
name that replaced it."""
import pytest
from flair import FlairInputDataError
from flair.flair_cli import build_parser
from flair.flair_transcriptome import parse_filter
from flair.bed_to_gtf import build_parser as build_bed_to_gtf_parser
from flair.plot_isoform_usage import build_parser as build_plot_isoform_usage_parser
from flair.cli_retired import (RETIRED_EVERYWHERE, RETIRED_BY_COMMAND, RETIRED_SUBCOMMANDS,
                               check_tables, subcommand_parsers)

PROGRAM_PARSERS = {'bed_to_gtf': build_bed_to_gtf_parser,
                   'plot_isoform_usage': build_plot_isoform_usage_parser}


def _error_message(capsys, argv):
    with pytest.raises(SystemExit):
        build_parser().parse_args(argv)
    return capsys.readouterr().err


def _retired_options():
    "every retired option of every subcommand, with the subcommand it is retired from"
    return [(subcommand, name)
            for subcommand in subcommand_parsers(build_parser())
            for name in sorted(RETIRED_EVERYWHERE | RETIRED_BY_COMMAND[subcommand])]


@pytest.mark.parametrize("subcommand,name", _retired_options(), ids=lambda value: value.lstrip('-'))
def test_every_retired_option_reports_what_replaced_it(capsys, subcommand, name):
    message = _error_message(capsys, [subcommand, name])
    assert f"flair {subcommand}: error: {name} " in message


@pytest.mark.parametrize("subcommand", sorted(RETIRED_SUBCOMMANDS))
def test_every_retired_subcommand_reports_what_replaced_it(capsys, subcommand):
    message = _error_message(capsys, [subcommand, '-o', 'out'])
    assert f"flair: error: flair {subcommand} " in message


def test_no_retired_name_is_also_a_live_option():
    check_tables(build_parser())


def test_a_value_given_with_an_equals_sign_is_part_of_the_name(capsys):
    message = _error_message(capsys, ['combine', '--minpercentusage=10'])
    assert 'flair combine: error: --minpercentusage is now --min_percent_usage' in message


def test_a_retired_name_is_reported_although_required_arguments_are_missing(capsys):
    # the point of checking before parsing: a command line being converted has
    # neither the new names nor, yet, the arguments they are required for
    message = _error_message(capsys, ['transcriptome', '--se_support', '3'])
    assert 'flair transcriptome: error: --se_support is now --single_exon_support' in message


def test_a_name_retired_from_one_subcommand_is_live_in_another():
    # --trust_ends is removed from transcriptome and is an option of quantify
    args = build_parser().parse_args(['quantify', '--manifest', 'm.tsv', '-g', 'g.fa',
                                      '--isoform_bed', 'i.bed', '--trust_ends'])
    assert args.trust_ends


def test_an_abbreviation_of_a_live_option_still_resolves():
    # the retired names are not options, so --out_dir and --output_prefix do not
    # make --out ambiguous
    args = build_parser().parse_args(['transcriptome', '-b', 'x.bam', '-g', 'g.fa',
                                      '--sample_name', 'sample', '--out', 'out'])
    assert args.output == 'out'


def _program_options():
    "every retired option of the programs with a parser of their own"
    return [(program, name)
            for program in PROGRAM_PARSERS
            for name in sorted(RETIRED_BY_COMMAND[program])]


@pytest.mark.parametrize("program,name", _program_options(), ids=lambda value: value.lstrip('-'))
def test_every_retired_option_of_a_program_reports_what_replaced_it(capsys, program, name):
    with pytest.raises(SystemExit):
        PROGRAM_PARSERS[program]().parse_args([name])
    assert f"{program}: error: {name} " in capsys.readouterr().err


def test_the_retired_value_of_filter_reports_what_replaced_it():
    # a value rather than a name, so it is reported where the value is checked
    with pytest.raises(FlairInputDataError) as raised:
        parse_filter('bysupport')
    assert '--filter bysupport is now --filter <N>X' in str(raised.value)
