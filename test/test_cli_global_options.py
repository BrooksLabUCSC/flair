"""The global options of flair, which are accepted only before the subcommand.  Given
after it, the error states the command written again, which must itself be a command
that runs."""
import pytest
from flair.flair_cli import build_parser


def _error_message(capsys, argv):
    with pytest.raises(SystemExit):
        build_parser().parse_args(argv)
    return capsys.readouterr().err


def test_the_correction_of_an_option_with_a_value_keeps_the_value(capsys):
    # without the value, the correction reads the subcommand name as the level
    message = _error_message(capsys, ['transcriptome', '--log-level=DEBUG'])
    assert 'flair --log-level DEBUG transcriptome ...' in message


def test_the_value_of_a_separate_argument_is_kept_too(capsys):
    message = _error_message(capsys, ['transcriptome', '--log-level', 'DEBUG'])
    assert 'flair --log-level DEBUG transcriptome ...' in message


def test_the_correction_of_a_flag_leaves_the_next_argument_alone(capsys):
    message = _error_message(capsys, ['quantify', '-z', 'reads.bam'])
    assert 'flair -z quantify ...' in message


def test_a_global_option_before_the_subcommand_is_accepted():
    args = build_parser().parse_args(['-z', 'transcriptome', '-b', 'x.bam', '-g', 'g.fa',
                                      '--sample_name', 's'])
    assert args.compress
    assert args.subcommand == 'transcriptome'
