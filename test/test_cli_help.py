import argparse
from flair.cli_help import FlairHelpFormatter, help_blocks

HELP = ('which ones to keep:\n'
        '- nosubset: none of them\n'
        '- ginormous: all of them\n'
        '(default: nosubset)')


def test_help_blocks_groups_items():
    assert help_blocks(HELP) == [('text', ['which ones to keep:']),
                                 ('items', ['- nosubset: none of them',
                                            '- ginormous: all of them']),
                                 ('text', ['(default: nosubset)'])]


def test_help_blocks_joins_prose_lines():
    assert help_blocks('a help string\nbroken over lines') == [('text', ['a help string broken over lines'])]


def _option_help(help_text, width=40):
    parser = argparse.ArgumentParser(prog='t', formatter_class=FlairHelpFormatter)
    parser.add_argument('--filter', help=help_text)
    return parser._get_formatter()._split_lines(help_text, width)


def test_items_are_their_own_lines():
    assert _option_help(HELP) == ['which ones to keep:',
                                  '- nosubset: none of them',
                                  '- ginormous: all of them',
                                  '(default: nosubset)']


def test_item_continuation_is_indented():
    lines = _option_help('- ginormous: every spliced subset isoform and every single exon one', width=30)
    assert lines[0] == '- ginormous: every spliced'
    assert all(line.startswith('  ') for line in lines[1:])
