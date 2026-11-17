"""Help text that can list an option's values.

A help string's lines are kept rather than run together, and a line starting with
`- ' is an item of a list.  FlairHelpFormatter wraps them that way in the terminal,
and dev/bin/flair-gen-cli-doc writes them as a bullet list on the command line pages,
both through help_blocks().
"""
import argparse
import re
import textwrap

_WHITESPACE_RE = re.compile(r'\s+')
_ITEM_PREFIX = '- '
_ITEM_INDENT = '  '


def _help_lines(text):
    "the help string's lines, each one's whitespace collapsed, the empty ones dropped"
    lines = (_WHITESPACE_RE.sub(' ', line).strip() for line in text.splitlines())
    return [line for line in lines if line != '']


def _block_kind(line):
    return 'items' if line.startswith(_ITEM_PREFIX) else 'text'


def help_blocks(text):
    """The help string as blocks of consecutive lines of one kind: ('items', lines)
    for a run of `- ' items, each its own line, and ('text', [paragraph]) for the
    prose around them, whose own line breaks are not the writer's"""
    blocks = []
    for line in _help_lines(text):
        if blocks and blocks[-1][0] == _block_kind(line):
            blocks[-1][1].append(line)
        else:
            blocks.append((_block_kind(line), [line]))
    return [(kind, lines if kind == 'items' else [' '.join(lines)])
            for kind, lines in blocks]


class FlairHelpFormatter(argparse.HelpFormatter):
    """Help that keeps the lines of a help string, so an option can list its values,
    an item's continuation indented under it"""

    def _split_lines(self, text, width):
        lines = []
        for _, block in help_blocks(text):
            lines.extend(self._wrap_block(block, width))
        return lines

    def _wrap_block(self, block, width):
        wrapped = []
        for line in block:
            indent = _ITEM_INDENT if line.startswith(_ITEM_PREFIX) else ''
            wrapped.extend(textwrap.wrap(line, width, subsequent_indent=indent))
        return wrapped
