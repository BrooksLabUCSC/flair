"""Column converters and header handling shared by the flair TSV formats.

A converter is the (parse, format) pair a typeMap entry takes, so a column whose
text contains several values is a list on the row and is joined again on write.
Callers then never split or join a column themselves, and the separator is
stated once, here and in the format that uses it.
"""
from flair import FlairInputDataError
from flair.pycbio.sys import fileOps


def _joined(values, sep):
    """The values of a list column, joined.  A plain string is refused rather than
    joined character by character, which wrote GENEA as G,E,N,E,A."""
    if isinstance(values, str):
        raise FlairInputDataError(f"a list column was given the string {values!r}; it takes a list or "
                                  f"tuple, whose values are joined with {sep!r}")
    return sep.join(str(value) for value in values)

def separated_list_type(sep):
    "a column of strings joined by sep; the row contains the list of them"
    return (lambda value: [] if value == '' else value.split(sep),
            lambda values: _joined(values, sep))

def separated_int_list_type(sep):
    "the same for a column whose values are numbers"
    return (lambda value: [] if value == '' else [int(v) for v in value.split(sep)],
            lambda values: _joined(values, sep))


def _format_count(value):
    "a count is summed as a float, so 3.0 is written 3, the integer that is read back"
    return str(round(value))


# a read count: parsed as an integer and written as one
countType = (int, _format_count)


# comma separated, the usual case: read names, isoform ids, filter tags
commaListType = separated_list_type(',')

# '; ' separated, used where the values themselves contain commas
semicolonSpaceListType = separated_list_type('; ')

# ';' separated counts, used by the columns that pack a pair of counts into one
# field
semicolonIntListType = separated_int_list_type(';')

# comma separated start-end ranges, as junctions or read ends; the row contains a tuple
# of (start, end) pairs
def _parse_ranges(value):
    return () if value == '' else tuple(tuple(int(x) for x in pair.split('-')) for pair in value.split(','))


rangeListType = (_parse_ranges, lambda ranges: ','.join(f'{start}-{end}' for start, end in ranges))


def read_header_columns(tsv):
    """The column names of a TSV, for a format whose columns are not known until
    the file is open, such as one with a column per sample."""
    for line in fileOps.iterLines(tsv):
        # the leading # of a UCSC style header is not part of the name, as
        # TsvReader._dehashHeader also takes it
        return line.lstrip('#').split('\t')
    raise FlairInputDataError(f"{tsv}: file is empty, expected a header line")
