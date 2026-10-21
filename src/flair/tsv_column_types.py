"""Column converters and header handling shared by the flair TSV formats.

A converter is the (parse, format) pair a typeMap entry takes, so a column whose
text holds several values is a list on the row and is joined again on write.
Callers then never split or join a column themselves, and the separator is
stated once, here and in the format that uses it.
"""
from flair import FlairInputDataError
from flair.pycbio.sys import fileOps


def separated_list_type(sep):
    "a column holding strings joined by sep; the row carries the list of them"
    return (lambda value: [] if value == '' else value.split(sep),
            lambda values: sep.join(values))

def separated_int_list_type(sep):
    "the same for a column whose values are numbers"
    return (lambda value: [] if value == '' else [int(v) for v in value.split(sep)],
            lambda values: sep.join(str(v) for v in values))


# comma separated, the usual case: read names, isoform ids, filter tags
commaListType = separated_list_type(',')

# '; ' separated, used where the values themselves contain commas
semicolonSpaceListType = separated_list_type('; ')

# ';' separated counts, used by the columns that pack a pair of counts into one
# field
semicolonIntListType = separated_int_list_type(';')

# comma separated counts
commaIntListType = separated_int_list_type(',')


def read_header_columns(tsv):
    """The column names of a TSV, for a format whose columns are not known until
    the file is open, such as one with a column per sample."""
    for line in fileOps.iterLines(tsv):
        return line.split('\t')
    raise FlairInputDataError(f"{tsv}: file is empty, expected a header line")
