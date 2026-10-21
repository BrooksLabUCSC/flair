"""Reads grouped by the gene and junction chain they support, written per region
by flair spliceevents and read back when the region's samples are joined.

Columns, in order:

  gene         the gene the reads were assigned to
  juncs        the junction chain, each junction written start-end as genome
               coordinates are; comma separated in the file, and on the row a
               tuple of Junc, the splice junction PosRange
  start        the read's start
  end          the read's end
  strand       the strand of the chain
  read         the read name

Coordinates are zero based, half open, as everywhere in flair.

One row per read, so a chain appears once per supporting read.
"""
from flair import FlairInputDataError
from flair.isoform_data import Junc
from flair.pycbio.tsv import TsvReader, TsvWriter

COLUMNS = ('gene', 'juncs', 'start', 'end', 'strand', 'read')

JUNC_SEP = ','
COORD_SEP = '-'

def parse_junc_chain(value):
    "the chain as a tuple of Junc"
    if value == '':
        return ()
    return tuple(parse_junc(junc) for junc in value.split(JUNC_SEP))

def parse_junc(junc):
    try:
        start, end = junc.split(COORD_SEP)
        return Junc(int(start), int(end))
    except ValueError as ex:
        raise FlairInputDataError(
            f"junction '{junc}' is not start{COORD_SEP}end") from ex

def format_junc_chain(juncs):
    return JUNC_SEP.join(f'{junc.start}{COORD_SEP}{junc.end}' for junc in juncs)


juncChainType = (parse_junc_chain, format_junc_chain)

TYPE_MAP = {'gene': str, 'juncs': juncChainType, 'start': int, 'end': int,
            'strand': str, 'read': str}

class GeneJuncsReader(TsvReader):
    def __init__(self, gene_juncs_txt):
        super().__init__(gene_juncs_txt, typeMap=TYPE_MAP)

class GeneJuncsWriter(TsvWriter):
    def __init__(self, gene_juncs_txt, *, outFh=None):
        super().__init__(gene_juncs_txt, columns=COLUMNS, typeMap=TYPE_MAP,
                         outFh=outFh)
