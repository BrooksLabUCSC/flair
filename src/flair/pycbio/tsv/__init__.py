# Copyright 2006-2025 Mark Diekhans
"""" TSV (Tab Separated File) parsing"""

from flair.pycbio import PycbioException


class TsvError(PycbioException):
    "Error from reading or parsing a TSV file"
    def __init__(self, msg, reader=None):
        if (reader is not None):
            msg = str(reader.fileName) + ":" + str(reader.lineNum) + ": " + msg
        super(TsvError, self).__init__(msg)


from flair.pycbio.tsv.tsvRow import TsvRow, tsvRowToDict
from flair.pycbio.tsv.tsvReader import TsvReader, strOrNoneType, intOrNoneType, printf_basic_dialect
from flair.pycbio.tsv.tabFile import TabFile
from flair.pycbio.tsv.tabFile import TabFileReader

__all__ = (TsvError.__name__, TsvRow.__name__, TsvReader.__name__,
           "strOrNoneType", "intOrNoneType", "floatOrNoneType",
           TabFile.__name__, TabFileReader.__name__,
           tsvRowToDict.__name__,
           printf_basic_dialect.__name__)
