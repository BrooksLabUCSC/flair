# Copyright 2006-2025 Mark Diekhans
"""TSV writing classes"""
import csv
from flair.pycbio.sys import fileOps
from flair.pycbio.tsv import TsvError
from flair.pycbio.tsv.tsvRow import TsvRow
from flair.pycbio.tsv.tsvColumns import ColumnSpecs, columnsSpecBuild


class TsvWriter:
    """Class for writing TSV files with column definitions and optional
    type formatting.  Supports writing from dicts, lists, tuples,
    namedtuples, TsvRow objects, or any object with attributes matching
    column names.

    If the input is a list or tuple, values are assumed to be in column order.
    For dicts and objects with attributes, values are looked up by column name.

    Supports writing to compressed files (.gz, .bz2) via opengz.
    Can be used as a context manager.
    """

    def __init__(self, fileName, *, columns=None, typeMap=None, defaultColType=None,
                 columnSpecs=None, outFh=None, dialect=csv.excel_tab, encoding=None,
                 errors=None):
        """
        :param fileName: Path to the output TSV file, unless `outFh` is provided.
        :param columns: List of column names.  If omitted but `typeMap` is given,
            columns are taken from `typeMap` keys in insertion order.
        :param typeMap: Dictionary mapping column names to either a type or a
            (parseFunc, formatFunc) tuple.  Only the formatFunc is used for writing.
        :param defaultColType: Type to use for columns not listed in `typeMap`.
        :param columnSpecs: A pre-built `ColumnSpecs` to use directly.  Mutually
            exclusive with `columns`, `typeMap`, and `defaultColType`.  A
            `ColumnSpecs` can be obtained from a `TsvReader` (via
            `reader.columnSpecs`) or a `TsvRow` (via `tsvRowGetColumnSpecs(row)`),
            making it easy to write a TSV with the same schema as one being read.
        :param outFh: An open file-like object to write to instead of `fileName`.
            It will not be closed automatically.
        :param dialect: A `csv.Dialect` instance or dialect name.
        :param encoding: Optional text encoding (e.g., 'utf-8').
        :param errors: Optional error handling strategy.
        """
        self.fileName = fileName
        self.columnSpecs = self._resolveColumnSpecs(columns, typeMap, defaultColType, columnSpecs)
        self.outFh = None
        self._shouldClose = False
        self._open(fileName, outFh, dialect, encoding, errors)
        self._writeHeader()

    @staticmethod
    def _resolveColumnSpecs(columns, typeMap, defaultColType, columnSpecs):
        if columnSpecs is not None:
            if columns is not None or typeMap is not None or defaultColType is not None:
                raise TsvError("columnSpecs cannot be combined with columns, "
                               "typeMap, or defaultColType")
            if not isinstance(columnSpecs, ColumnSpecs):
                raise TsvError("columnSpecs must be a ColumnSpecs instance")
            return columnSpecs
        if columns is None:
            if typeMap is None:
                raise TsvError("must specify columns, typeMap, or columnSpecs")
            columns = list(typeMap.keys())
        return columnsSpecBuild(columns, typeMap, defaultColType)

    @property
    def columns(self):
        "column names"
        return self.columnSpecs.columns

    def _open(self, fileName, outFh, dialect, encoding, errors):
        if outFh is not None:
            self.outFh = outFh
            self._shouldClose = False
        else:
            self.outFh = fileOps.opengz(fileName, "w", encoding=encoding, errors=errors)
            self._shouldClose = True
        self._writer = csv.writer(self.outFh, dialect=dialect, lineterminator='\n')

    def close(self):
        """Close the file if this object opened it."""
        if self._shouldClose and self.outFh is not None:
            self.outFh.close()
            self.outFh = None

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        self.close()

    def _writeHeader(self):
        """Write the column header line."""
        self._writer.writerow(self.columns)

    def _rowFromSequence(self, row):
        """Format a list/tuple, assuming values are in column order."""
        return [self.columnSpecs.fmtValue(i, row[i]) for i in range(len(self.columns))]

    def _rowFromMapping(self, row):
        """Format a dict-like object, looking up values by column name."""
        return [self.columnSpecs.fmtValue(i, row[col]) for i, col in enumerate(self.columns)]

    def _rowFromObject(self, row):
        """Format an object with attributes matching column names."""
        return [self.columnSpecs.fmtValue(i, getattr(row, col)) for i, col in enumerate(self.columns)]

    def _formatRow(self, row):
        """Convert a row to a list of formatted strings."""
        if isinstance(row, dict):
            return self._rowFromMapping(row)
        elif isinstance(row, TsvRow):
            return self._rowFromObject(row)
        elif isinstance(row, tuple) and hasattr(row, '_fields'):
            # namedtuple — use attribute names, not positions
            return self._rowFromObject(row)
        elif isinstance(row, (list, tuple)):
            return self._rowFromSequence(row)
        else:
            # object with attributes
            return self._rowFromObject(row)

    def writeRow(self, row):
        """Write a single row.  Row can be a list, tuple, dict, TsvRow,
        or any object with attributes matching column names."""
        self._writer.writerow(self._formatRow(row))

    def writeColumns(self, **kwargs):
        """Write a row using keyword arguments for column values.
        Columns not present in kwargs are written as empty strings (the type
        formatter is not invoked).  An explicitly-passed None is still routed
        through the formatter."""
        row = []
        for i, col in enumerate(self.columns):
            if col in kwargs:
                row.append(self.columnSpecs.fmtValue(i, kwargs[col]))
            else:
                row.append("")
        self._writer.writerow(row)

    def writeRows(self, rows):
        """Write multiple rows."""
        for row in rows:
            self.writeRow(row)
