# Copyright 2006-2025 Mark Diekhans
"""TSV column and type information shared by TsvReader and TsvWriter."""
from flair.pycbio.tsv import TsvError


class ColumnSpecs:
    """Column names, types, and formatting information for TSV data.
    Disconnected from reader/writer to allow rows to be pickled."""
    __slots__ = ("columns", "extColumns", "columnMap", "types")

    def __init__(self, columns, extColumns, columnMap, types):
        self.columns = columns
        self.extColumns = extColumns
        self.columnMap = columnMap
        self.types = types

    def fmtValue(self, iCol, value):
        """Format a single value using the column type if available."""
        if self.types is not None:
            ct = self.types[iCol]
            if ct is not None:
                if isinstance(ct, tuple):
                    return ct[1](value)
                elif value is None:
                    return ""
                else:
                    return str(value)
        if value is None:
            return ""
        return str(value)

    def formatRow(self, row):
        "format row to a list of strings given specified types"
        return [self.fmtValue(iCol, row[iCol]) for iCol in range(len(self.columns))]


def _getColTypes(columns, typeMap, defaultColType):
    "save col types as column indexed list"
    if typeMap is not None:
        return [typeMap.get(col, defaultColType) for col in columns]
    elif defaultColType is not None:
        return len(columns) * [defaultColType]
    else:
        return None

def _setupColumns(columnNames, columnNameMapper):
    # n.b. columns could be passed in from client, must copy
    i = 0
    columns = []
    extColumns = [] if columnNameMapper is not None else columns
    colMap = {}
    for col in columnNames:
        if columnNameMapper is not None:
            extColumns.append(col)
            col = columnNameMapper(col)
        columns.append(col)
        if col in colMap:
            raise TsvError("Duplicate column name: '{}'".format(col))
        colMap[col] = i
        i += 1
    return columns, extColumns, colMap

def columnsSpecBuild(columnNames, typeMap, defaultColType, columnNameMapper=None):
    """Build a ColumnSpecs from column names, optional typeMap, and optional name mapper."""
    columns, extColumns, columnMap = _setupColumns(columnNames, columnNameMapper)
    colTypes = _getColTypes(columns, typeMap, defaultColType)
    return ColumnSpecs(columns, extColumns, columnMap, colTypes)
