# Copyright 2006-2026 Mark Diekhans
"""support classes for parsing autoSql generated objects"""
from flair.pycbio import PycbioDataError


def _arrayElemStr(value, fmt=None):
    """a value as it appears in a comma-separated list; a missing value is an empty
    element, never the string "None", which would land in the data file"""
    if value is None:
        return ""
    return fmt.format(value) if fmt is not None else str(value)

##
# string array
##
def strArraySplit(commaStr):
    """parser for comma-separated string list into a list.  This follows kent
    sqlStringDynamicArray: a trailing comma terminates the list rather than adding an
    element, and an empty element is an empty string, so a lone ',' is a list of one
    empty string.  Use strArraySplitNone when an empty element means a missing value."""
    if len(commaStr) == 0:
        return []
    # autosql uses longblob, so if this came from a mysql database, we need to convert to bytes
    if isinstance(commaStr, bytes):
        commaStr = commaStr.decode('utf-8')
    strs = commaStr.split(",")
    if commaStr.endswith(","):
        strs = strs[0:-1]
    return strs


def strArrayJoin(strs):
    """formatter for a list of values into a comma separated string, not-str values are
    converted to a string"""
    if (strs is None) or (len(strs) == 0):
        return ""
    return ",".join([_arrayElemStr(s) for s in strs]) + ","


def strArraySplitNone(commaStr):
    """parser for a comma-separated string list where an empty element is a missing
    value rather than an empty string, so it parses to None.  This is the inverse of
    strArrayJoin, which writes a None as an empty element; plain strArraySplit follows
    autoSql and keeps the empty string."""
    return [s if s != '' else None for s in strArraySplit(commaStr)]


# TSV typeMap tuple for str arrays
strArrayType = (strArraySplit, strArrayJoin)

# TSV typeMap tuple for str arrays whose empty elements are missing values
strArrayNoneType = (strArraySplitNone, strArrayJoin)

##
# int arrays
##
def intArraySplit(commaStr):
    "parser for comma-separated string list into a list of ints"
    strs = strArraySplit(commaStr)
    try:
        return [int(s) for s in strs]
    except (TypeError, ValueError) as ex:
        raise PycbioDataError("not a comma-separated list of integers: '{}'".format(
            commaStr)) from ex


def intArrayJoin(ints):
    "formatter for a list of ints into a comma seperated string"
    if (ints is None) or (len(ints) == 0):
        return ""
    return ",".join([_arrayElemStr(i) for i in ints]) + ","


def intArraySplitNone(commaStr):
    """parser for comma-separated ints where an empty element is a missing value
    rather than an error, so it parses to None.  This is the inverse of intArrayJoin,
    which writes a None as an empty element."""
    strs = strArraySplitNone(commaStr)
    try:
        return [int(s) if s is not None else None for s in strs]
    except (TypeError, ValueError) as ex:
        raise PycbioDataError("not a comma-separated list of integers: '{}'".format(
            commaStr)) from ex


# TSV typeMap tuple for int arrays
intArrayType = (intArraySplit, intArrayJoin)

# TSV typeMap tuple for int arrays whose empty elements are missing values
intArrayNoneType = (intArraySplitNone, intArrayJoin)

##
# float arrays
##
def floatArraySplit(commaStr):
    "parser for comma-separated string list into a list of floats"
    strs = strArraySplit(commaStr)
    try:
        return [float(s) for s in strs]
    except (TypeError, ValueError) as ex:
        raise PycbioDataError("not a comma-separated list of numbers: '{}'".format(
            commaStr)) from ex


def floatArrayJoin(floats, fmt=None):
    "formatter for a list of floats a comma seperated string"
    if (floats is None) or (len(floats) == 0):
        return ""
    return ",".join([_arrayElemStr(f, fmt) for f in floats]) + ","


def floatArraySplitNone(commaStr):
    """parser for comma-separated floats where an empty element is a missing value
    rather than an error, so it parses to None.  This is the inverse of floatArrayJoin,
    which writes a None as an empty element."""
    strs = strArraySplitNone(commaStr)
    try:
        return [float(s) if s is not None else None for s in strs]
    except (TypeError, ValueError) as ex:
        raise PycbioDataError("not a comma-separated list of numbers: '{}'".format(
            commaStr)) from ex


# TSV typeMap tuple for float arrays
floatArrayType = (floatArraySplit, floatArrayJoin)

# TSV typeMap tuple for float arrays whose empty elements are missing values
floatArrayNoneType = (floatArraySplitNone, floatArrayJoin)
