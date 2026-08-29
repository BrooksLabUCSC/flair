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
    "parser for comma-separated string list into a list"
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


# TSV typeMap tuple for str arrays
strArrayType = (strArraySplit, strArrayJoin)

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


# TSV typeMap tuple for str arrays
intArrayType = (intArraySplit, intArrayJoin)

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


# TSV typeMap tuple for str arrays
floatArrayType = (floatArraySplit, floatArrayJoin)
