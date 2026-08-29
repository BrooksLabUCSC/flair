# Copyright 2006-2026 Mark Diekhans

# FIXME: needed for faster readings, but needs cleaned up, need reader/writer
# classes
from flair.pycbio import PycbioDataError
from flair.pycbio.sys import fileOps


class TabFile(list):
    """Class for reading and hold tab-separated files.
    """

    def __init__(self, fileName, rowClass=None, hashAreComments=False, skipBlankLines=False):
        """Read tab file into the object
        """
        self.fileName = fileName
        self.rowClass = rowClass
        for row in TabFileReader(self.fileName, rowClass=rowClass, hashAreComments=hashAreComments, skipBlankLines=skipBlankLines):
            self.append(row)


def TabFileReader(fspec, rowClass=None, hashAreComments=False, skipBlankLines=False):
    """generator over tab file rows"""
    def buildRow(row, lineNum):
        try:
            return rowClass(row)
        except Exception as ex:
            raise PycbioDataError("{}:{}: can not build row from: {}".format(
                fileOps.fileSpecName(fspec), lineNum, row)) from ex

    def processLine(line, lineNum):
        # rstrip rather than line[:-1]: the last line of a file need not have a
        # newline, and dropping its final character loses a character from the
        # last column
        line = line.rstrip("\n")
        if hashAreComments and line.startswith("#"):
            return None
        if skipBlankLines and (line == ""):
            return None
        row = line.split('\t')
        return buildRow(row, lineNum) if rowClass is not None else row

    lineNum = 0
    with fileOps.FileAccessor(fspec) as fh:
        for line in fh:
            lineNum += 1
            row = processLine(line, lineNum)
            if row is not None:
                yield row
