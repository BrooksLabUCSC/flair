# Copyright 2006-2025 Mark Diekhans

# FIXME: needed for faster readings, but needs cleaned up, need reader/writer
# classes
# FIXME add try and error msg with file/line num, move to row reader class; see fileOps.iterRows
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
    def processLine(line):
        if hashAreComments and line.startswith("#"):
            return None
        if skipBlankLines and (len(line) <= 1):  # include \n
            return None
        row = line[0:-1].split('\t')
        return rowClass(row) if rowClass is not None else row

    lineNum = -1
    with fileOps.FileAccessor(fspec) as fh:
        for line in fh:
            lineNum += 1
            row = processLine(line)
            if row is not None:
                yield row
