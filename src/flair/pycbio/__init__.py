# Copyright 2006-2025 Mark Diekhans
import traceback
from io import StringIO

class NoStackError:
    """When added as an additional base class to an Exception, stack traces
    should not be printed to users expect for debug logging.  These are
    normally caused by user input where print the stack is not helpful.
    """
    __slots__ = ()

class PycbioException(Exception):
    """Base class for Pycbio exceptions."""
    pass

class PycbioOptionalFeatureException(PycbioException, NoStackError):
    """Thrown to indicate an optional feature is missing"""
    def __init__(self, feature, reason):
        super(PycbioOptionalFeatureException, self).__init__(f"{feature} is disabled due to {reason}")
        self.feature = feature
        self.reason = reason

def _exceptionPrintNoTraceback(exc, file, indent):
    depth = 0
    while exc:
        if depth > 0:
            prefix = depth * '  ' if indent else ''
            file.write(f"{prefix}Caused by: ")
        file.write(f"{type(exc).__name__}: {exc}\n")
        exc = exc.__cause__ or exc.__context__
        depth += 1

def exceptionPrint(exc, *, file=None, showTraceback=True, indent=True):
    """print a chained exception following causal chain"""
    if showTraceback:
        traceback.print_exception(exc, file=file)
    else:
        _exceptionPrintNoTraceback(exc, file, indent)

def exceptionFormat(exc, *, showTraceback=True, indent=True):
    """format a chained exception following causal chain"""
    fh = StringIO()
    exceptionPrint(exc, file=fh, showTraceback=showTraceback, indent=indent)
    return fh.getvalue()
