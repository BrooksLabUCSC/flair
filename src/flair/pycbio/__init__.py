# Copyright 2006-2025 Mark Diekhans

class NoStackError:
    """When added as an additional base class to an Exception, stack traces
    should not be printed to users expect for debug logging.  These are
    normally caused by user input errors where print the stack is not helpful.
    """
    __slots__ = ()

class PycbioException(Exception):
    """Base class for Pycbio exceptions."""
    pass

class PycbioDataError(PycbioException, NoStackError):
    """Indicates an data format error or related error.  The pycbio.sys.cli.ErrorHandler
    facility will not generate stacktrace by default."""
    pass

class PycbioOptionalFeatureException(PycbioException, NoStackError):
    """Thrown to indicate an optional feature is missing"""
    def __init__(self, feature, reason):
        super(PycbioOptionalFeatureException, self).__init__(f"{feature} is disabled due to {reason}")
        self.feature = feature
        self.reason = reason
