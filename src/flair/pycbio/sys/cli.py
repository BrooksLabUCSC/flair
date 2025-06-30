# Copyright 2006-2025 Mark Diekhans
"""Support for command line parsing"""
import argparse
import logging
from flair.pycbio import NoStackError, exceptionFormat
from flair.pycbio.sys.objDict import ObjDict
from flair.pycbio.sys import loggingOps

def getOptionalArgs(parser, args):
    """Get the parse command line option arguments (-- or - options) as an
    object where the options are fields in the object.  Useful for packaging up
    a large number of options to pass around."""

    opts = ObjDict()
    for act in parser._actions:
        if (len(act.option_strings) > 0) and (act.dest != 'help'):
            setattr(opts, act.dest, getattr(args, act.dest))
    return opts

class ArgumentParserExtras(argparse.ArgumentParser):
    """Wrapper around ArgumentParser that adds logging
    related options.  Also can parse splitting options and arguments
    into separate objects.
    """
    def __init__(self, *args, add_profiler=False, **kwargs):
        super().__init__(*args, **kwargs)
        self.add_profiler = add_profiler

    def add_extras(self):
        """Add extra options to parser, done just before parse
        so they are last in the help list.  Override to add other options.
        """
        loggingOps.addCmdOptions(self)

    def process_extras(self, args):
        """Process extras after argments have been parsed"""
        loggingOps.setupFromCmd(args, prog=self.prog)

    def parse_known_args(self, *args, **kwargs):
        # parse args goes through parse_known_args
        self.add_extras()
        cmdargs, cmdargv = super().parse_known_args(*args, **kwargs)
        self.process_extras(cmdargs)
        return cmdargs, cmdargv

    def parse_opts_args(self):
        """Get the parse command line option arguments (-- or - arguments) as an
        object where the options are fields in the object.  Useful for packaging up
        a large number of options to pass around without the temptation to pass
        the positional args as well. Returns (opts, args).  The options are not
        removed from args.
        """
        args = self.parse_args()
        return getOptionalArgs(self, args), args

class ErrorHandler:
    """Command line error handling. This is a context manager that wraps execution
    of a program.  This is called after command line parsing has been done.
    If an exception is caught, it handles printing the error via logging, normally
    to stderr, and exits non-zero.

    If DEBUG log level is set, the call stack traceback will always be logged.
    Otherwise, the stack is not print if the Exception inherits from
    NoStackError, or its class is in the noStackExcepts, or printStackFilter
    indicates is should not be printed.

    noStackExcepts  - if not None, exceptions classes to not print, defaults to (OSError, ImportError).
    printStackFilter - a function that is called with the exception object. If it returns
    True, the stack is printed, False it is not printed, and None, it defers to the noStackExcepts
    checks.

    Example:
        def main():
           opts, args = cmd_line_parse()
           with cli.ErrorHandler():
              real_prog(opts, args.one, args.two)

    """
    DEFAULT_NO_STACK_EXCEPTS = (OSError, ImportError)

    def __init__(self, *, noStackExcepts=DEFAULT_NO_STACK_EXCEPTS, printStackFilter=None):
        self.noStackExcepts = tuple(noStackExcepts) if noStackExcepts is not None else None
        self.printStackFilter = printStackFilter

    def __enter__(self):
        pass

    def _showTraceBack(self, logger, exc_val):
        if logger.getEffectiveLevel() <= logging.DEBUG:
            return True
        if self.printStackFilter is not None:
            filt = self.printStackFilter(exc_val)
            if filt is not None:
                return filt
        if isinstance(exc_val, NoStackError):
            return False
        if self.noStackExcepts is None:
            return True
        return not isinstance(exc_val, self.noStackExcepts)

    def _handleError(self, exc_val):
        logger = logging.getLogger()
        showTraceback = self._showTraceBack(logger, exc_val)
        msg = exceptionFormat(exc_val, showTraceback=showTraceback)
        if not showTraceback:
            msg += "Specify --log-debug for details"
        logger.error(msg.rstrip('\n'))

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is not None:
            self._handleError(exc_val)
            exit(1)
