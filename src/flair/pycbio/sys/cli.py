# Copyright 2006-2025 Mark Diekhans
"""Support for command line parsing"""
import argparse
import logging
import traceback
from io import StringIO
from flair.pycbio import NoStackError
from flair.pycbio.sys.objDict import ObjDict
from flair.pycbio.sys import loggingOps

def splitOptionsArgs(parser, inargs):
    """Split command line arguments into two objects one of option arguments
    (-- or - options) and one of positional arguments.  Useful for packaging up
    a large number of options to pass around."""

    opts = ObjDict()
    args = ObjDict()
    optnames = set([a.dest for a in parser._actions if a.option_strings])

    for name, value in vars(inargs).items():
        if name in optnames:
            opts[name] = value
        else:
            args[name] = value
    return opts, args

class ArgumentParserExtras(argparse.ArgumentParser):
    """Wrapper around ArgumentParser that adds logging
    related options.  Also can parse splitting options and positional
    arguments into separate objects.
    """
    def __init__(self, *args, incl_syslog=False, **kwargs):
        super().__init__(*args, **kwargs)
        self.incl_syslog = incl_syslog

    def add_extras(self):
        """Add extra options to parser, done just before parse
        so they are last in the help list.  Override to add other options.
        """
        loggingOps.addCmdOptions(self, inclSyslog=self.incl_syslog)

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
        the positional args as well. Returns (opts, args).
        """
        args = self.parse_args()
        return splitOptionsArgs(self, args)

def _exceptionPrintNoTraceback(exc, file, indent):
    depth = 0
    while exc:
        if depth > 0:
            prefix = depth * '  ' if indent else ''
            file.write(f"{prefix}Caused by: ")
        file.write(f"{type(exc).__name__}: {str(exc)}\n")
        exc = exc.__cause__ or exc.__context__
        depth += 1

def _exceptionPrint(exc, *, file=None, showTraceback=True, indent=True):
    """print a chained exception following causal chain"""
    if showTraceback:
        traceback.print_exception(exc, file=file)
    else:
        _exceptionPrintNoTraceback(exc, file, indent)

def _exceptionFormat(exc, *, showTraceback=True, indent=True):
    """format a chained exception following causal chain"""
    fh = StringIO()
    _exceptionPrint(exc, file=fh, showTraceback=showTraceback, indent=indent)
    return fh.getvalue()


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
        msg = _exceptionFormat(exc_val, showTraceback=showTraceback)
        if not showTraceback:
            msg += "Specify --log-debug for details"
        logger.error(msg.rstrip('\n'))

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is not None:
            if isinstance(exc_val, SystemExit):
                raise exc_val
            else:
                self._handleError(exc_val)
                exit(1)
