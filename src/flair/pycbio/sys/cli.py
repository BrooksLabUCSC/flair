# Copyright 2006-2025 Mark Diekhans
"""Support for command line parsing"""
import logging
import traceback
import io
import argparse
from flair.pycbio import PycbioException, NoStackError
from flair.pycbio.sys import loggingOps
from flair.pycbio.sys.objDict import ObjDict

###
# option parsing
###
def _find_subparser_action(parser):
    """Find subparser action in a parser"""
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            return action
    return None

def _findParserUsed(parser, parsed_args):
    """Find the parser object used to parse the arguments; required for subparsers.
    """
    # Traverse through the subparser hierarchy using parsed argument values
    while parser:
        subparsers_action = _find_subparser_action(parser)
        if subparsers_action is None:
            break

        # Find the next command in parsed args
        cmd_name = getattr(parsed_args, subparsers_action.dest, None)
        if (cmd_name is None) or (cmd_name not in subparsers_action.choices):
            break
        parser = subparsers_action.choices[cmd_name]
    return parser

def _getOptionDests(parser):
    return {a.dest for a in parser._actions if a.option_strings}

def _findOptNames(parser, parsed_args):
    "Find option names"
    # This is ugly, there isn't a good way to walk back the
    # chain of parsers used.
    usedParser = _findParserUsed(parser, parsed_args)
    return _getOptionDests(usedParser) | _getOptionDests(parser)

def splitOptionsArgs(parser, parsed_args):
    """Split command line arguments into two objects one of option arguments
    (-- or - options) and one of positional arguments.  Useful for packaging up
    a large number of options to pass around.  This does not handle subparsers
    unless the subparsers is the supplied parser, not the top-level one used
    in parsing"""

    opts = ObjDict()
    args = ObjDict()
    optnames = _findOptNames(parser, parsed_args)
    for name, value in vars(parsed_args).items():
        if name in optnames:
            opts[name] = value
        else:
            args[name] = value
    return opts, args

def parseOptsArgs(parser, *, args=None, namespace=None):
    """Call argparse parse_args and return (opts, args)"""
    return splitOptionsArgs(parser, parser.parse_args(args, namespace))

def parseArgsWithLogging(parser, *, args=None, namespace=None, defaultLevel=logging.WARNING,
                         inclSyslog=False):
    """Call argparse.parse_args and return args. Add logging command options
    if they are not already there and configuring logging after parsing.  This
    handles common cases.

    WARNING: this does not work with sub-parsers
    """
    if not loggingOps.haveCmdOptions(parser):
        loggingOps.addCmdOptions(parser, defaultLevel=defaultLevel, inclSyslog=inclSyslog)
    args = parser.parse_args(args, namespace)
    loggingOps.setupFromCmd(args)
    return args

def parseOptsArgsWithLogging(parser, args=None, namespace=None):
    """Call argparse.parse_args and return (opts, args). Add logging command options
    if they are not already there and configuring logging after parsing.  This
    handles common cases.

    WARNING: this does not work with sub-parsers
    """
    return splitOptionsArgs(parser, parseArgsWithLogging(parser, args=args, namespace=namespace))

###
# error handling
###
class CliError(PycbioException, NoStackError):
    "Error with command line arguments"
    pass

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
    fh = io.StringIO()
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
    DEFAULT_NO_STACK_EXCEPTS = (OSError, ImportError, KeyboardInterrupt)

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
