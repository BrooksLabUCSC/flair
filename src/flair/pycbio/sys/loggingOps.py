"""
Operations associated with logging
"""
import logging
import os
import sys
from logging.handlers import SysLogHandler


def getFacilityNames():
    return tuple(SysLogHandler.facility_names.keys())


def getLevelNames():
    return tuple(sorted(logging._nameToLevel.keys()))


def parseFacility(facilityStr):
    """Convert case-insensitive facility string to a facility number."""
    facility = SysLogHandler.facility_names.get(facilityStr.lower())
    if facility is None:
        raise ValueError("invalid syslog facility: \"{}\"".format(facilityStr))
    return facility


def parseLevel(levelStr):
    "convert a log level string to numeric value"
    level = logging._nameToLevel.get(levelStr.upper())
    if level is None:
        raise ValueError("invalid logging level: \"{}\"".format(levelStr))
    return level


def _convertFacility(facility):
    """convert facility from string to number, if not already a number"""
    return facility if isinstance(facility, int) else parseFacility(facility)


def _convertLevel(level):
    """convert level from string to number, if not already a number"""
    return level if isinstance(level, int) else parseLevel(level)

def _loggerBySpec(logger):
    """Logger maybe a logger, logger name, or None for default logger, returns the
    logger."""
    if not isinstance(logger, logging.Logger):
        logger = logging.getLogger(logger)
    return logger

def setupLogger(logger, handler, level=logging.INFO, *, formatter=None):
    """add handle to logger and set logger level to the minimum of it's current
    and the handler level.  Logger maybe a logger, logger name, or None for
    default logger, returns the logger.

    """
    logger = _loggerBySpec(logger)
    if level is not None:
        logger.setLevel(_convertLevel(level))
    if handler.level is not None:
        logger.setLevel(min(handler.level, logger.level))
    logger.addHandler(handler)
    if formatter is not None:
        handler.setFormatter(formatter)
    return logger

def setupStreamLogger(logger, fh, level, *, formatter=None):
    """Configure logging to a specified open file.  Logger maybe a logger or
    logger name, returns the logger."""
    level = _convertLevel(level)
    handler = logging.StreamHandler(stream=fh)
    handler.setLevel(level)
    return setupLogger(logger, handler, formatter=formatter, level=level)


def setupStderrLogger(logger=None, level=logging.INFO, *, formatter=None):
    """configure logging to stderr.  Logger maybe a logger, logger name, or
    None for default logger, returns the logger."""
    return setupStreamLogger(logger, sys.stderr, _convertLevel(level), formatter=formatter)


def getSyslogAddress():
    """find the address to use for syslog"""
    for dev in ("/dev/log", "/var/run/syslog"):
        if os.path.exists(dev):
            return dev
    return ("localhost", 514)


def setupSyslogLogger(logger, facility, level, *, prog=None, address=None, formatter=None):
    """configure logging to syslog based on the specified facility.  If prog
    specified, each line is prefixed with the name.  Logger maybe a logger or
    logger name, returns the logger."""
    if address is None:
        address = getSyslogAddress()
    handler = SysLogHandler(address=address, facility=facility)
    # add a formatter that includes the program name as the syslog ident
    if prog is not None:
        handler.setFormatter(logging.Formatter(fmt="{} %(message)s".format(prog)))
    handler.setLevel(level)
    return setupLogger(logger, handler, formatter=formatter)


def setupNullLogger(logger, level=logging.INFO):
    "configure discard logging.  Returns logger."
    handler = logging.NullHandler()
    if level is not None:
        handler.setLevel(_convertLevel(level))
    return setupLogger(logger, handler)


def addCmdOptions(parser, *, defaultLevel=logging.INFO, inclSyslog=False):
    """
    Add command line options related to logging.  None of these are defaulted,
    as one might need to determine if they were explicitly set. The use case
    being getting a value from a configuration file if it is not specified on
    the command line.
    """
    # want to validate name, but want to store the string in the arguments vector
    # rather than the numeric value.
    def validateFacility(facilityStr):
        parseFacility(facilityStr)
        return facilityStr

    def validateLevel(levelStr):
        parseLevel(levelStr)
        return levelStr
    if inclSyslog:
        parser.add_argument("--syslog-facility", type=validateFacility,
                            help="Set syslog facility to case-insensitive symbolic value, if not specified, logging is not done to stderr, "
                            " one of {}".format(", ".join(getFacilityNames())))
    parser.add_argument("--log-stderr", action="store_true",
                        help="also log to stderr, even when logging to syslog")
    parser.add_argument("--log-level", type=validateLevel, default=defaultLevel,
                        help="Set level to case-insensitive symbolic value, one of {}".format(", ".join(getLevelNames())))
    parser.add_argument("--log-conf",
                        help="Python logging configuration file, see logging.config.fileConfig()")
    parser.add_argument("--log-debug", action="store_true",
                        help="short-cut that that sets --log-stderr and --log-level=DEBUG")

def setupFromCmd(args, *, logger=None, prog=None):
    """configure logging based on command options. Prog is used it to set the
    syslog program name. If prog is not specified, it is obtained from sys.arg.
    Logger maybe a logger, logger name, or None for default logger, returns the logger.

    Will log to stderr if not other login option is specified.

    N.B: logging must be initialized after daemonization
    """
    if args.log_debug:
        args.log_stderr = True
        args.log_level = logging.DEBUG
    if prog is None:
        prog = os.path.basename(sys.argv[0])
    logger = _loggerBySpec(logger)
    level = _convertLevel(args.log_level) if args.log_level is not None else logging.WARN
    haveFacility = ('syslog_facility' in args) and (args.syslog_facility is not None)
    if haveFacility:
        setupSyslogLogger(logger, args.syslog_facility, level, prog=prog)
    if (not haveFacility) or args.log_stderr:
        setupStderrLogger(logger, level)
    if args.log_conf is not None:
        logging.config.fileConfig(args.log_conf)
    return logger


class StreamToLogger:
    """
    File-like stream object that redirects writes to a logger instance.
    """
    # taken from
    # https://stackoverflow.com/questions/19425736/how-to-redirect-stdout-and-stderr-to-logger-in-python

    def __init__(self, logger, level):
        self.logger = logger
        self.level = level
        self.linebuf = ''

    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.level, line.rstrip())

    def flush(self):
        pass
