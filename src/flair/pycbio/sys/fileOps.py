"""Miscellaneous file operations"""
# Copyright 2006-2026 Mark Diekhans

import os
import os.path as osp
from pathlib import Path
import sys
import re
import glob
import socket
import threading
import tempfile
import pipettor
from shutil import which
from contextlib import contextmanager
from flair.pycbio import PycbioException

# FIXME: normalize file line read routines to all take fh or name, remove redundant code.

def isFilePath(fspec):
    return isinstance(fspec, str) or isinstance(fspec, Path)

def ensureDir(dir):
    """Ensure that a directory exists, creating it (and parents) if needed, avoiding
    race conditions"""
    os.makedirs(dir, exist_ok=True)


def ensureFileDir(fname):
    """Ensure that the directory for a file exists, creating it (and parents) if needed.
    Returns the directory path"""
    dir = osp.dirname(fname)
    if len(dir) > 0:
        ensureDir(dir)
        return dir
    else:
        return "."

def unlinkIfExists(path, *, dir_fd=None):
    """unlink a file if it exists, ignoring FileNotFoundError, which prevents
    race conditions."""
    try:
        os.unlink(path, dir_fd=dir_fd)
    except FileNotFoundError:
        pass

def rmdirIfExists(path, *, dir_fd=None):
    """remove a directory if it exists, ignoring FileNotFoundError, which prevents
    race conditions."""
    try:
        os.rmdir(path, dir_fd=dir_fd)
    except FileNotFoundError:
        pass

def rmFiles(*files):
    """Remove one or more files if they exist File paths of None are skipped.
    Missing files don't generate an error"""
    for f in files:
        if f is not None:
            unlinkIfExists(f)


def _rmTreeEntries(dir, subdirs, files):
    """unlink the files and the symlinked directories of one directory of a walk;
    os.walk reports a symlink to a directory in subdirs and does not descend into
    it, so it must be unlinked here or its parent is not empty."""
    dir_fd = os.open(dir, os.O_DIRECTORY)
    try:
        for f in files:
            unlinkIfExists(f, dir_fd=dir_fd)
        for d in subdirs:
            if os.path.islink(osp.join(dir, d)):
                unlinkIfExists(d, dir_fd=dir_fd)
    finally:
        os.close(dir_fd)

def rmTree(root):
    """remove a file hierarchy, root can be a file, a symlink, or a directory,
    missing files don't generate errors.  Symlinks are removed, not followed."""
    if osp.isdir(root) and not osp.islink(root):
        for dir, subdirs, files in os.walk(root, topdown=False):
            _rmTreeEntries(dir, subdirs, files)
            rmdirIfExists(dir)
    else:
        unlinkIfExists(root)


FAST_REMOVE_SUFFIX = ".drop"    # marks a tree renamed aside, awaiting removal
_GLOB_CHARS = "*?["             # a path spec containing any of these is a pattern


def _glob_matches(spec):
    "the existing paths a spec refers to; a glob pattern may match several"
    spec = os.fspath(spec)
    if any(c in spec for c in _GLOB_CHARS):
        return sorted(glob.glob(spec))
    return [spec] if osp.lexists(spec) else []


def _partition_dirs(targets):
    "(directories, other paths) among targets; a symlink counts as an other path"
    dirs, others = [], []
    for target in targets:
        isdir = osp.isdir(target) and not osp.islink(target)
        (dirs if isdir else others).append(target)
    return dirs, others


def _rename_aside(path):
    "rename a directory to a unique <name>.<uniq>.drop sibling, returning the new name"
    parent = osp.dirname(osp.abspath(path))
    base = osp.basename(osp.normpath(path))
    drop = tempfile.mktemp(prefix=base + ".", suffix=FAST_REMOVE_SUFFIX, dir=parent)
    os.rename(path, drop)
    return drop


def _background_rm(paths):
    """fire off a detached `rm -rf` of paths; the caller neither waits nor checks it.
    setsid puts it in its own session, so it is not killed by a hangup on the
    terminal that started it."""
    pipettor.Pipeline(["setsid", "rm", "-rf", *paths], stdout=1, stderr=2).start()


def _report_removed(label, dropped, unlinked):
    "report what was renamed aside for background removal and what was unlinked"
    for path in dropped:
        prfErr(f"{label}: renamed aside, removing in background: {path}")
    for path in unlinked:
        prfErr(f"{label}: removed {path}")
    if not (dropped or unlinked):
        prfErr(f"{label}: nothing to remove")


def fast_remove(specs, label="clean"):
    """Discard paths without waiting for the (slow) recursive delete: each existing
    directory is renamed aside to a unique .drop sibling and one detached `rm -rf`
    is fired off for them all, so the caller returns immediately and the space is
    reclaimed on its own; other paths (plain files, symlinks) are unlinked outright.
    Specs may be str or Path, and may be glob patterns.  Renaming is a rename within
    the parent directory, so a spec's parent must be on one filesystem.  Returns the
    renamed-aside names.
    """
    targets = [match for spec in specs for match in _glob_matches(spec)]
    dirs, others = _partition_dirs(targets)
    dropped = [_rename_aside(path) for path in dirs]
    for path in others:
        os.unlink(path)
    _report_removed(label, dropped, others)
    if dropped:
        _background_rm(dropped)
    return dropped


def isCompressed(path):
    "determine if a file appears to be compressed by extension"
    path = os.fspath(path)
    return path.endswith(".gz") or path.endswith(".bgz") or path.endswith(".bz2") or path.endswith(".Z")


def compressCmd(path, *, bgzip=False):
    """return the command to compress the path, or default if not compressed, which defaults
    to the `cat' command, so that it just gets written through"""
    path = os.fspath(path)
    if path.endswith(".Z"):
        raise PycbioException("writing compress .Z files not supported")

    if path.endswith(".bgz") or bgzip:
        return ["bgzip"]
    if path.endswith(".gz"):
        if which("pigz"):
            return ["pigz"]
        else:
            return ["gzip"]
    if path.endswith(".bz2"):
        return ["bzip2"]
    else:
        return ["cat"]

def compressBaseName(path):
    """if a file is compressed, return the path without the compressed extension"""
    path = os.fspath(path)
    if isCompressed(path):
        return osp.splitext(path)[0]
    else:
        return path

def decompressCmd(path):
    """"return the command to decompress the file to stdout, or default if not compressed, which defaults
    to the `cat' command, so that it just gets written through"""
    # FIXME: default MacOS zcat doesn't recongize .gz
    path = os.fspath(path)
    if path.endswith(".gz") or path.endswith(".bgz") or path.endswith(".Z"):
        return ["zcat"]
    elif path.endswith(".bz2"):
        return ["bzcat"]
    else:
        return ["cat"]


def opengz(fileName, mode="r", *, buffering=-1, encoding=None, errors=None, newline=None, bgzip=False):
    """open a file, if it ends in an extension indicating compression, open
    with a compression or decompression pipe.  If bgzip is specified for write,
    it is used to writing"""
    fileName = os.fspath(fileName)
    if not isCompressed(fileName):
        return open(fileName, mode, buffering=buffering, encoding=encoding, errors=errors, newline=newline)
    elif mode.startswith("r"):
        cmd = decompressCmd(fileName)
        return pipettor.Popen(cmd + [fileName], mode=mode, buffering=buffering, encoding=encoding, errors=errors, newline=newline)
    elif mode.startswith("w"):
        cmd = compressCmd(fileName, bgzip=bgzip)
        return pipettor.Popen(cmd, mode=mode, stdout=fileName, buffering=buffering, encoding=encoding, errors=errors, newline=newline)
    else:
        raise PycbioException("mode {} not support with compression for {}".format(mode, fileName))

# FIXME: make these consistent and remove redundant code.  Maybe use
# keyword for flush. Do we even need them with print function?


def prLine(fh, *objs):
    "write each str(obj) followed by a newline"
    for o in objs:
        fh.write(str(o))
    fh.write("\n")


def prsLine(fh, *objs):
    "write each str(obj), seperated by a space followed by a newline"
    n = 0
    for o in objs:
        if n > 0:
            fh.write(' ')
        fh.write(str(o))
        n += 1
    fh.write("\n")


def prOut(*objs):
    "write each str(obj) to stdout followed by a newline"
    for o in objs:
        sys.stdout.write(str(o))
    sys.stdout.write("\n")


def prErr(*objs):
    "write each str(obj) to stderr followed by a newline"
    for o in objs:
        sys.stderr.write(str(o))
    sys.stderr.write("\n")


def prsOut(*objs):
    "write each str(obj) to stdout, separating with spaces and followed by a newline"
    n = 0
    for o in objs:
        if n > 0:
            sys.stdout.write(' ')
        sys.stdout.write(str(o))
        n += 1
    sys.stdout.write("\n")


def prsfErr(*objs):
    "write each str(obj) to stderr, separating with spaces and followed by a newline and a flush"
    n = 0
    for o in objs:
        if n > 0:
            sys.stderr.write(' ')
        sys.stderr.write(str(o))
        n += 1
    sys.stderr.write("\n")
    sys.stderr.flush()


def prfErr(*objs):
    "write each str(obj) to stderr followed by a newline and a flush"
    for o in objs:
        sys.stderr.write(str(o))
    sys.stderr.write("\n")
    sys.stderr.flush()


def prsErr(*objs):
    "write each str(obj) to stderr, separating with spaces and followed by a newline"
    n = 0
    for o in objs:
        if n > 0:
            sys.stderr.write(' ')
        sys.stderr.write(str(o))
        n += 1
    sys.stderr.write("\n")


def prStrs(fh, *objs):
    "write each str(obj), with no newline"
    for o in objs:
        fh.write(str(o))


def prRow(fh, row):
    """Print a row (list or tupe) to a tab file.
    Does string conversion on each columns, None is written as empty"""
    for i in range(len(row)):
        if i > 0:
            fh.write("\t")
        fh.write(str(row[i]) if row[i] is not None else '')
    fh.write("\n")


def prRowv(fh, *objs):
    """Print a row from each argument to a tab file.
    Does string conversion on each columns,   None is written as empty"""
    prRow(fh, objs)


def fileSpecName(fspec):
    "name to use for a file spec in an error message; file objects may have one"
    if isFilePath(fspec):
        return os.fspath(fspec)
    return getattr(fspec, "name", "<file>")


class FileAccessor:
    """Context manager that opens a file (possibly compressed) if specified as
    a string, otherwise assume it is file-like and don't open/close"""
    def __init__(self, fspec, mode="r"):
        self.fspec = fspec
        self.mode = mode
        self.fh = None

    def __enter__(self):
        self.fh = opengz(self.fspec, self.mode) if isFilePath(self.fspec) else self.fspec
        return self.fh

    def __exit__(self, typ, value, traceback):
        if isFilePath(self.fspec):
            self.fh.close()


def iterLines(fspec):
    """generator over lines in file, dropping newlines.  If fspec is a string,
    open the file and close at end. Otherwise it is file-like object and will
    not be closed."""
    with FileAccessor(fspec) as fh:
        for line in fh:
            yield line.rstrip("\n")


def iterRows(fspec):
    """generator over rows in a tab-separated file.  Each line of the file is
    parsed, split into columns and returned.  If fspec is a string, open the
    file and close at end. Otherwise it is file-like object and will not be
    closed."""
    with FileAccessor(fspec) as fh:
        for line in fh:
            yield line.rstrip("\n").split("\t")


def readFileLines(fspec):
    "read lines from a open file or a file name into a list, removing the newlines"
    return [l for l in iterLines(fspec)]


def readNonCommentLines(fspec):
    """read lines from an open file or file by name into a list, removing the
    newlines, striping leading and training white space, and skipping blank
    lines and those with the first non-space character is '#'."""
    lines = []
    for line in iterLines(fspec):
        line = line.strip()
        if (len(line) > 0) and (line[0] != '#'):
            lines.append(line)
    return lines


def readLine(fh):
    "read a line from a file, dropping a newline; None on eof"
    # FIXME: delete?
    line = fh.readline()
    if len(line) == 0:
        return None
    if line[-1:] == "\n":
        line = line[:-1]
    return line

def writeLines(fspec, lines):
    "write each line, followed by a newline"
    with FileAccessor(fspec, 'w') as fh:
        for l in lines:
            fh.write(str(l))
            fh.write('\n')

def writeRows(fspec, rows):
    "write each row, joined by tabs, followed by a newline"
    with FileAccessor(fspec, 'w') as fh:
        for r in rows:
            fh.write('\t'.join([str(c) for c in r]))
            fh.write('\n')

def findTmpDir(tmpDir=None):
    """find the temporary directory to use, if tmpDir is not None, it is use"""
    if tmpDir is not None:
        return tmpDir
    tmpDir = os.getenv("TMPDIR")
    if tmpDir is not None:
        return tmpDir
    # UCSC special checks
    for tmpDir in ("/data/tmp", "/scratch/tmp", "/var/tmp", "/tmp"):
        if osp.exists(tmpDir):
            return tmpDir
    raise PycbioException("can't find a tmp directory")


def setTmpEnv(tmpDir=None):
    """Setup TMPDIR env. If tmpDir arg is not None, set TMPDIR to this value.
    If tmpDir arg is not, keep TMPDIR if set, otherwise set to value returned
    by findTmpDir"""
    os.environ["TMPDIR"] = findTmpDir(tmpDir)


def tmpFileGet(prefix=None, suffix=".tmp", tmpDir=None):
    """Obtain a tmp file with a unique name in a secure way. File
    will only be accessible to user."""
    fh = tempfile.NamedTemporaryFile(prefix=prefix, suffix=suffix,
                                     dir=findTmpDir(tmpDir), delete=False)
    fh.close()
    return fh.name


def tmpDirGet(prefix=None, suffix=".tmp", tmpDir=None):
    """Obtain a tmp directory with a unique name in a secure way.  Directory
    will only be accessible to user."""
    return tempfile.mkdtemp(prefix=prefix, suffix=suffix, dir=findTmpDir(tmpDir))

def atomicTmpFile(finalPath):
    """Return a tmp file name to use with atomicInstall.  This will be in the
    same directory as finalPath. The temporary file will have the same
    extension as finalPath.  In final path is in /dev (/dev/null,
    /dev/stdout), it is returned unchanged and atomicTmpInstall will do
    nothing.  The output directory will be created if it doesn't exist.
    Thread-safe."""
    # FIXME: not that test name for function: atomicCreate ?
    # note: this can't use tmpFileGet, since file should not be created or be private
    finalDir = osp.dirname(osp.normpath(finalPath))  # maybe empty
    if finalDir == '/dev':
        return finalPath
    finalBasename = osp.basename(finalPath)
    finalExt = osp.splitext(finalPath)[1]
    # the thread id is part of the name for the thread-safety the docstring claims:
    # host and pid alone gave two threads of one process the same temporary file
    tmpBasename = "{}.{}.{}.{}.tmp{}".format(finalBasename, socket.gethostname(), os.getpid(),
                                             threading.get_ident(), finalExt)
    tmpPath = osp.join(finalDir, tmpBasename)
    if osp.exists(tmpPath):
        os.unlink(tmpPath)
    elif finalDir != "":
        ensureDir(finalDir)
    return tmpPath

def atomicInstall(tmpPath, finalPath):
    "atomic install of tmpPath as finalPath"
    if osp.dirname(osp.normpath(finalPath)) != '/dev':
        os.rename(tmpPath, finalPath)


@contextmanager
def AtomicFileCreate(finalPath, *, keep=False):
    """Context manager to create a temporary file.  Entering returns path to
    the temporary file in the same directory as finalPath.  If the code in
    context succeeds, the file renamed to its actually name.  If an error
    occurs, the file is not installed and is removed unless keep is specified.
    The output directory will be created if it doesn't exist.  Thread-safe.
    """
    tmpPath = atomicTmpFile(finalPath)
    try:
        yield tmpPath
        atomicInstall(tmpPath, finalPath)
    except Exception:
        if not keep:
            try:
                os.unlink(tmpPath)
            except Exception:
                pass
        raise

@contextmanager
def AtomicFileOpen(finalPath, mode='w', *, buffering=-1, encoding=None,
                   errors=None, newline=None, bgzip=False, keep=False):
    """Context manager to open a temporary file.  Entering returns an open file
    object on a temporary file in the same directory as finalPath.  If the code
    in context succeeds, the file renamed to its actually name.  If an error
    occurs, the file is not installed and is removed unless keep is specified.
    The output directory will be created if it doesn't exist.  Thread-safe.
    A compression extension on finalPath is honored, as with opengz.
    """
    with AtomicFileCreate(finalPath, keep=keep) as tmpFileName:
        with opengz(tmpFileName, mode=mode, buffering=buffering, encoding=encoding,
                    errors=errors, newline=newline, bgzip=bgzip) as fh:
            yield fh

def uncompressedBase(path):
    "return the file path, removing a compression extension if it exists"
    return compressBaseName(path)


_devNullFh = None


def getDevNull():
    "get a file object open to /dev/null, caching only one instance"
    global _devNullFh
    if _devNullFh is None:
        _devNullFh = open("/dev/null", "r+")
    return _devNullFh


def _parseMd5(line):
    "parse output of openssl md5"
    m = re.match("^MD5\\((.+)\\)= ([a-f0-9]+)\\n?", line)
    return m.group(1), m.group(2)


def md5sum(filePath):
    "compute md5 on a file"
    return _parseMd5(pipettor.runout(["openssl", "md5", filePath]))[1]


_sc_arg_max = None


def getArgMax():
    global _sc_arg_max
    if _sc_arg_max is None:
        _sc_arg_max = os.sysconf("SC_ARG_MAX")
    return _sc_arg_max


def _mkMd5SumCmd(filePaths, i, maxCmdLen):
    cmd = ["openssl", "md5"]
    cmdLen = 0
    while (cmdLen < maxCmdLen) and (i < len(filePaths)):
        cmd.append(filePaths[i])
        cmdLen += len(filePaths[i])
        i += 1
    return i, cmd


def _runMd5SumCmd(cmd):
    return [_parseMd5(line) for line in pipettor.runout(cmd)[0:-1].split('\n')]


def md5sums(filePaths):
    "compute md5 on a list of files, returning list of list of (path, sum)"
    maxCmdLen = getArgMax() - 1024  # a little padding
    i = 0
    results = []
    while i < len(filePaths):
        i, cmd = _mkMd5SumCmd(filePaths, i, maxCmdLen)
        results.extend(_runMd5SumCmd(cmd))
    return results
