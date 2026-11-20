"""Reusing a file an earlier run left.

A file is reused only when the inputs and the options that produced it are those of
this run, recorded beside it as it was written.  Reuse on the file's existence
alone gives the earlier run's result under this run's options, with nothing in the
output to state that it happened.

A file is reusable when its record is there and matches.  A file an interrupted run
left has no record beside it, since the record is written last, and a file made from
other inputs has one that differs, so both are made again.
"""

import json
import logging
import os
from flair import VERSION
from flair.pycbio.sys import fileOps

PROVENANCE_SUFFIX = '.provenance.json'


def _file_stamp(path):
    """a file's identity, for a provenance record: its path, its size and the time it
    was last modified.  A file replaced by one of the same size within the same
    second is not told apart, which is the limit of a stamp of this kind."""
    return [os.path.abspath(path), os.path.getsize(path), int(os.path.getmtime(path))]


def _stamps(files):
    "the stamps of the input files, each a path or a list of them, and None for none"
    paths = []
    for entry in files:
        paths += [entry] if isinstance(entry, str) else list(entry or ())
    return sorted(_file_stamp(path) for path in paths if path)


def provenance(*, files=(), options):
    """What a file was produced from: the inputs, the options that change its
    contents, and the version of flair that wrote it."""
    return {'version': VERSION, 'files': _stamps(files), 'options': dict(options)}


def _provenance_path(path):
    return path + PROVENANCE_SUFFIX


def _recorded_provenance(path):
    "the record beside a file, or None when it has none or cannot be read"
    try:
        with open(_provenance_path(path)) as fh:
            return json.load(fh)
    except (OSError, ValueError):
        return None


def record_provenance(path, prov):
    """Record what produced a file, once the file itself is complete: the record is
    what a later run tests, so it must not be there before the file it describes."""
    with fileOps.AtomicFileCreate(_provenance_path(path)) as tmp_path:
        with open(tmp_path, 'w') as fh:
            json.dump(prov, fh, indent=2, sort_keys=True)
            fh.write('\n')


def reusable(path, prov):
    "whether a file an earlier run left was produced from these inputs and options"
    return os.path.exists(path) and (_recorded_provenance(path) == prov)


def reuse_file(path, prov, *, what):
    """Whether to reuse the file an earlier run left, reporting what is reused and
    what is made again: a run whose result is partly an earlier run's states which
    parts, and one that redoes work the user expected to be reused states why."""
    if not os.path.exists(path):
        return False
    if _recorded_provenance(path) != prov:
        logging.info('%s: from other inputs or options, or from an interrupted run, '
                     'so made again: %s', what, path)
        return False
    logging.info('reusing %s of an earlier run: %s', what, path)
    return True
