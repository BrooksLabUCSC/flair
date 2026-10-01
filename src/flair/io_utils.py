"""Filesystem / I/O helpers shared across FLAIR pipelines."""

import os
import tempfile


def make_temp_dir(out_prefix):
    # FIXME: use TMPDIR unless directory explicitly specified
    temp_dir = out_prefix + ".intermediate"
    try:
        os.makedirs(temp_dir, exist_ok=True)
    except OSError as exc:
        raise OSError(f"Creation of the directory `{temp_dir}' failed") from exc
    return temp_dir + '/'


def make_run_temp_dir(out_prefix, temp_dir=None):
    """Create a new directory for one run's temporary files, named for the output
    prefix, in temp_dir, or in the system temporary directory ($TMPDIR) if it is
    None.  Each run gets its own directory, so runs can share temp_dir."""
    base = temp_dir if temp_dir is not None else tempfile.gettempdir()
    try:
        os.makedirs(base, exist_ok=True)
        run_dir = tempfile.mkdtemp(prefix=os.path.basename(out_prefix) + '.', suffix='.intermediate', dir=base)
    except OSError as exc:
        raise OSError(f"Creation of a temporary directory in `{base}' failed") from exc
    return run_dir + '/'
