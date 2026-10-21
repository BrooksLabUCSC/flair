"""Filesystem / I/O helpers shared across FLAIR pipelines."""

import os
import tempfile
from flair import FlairInputDataError


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


def read_manifest_rows(path, field_counts, line_format, file_fields):
    """A manifest's rows, each a list of its tab-separated fields, for a line each but
    blank lines and '#' comments.  A line has one of field_counts fields, none of them
    empty, as line_format describes for errors; the fields at the indexes of
    file_fields, index -> what the file is, name files that must exist; and the first
    field, the sample's name, is unique."""
    rows = []
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip('\r\n')
            if line.strip() == '' or line.startswith('#'):
                continue
            fields = line.split('\t')
            where = f'{path}:{lineno}'
            if len(fields) not in field_counts or '' in fields:
                spaces = '; fields are separated by tabs, not spaces' if len(fields) == 1 and ' ' in line.strip() else ''
                raise FlairInputDataError(f'{where}: a manifest line is {line_format}, tab-separated{spaces}: {line!r}')
            for index, what in file_fields.items():
                if index < len(fields) and not os.path.exists(fields[index]):
                    raise FlairInputDataError(f'{where}: {what} {fields[index]} does not exist')
            rows.append(fields)
    if len(rows) == 0:
        raise FlairInputDataError(f'{path}: manifest has no samples')
    names = [row[0] for row in rows]
    duplicated = sorted({name for name in names if names.count(name) > 1})
    if duplicated:
        raise FlairInputDataError(f'{path}: samples named more than once: {", ".join(duplicated)}')
    return rows
