"""The sample columns of a flair quantify counts matrix.

flair quantify writes a sample info TSV beside the counts matrix, one row per counts
column, naming each sample's condition and batch.  Readers take the fields from there
rather than from the column names.

Counts matrices written before that file existed name each column
sample_condition_batch, joining the first three manifest fields with an underscore,
so condition and batch have to be recovered by splitting the name.  That is why the
manifest documentation forbids underscores in those fields.  Such a matrix is still
read, by falling back to that parsing.

This is also the one place that decides which two conditions are being compared and
which of them is the reference.
"""
import os
import csv
import logging
from collections import namedtuple
from flair import FlairInputDataError

CONDITION_FIELD = 1
BATCH_FIELD = -1

SAMPLE_INFO_COLUMNS = ('sample_id', 'condition', 'batch')

class SampleInfo(namedtuple('SampleInfo', SAMPLE_INFO_COLUMNS)):
    "one counts matrix column: which sample it holds and how that sample was grouped"
    __slots__ = ()

def sample_info_path(counts_matrix_tsv):
    "flair quantify writes <prefix>.counts.tsv beside <prefix>.sample_info.tsv"
    for suffix in ('.counts.tsv', '.tsv'):
        if counts_matrix_tsv.endswith(suffix):
            return counts_matrix_tsv[:-len(suffix)] + '.sample_info.tsv'
    return counts_matrix_tsv + '.sample_info.tsv'

def write_sample_info(sample_info_tsv, sample_infos):
    with open(sample_info_tsv, 'w') as fh:
        writer = csv.writer(fh, delimiter='\t', dialect='unix', quoting=csv.QUOTE_NONE)
        writer.writerow(SAMPLE_INFO_COLUMNS)
        writer.writerows(sample_infos)

def _read_sample_info_rows(sample_info_tsv):
    with open(sample_info_tsv) as fh:
        reader = csv.reader(fh, delimiter='\t')
        columns = tuple(next(reader))
        if columns != SAMPLE_INFO_COLUMNS:
            raise FlairInputDataError(
                f"{sample_info_tsv}: expected columns {', '.join(SAMPLE_INFO_COLUMNS)}, "
                f"found: {', '.join(columns)}")
        return [SampleInfo(*row) for row in reader]

def _check_sample_info(sample_infos, sample_columns, sample_info_tsv, counts_matrix_tsv):
    "the sample info rows must describe the counts columns, in the same order"
    named = [si.sample_id for si in sample_infos]
    if named != sample_columns:
        raise FlairInputDataError(
            f"{sample_info_tsv} does not describe the columns of {counts_matrix_tsv}: "
            f"it names {', '.join(named)}, the counts columns are {', '.join(sample_columns)}")

def _sample_info_from_columns(sample_columns, counts_matrix_tsv):
    "recover the fields from sample_condition_batch column names"
    conditions, batches = parse_sample_fields(sample_columns, counts_matrix_tsv)
    return [SampleInfo(col, condition, batch)
            for col, condition, batch in zip(sample_columns, conditions, batches)]

def read_sample_info(counts_matrix_tsv):
    """The sample of each counts column and how it was grouped, taken from the sample
    info file that flair quantify writes, or from the column names when a counts
    matrix predates that file."""
    sample_columns = read_sample_columns(counts_matrix_tsv)
    sample_info_tsv = sample_info_path(counts_matrix_tsv)
    if not os.path.exists(sample_info_tsv):
        return _sample_info_from_columns(sample_columns, counts_matrix_tsv)
    sample_infos = _read_sample_info_rows(sample_info_tsv)
    _check_sample_info(sample_infos, sample_columns, sample_info_tsv, counts_matrix_tsv)
    return sample_infos

def read_sample_columns(counts_matrix_tsv):
    "the sample column names, in column order, without the leading id column"
    with open(counts_matrix_tsv) as fh:
        return fh.readline().split()[1:]

def parse_sample_fields(sample_columns, counts_matrix_tsv):
    "the condition and the batch of each sample column"
    try:
        conditions = [col.split('_')[CONDITION_FIELD] for col in sample_columns]
        batches = [col.split('_')[BATCH_FIELD] for col in sample_columns]
    except IndexError as ex:
        raise FlairInputDataError(
            f"{counts_matrix_tsv}: counts columns must be named sample_condition_batch, "
            f"found: {' '.join(sample_columns)}") from ex
    return conditions, batches

def condition_column_indexes(conditions, condition):
    "indexes of the sample columns belonging to one condition"
    return [i for i, c in enumerate(conditions) if c == condition]

def _check_named_conditions(present, condition_a, condition_b, counts_matrix_tsv):
    for opt, name in (('--condition_a', condition_a), ('--condition_b', condition_b)):
        if name not in present:
            raise FlairInputDataError(
                f"{opt} {name} is not a condition in {counts_matrix_tsv}, which has: "
                f"{', '.join(present)}")
    if condition_a == condition_b:
        raise FlairInputDataError("--condition_a and --condition_b must name different conditions")

def _default_conditions(present, counts_matrix_tsv):
    if len(present) != 2:
        raise FlairInputDataError(
            f"{counts_matrix_tsv} has {len(present)} conditions ({', '.join(present)}); "
            "name the two to compare with --condition_a and --condition_b")
    return present[0], present[1]

def select_condition_pair(conditions, condition_a, condition_b, counts_matrix_tsv):
    """The two conditions to compare, with condition_a the reference that fold
    changes are measured against.  With neither named, the two conditions are taken
    in sorted order rather than in column order, so reordering the columns of the
    counts matrix cannot change the result."""
    present = sorted(set(conditions))
    if condition_a and condition_b:
        _check_named_conditions(present, condition_a, condition_b, counts_matrix_tsv)
    elif condition_a or condition_b:
        raise FlairInputDataError("--condition_a and --condition_b must both be given, "
                                  "or both left out to take the two conditions in sorted order")
    else:
        condition_a, condition_b = _default_conditions(present, counts_matrix_tsv)
        logging.info(f"comparing {condition_b} against reference {condition_a}, "
                     "the two conditions in sorted order; name them with "
                     "--condition_a and --condition_b to compare in the other direction")
    return condition_a, condition_b
