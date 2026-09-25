"""The sample columns of a flair quantify counts matrix.

flair quantify names each counts column sample_condition_batch, joining the first
three manifest fields with an underscore, so condition and batch are not columns of
their own and every reader has to take them apart again.  This module is the one
place that does that, and the one place that decides which two conditions are being
compared and which of them is the reference.
"""
import logging
from flair import FlairInputDataError

CONDITION_FIELD = 1
BATCH_FIELD = -1

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
