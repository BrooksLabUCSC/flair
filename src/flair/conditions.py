"""Which two conditions a differential comparison uses, and which of them is the
reference that changes are measured against.

Conditions come from the sample info of a counts matrix; see sample_info_tsv.
This is the one place that chooses the pair, so diffexp and diffsplice choose it
the same way.
"""
import logging
from flair import FlairInputDataError

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
