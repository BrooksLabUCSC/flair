"""Checks on a matrix before it is handed to DRIMSeq.

DRIMSeq indexes its rows by feature_id and refuses a repeat, reporting it as a
failed stopifnot that names neither the file nor a value.  Both commands that run
DRIMSeq check here first: flair diffsplice with an events quant TSV, flair diffexp
with its DRIMSeq counts matrix.
"""
from flair import FlairError
from flair.pycbio.tsv import TsvReader

FEATURE_ID_COLUMN = 'feature_id'

def read_feature_ids(drimseq_matrix_tsv):
    "the feature_id column, in file order"
    with TsvReader(drimseq_matrix_tsv, defaultColType=str) as reader:
        return [row.feature_id for row in reader]

def duplicate_feature_ids(feature_ids):
    "the ids that appear more than once, sorted"
    seen = set()
    return sorted({fid for fid in feature_ids if fid in seen or seen.add(fid)})

def name_examples(names, count=3):
    return ', '.join(sorted(set(names))[:count])

def check_feature_ids_unique(drimseq_matrix_tsv):
    "raise if any feature_id repeats, naming the file and the first few repeats"
    dups = duplicate_feature_ids(read_feature_ids(drimseq_matrix_tsv))
    if dups:
        raise FlairError(
            f"{drimseq_matrix_tsv} has {len(dups)} repeated {FEATURE_ID_COLUMN}: "
            f"{name_examples(dups)}. Each one must appear once; please report this "
            "with the input that produced it")
