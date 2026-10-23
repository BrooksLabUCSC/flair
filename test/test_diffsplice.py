"""
Tests for flair_diffSplice module.
"""
import pytest
from flair import FlairInputDataError
from flair.counts_matrix_tsv import (MAX_REPORTED_IDS, CountsRow,
                                     write_counts_matrix as write_counts_rows)
from flair.flair_diffSplice import check_isoform_ids

def write_bed(path, names):
    with open(path, 'w') as fh:
        for num, name in enumerate(names):
            start = 1000 + (100 * num)
            fh.write('\t'.join(['chr1', str(start), str(start + 50), name, '0', '+']) + '\n')

def write_counts_matrix(path, iso_ids):
    write_counts_rows(path, ['s1_ctl_b1', 's2_test_b1'],
                      [CountsRow('gene1', iso_id, ['1', '2']) for iso_id in iso_ids])

def write_inputs(tmp_path, bed_names, counts_ids):
    bed = str(tmp_path / 'isoforms.bed')
    counts = str(tmp_path / 'x.counts.tsv')
    write_bed(bed, bed_names)
    write_counts_matrix(counts, counts_ids)
    return bed, counts

def test_check_isoform_ids_all_present(tmp_path):
    bed, counts = write_inputs(tmp_path, ['iso1_gene1', 'iso2_gene1'],
                               ['iso1_gene1', 'iso2_gene1'])
    check_isoform_ids(bed, counts)

def test_check_isoform_ids_extra_bed_records(tmp_path):
    "a BED isoform with no counts is not an error, it contributes no counts"
    bed, counts = write_inputs(tmp_path, ['iso1_gene1', 'iso2_gene1'], ['iso1_gene1'])
    check_isoform_ids(bed, counts)

def test_check_isoform_ids_missing_from_bed(tmp_path):
    bed, counts = write_inputs(tmp_path, ['iso1_gene1'], ['iso1_gene1', 'iso2_gene1'])
    with pytest.raises(FlairInputDataError, match=r"1 isoform ids .* iso2_gene1"):
        check_isoform_ids(bed, counts)

def test_check_isoform_ids_missing_list_elided(tmp_path):
    missing = ['iso%d_gene1' % num for num in range(MAX_REPORTED_IDS + 5)]
    bed, counts = write_inputs(tmp_path, ['iso0_gene1'], missing)
    with pytest.raises(FlairInputDataError, match=r"14 isoform ids .*, \.\.\.;"):
        check_isoform_ids(bed, counts)
