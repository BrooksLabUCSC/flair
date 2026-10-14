"""
Fusions between paralogs are rejected using the bundled paralog table, which
has unversioned Ensembl gene ids, whatever the annotation's id style.
"""
import os
import flair
from flair.convert_synthetic_to_genome_bed import (get_paralog_ref, paralog_groups_for_genes,
                                                   identify_fusion_problems)

PARALOG_TABLE = os.path.join(os.path.dirname(flair.__file__), 'dgd_Hsa_all_v71.tsv')

# RHD and RHCE are paralogs in the table, ACTB is in no paralog group
RHD, RHCE, ACTB = 'ENSG00000187010', 'ENSG00000188672', 'ENSG00000075624'
NAMES = {'RHD': RHD, 'RHCE': RHCE, 'ACTB': ACTB}


def _is_good_fusion(gene_a, gene_b, gene_to_paralogs, genetoname):
    fgenes = (gene_a, gene_b)
    locustopartners = {g: {frozenset([o])} for g in fgenes for o in fgenes if o != g}
    return identify_fusion_problems(fgenes, locustopartners, 4, genetoname, gene_to_paralogs,
                                    ['chr1', 'chr1'], 5)


def _gene_ids(version):
    "gene name -> id, with a version suffix as GENCODE writes them, or without as Ensembl does"
    return {name: f'{gid}.{version}' if version else gid for name, gid in NAMES.items()}


def test_versioned_ids_find_paralog_group():
    table = get_paralog_ref(PARALOG_TABLE)
    ids = _gene_ids(7)
    groups = paralog_groups_for_genes(table, ids.values())
    assert groups[ids['RHD']] == groups[ids['RHCE']] == table[RHD]
    assert ids['ACTB'] not in groups


def test_paralog_fusion_rejected_with_versioned_and_unversioned_ids():
    table = get_paralog_ref(PARALOG_TABLE)
    for version in (None, 7):
        ids = _gene_ids(version)
        groups = paralog_groups_for_genes(table, ids.values())
        genetoname = {gid: name for name, gid in ids.items()}
        assert not _is_good_fusion(ids['RHD'], ids['RHCE'], groups, genetoname)
        assert _is_good_fusion(ids['RHD'], ids['ACTB'], groups, genetoname)


def test_versioned_ids_unmatched_without_remapping():
    "the table used as is, as before, misses GENCODE ids and lets a paralog fusion through"
    table = get_paralog_ref(PARALOG_TABLE)
    ids = _gene_ids(7)
    genetoname = {gid: name for name, gid in ids.items()}
    assert _is_good_fusion(ids['RHD'], ids['RHCE'], table, genetoname)


def test_genes_sharing_a_name_stem_are_not_rejected():
    "BCAS4--BCAS3, a known MCF7 fusion, share the first 4 letters of their names but aren't paralogs"
    bcas4, bcas3 = 'ENSG00000124243.20', 'ENSG00000141376.24'
    genetoname = {bcas4: 'BCAS4', bcas3: 'BCAS3'}
    assert _is_good_fusion(bcas4, bcas3, {}, genetoname)


def test_genes_of_one_name_are_rejected():
    "two gene ids of one name, as a gene's X and Y copies, are not a fusion"
    x_copy, y_copy = 'ENSG00000169100.15', 'ENSG00000169100.15_PAR_Y'
    genetoname = {x_copy: 'SLC25A6', y_copy: 'SLC25A6'}
    assert not _is_good_fusion(x_copy, y_copy, {}, genetoname)
