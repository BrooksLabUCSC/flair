"""
flair fusion recognises nuclear chromosomes and non-genic loci by their names in
UCSC/GENCODE, RefSeq and Ensembl genomes, not only by a chr prefix.
"""
from flair.convert_synthetic_to_genome_bed import is_reference_chromosome, is_nongenic_locus


def test_nuclear_chromosomes_in_each_naming():
    for chrom in ('chr12', 'chrX', 'NC_000012.12', 'NC_000023.11', '12', 'X', 'Y'):
        assert is_reference_chromosome(chrom), chrom


def test_mitochondria_and_scaffolds_are_not():
    for chrom in ('chrM', 'chrMT', 'MT', 'M', 'NC_012920.1', 'GL000009.2', 'KI270706.1',
                  'NT_187361.1', 'NW_025791756.1'):
        assert not is_reference_chromosome(chrom), chrom


def test_nongenic_loci_in_each_naming():
    for name in ('chr12:25205245-25227413', 'NC_000012.12:25205245-25227413', '12:100-200'):
        assert is_nongenic_locus(name), name
    for name in ('ENSG00000133703.12', 'KRAS', 'ENSG00000125991.19_PAR_Y'):
        assert not is_nongenic_locus(name), name
