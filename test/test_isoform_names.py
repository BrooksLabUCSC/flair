"""
Isoform names are hashes of their structure, and isoform ids are given by name, so
isoforms of the same coordinates on different chromosomes need different names: the
single-exon ERCC spike-ins starting at 1 all had one id.
"""
import os
import subprocess
import sys
from flair.isoform_data import Isoform, Junc

JUNCS = (Junc(1100, 1500), Junc(1700, 2000))


def test_single_exon_isoforms_on_different_chromosomes():
    a = Isoform('ERCC-00003', '+', (), 0, 1006)
    b = Isoform('ERCC-00042', '+', (), 0, 1006)
    assert a.name != b.name
    assert a.name == Isoform('ERCC-00003', '+', (), 0, 1006).name
    assert a.name != Isoform('ERCC-00003', '-', (), 0, 1006).name


def test_spliced_isoforms_on_different_chromosomes():
    assert Isoform('chr1', '+', JUNCS, 800, 2200).name != Isoform('chr2', '+', JUNCS, 800, 2200).name
    # before the ends are known, by the junctions
    assert Isoform('chr1', '+', JUNCS).name != Isoform('chr2', '+', JUNCS).name


def test_names_are_the_same_between_runs():
    "Python's string hash changes between runs, so names must not depend on it"
    code = ("from flair.isoform_data import Isoform, Junc; "
            "print(Isoform('chr1', '+', (Junc(1100, 1500),), 800, 2200).name, Isoform('ERCC-00003', '+', (), 0, 1006).name)")
    names = set()
    for seed in ('1', '2'):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        names.add(subprocess.run([sys.executable, '-c', code], env=env, capture_output=True, text=True,
                                 check=True).stdout)
    assert len(names) == 1
