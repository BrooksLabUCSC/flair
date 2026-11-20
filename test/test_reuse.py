"""Reusing a file an earlier run left: only where the inputs and the options that
produced it are this run's."""
import logging
import os
import pytest
from flair.reuse import provenance, record_provenance, reusable, reuse_file


@pytest.fixture
def made_file(tmp_path):
    "a file an earlier run wrote from one input, with its provenance recorded"
    source = tmp_path / 'reads.bam'
    source.write_text('reads')
    made = tmp_path / 'chimeric.bam'
    made.write_text('result')
    prov = provenance(files=(str(source),), options={'min_fragment_size': 40})
    record_provenance(str(made), prov)
    return str(source), str(made), prov


def test_the_same_inputs_and_options_are_reused(made_file):
    _, made, prov = made_file
    assert reusable(made, prov)


def test_other_options_are_not_reused(made_file):
    source, made, _ = made_file
    other = provenance(files=(str(source),), options={'min_fragment_size': 20})
    assert not reusable(made, other)


def test_another_input_file_is_not_reused(made_file, tmp_path):
    _, made, _ = made_file
    other_source = tmp_path / 'other.bam'
    other_source.write_text('other reads')
    other = provenance(files=(str(other_source),), options={'min_fragment_size': 40})
    assert not reusable(made, other)


def test_an_input_written_again_is_not_reused(made_file):
    "the input is the same path, with other contents: its size and time are recorded"
    source, made, _ = made_file
    with open(source, 'w') as fh:
        fh.write('other reads, of another length')
    os.utime(source, (0, 0))
    remade = provenance(files=(source,), options={'min_fragment_size': 40})
    assert not reusable(made, remade)


def test_a_file_of_an_interrupted_run_is_not_reused(made_file):
    "the record is written after the file, so a file without one was never finished"
    _, made, prov = made_file
    os.remove(made + '.provenance.json')
    assert not reusable(made, prov)


def test_a_file_that_is_not_there_is_not_reused(made_file, tmp_path):
    _, _, prov = made_file
    assert not reusable(str(tmp_path / 'nothing.bam'), prov)


def test_the_reuse_is_reported(made_file, caplog):
    _, made, prov = made_file
    with caplog.at_level(logging.INFO):
        assert reuse_file(made, prov, what='the chimeric reads')
    assert 'reusing the chimeric reads of an earlier run' in caplog.text


def test_making_it_again_is_reported(made_file, caplog):
    source, made, _ = made_file
    other = provenance(files=(source,), options={'min_fragment_size': 20})
    with caplog.at_level(logging.INFO):
        assert not reuse_file(made, other, what='the chimeric reads')
    assert 'from other inputs or options' in caplog.text


def test_the_order_of_the_input_files_does_not_matter(tmp_path):
    first, second = tmp_path / 'a.bam', tmp_path / 'b.gtf'
    first.write_text('a')
    second.write_text('b')
    assert (provenance(files=(str(first), str(second)), options={}) ==
            provenance(files=(str(second), str(first)), options={}))


def test_options_that_were_not_given_are_left_out(tmp_path):
    source = tmp_path / 'a.bam'
    source.write_text('a')
    prov = provenance(files=(str(source), None, ()), options={})
    assert len(prov['files']) == 1
