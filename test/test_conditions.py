"""
Tests for conditions module.
"""
import pytest
from flair import FlairInputDataError
from flair.conditions import select_condition_pair, condition_column_indexes

def test_select_condition_pair_sorted_default():
    assert select_condition_pair(['test', 'ctl', 'ctl'], '', '', 'x.tsv') == ('ctl', 'test')

def test_select_condition_pair_ignores_column_order():
    "the whole point: the same conditions in any column order give the same pair"
    assert (select_condition_pair(['ctl', 'test', 'test', 'ctl'], '', '', 'x.tsv') ==
            select_condition_pair(['test', 'test', 'ctl', 'ctl'], '', '', 'x.tsv'))

def test_select_condition_pair_named():
    assert select_condition_pair(['ctl', 'test'], 'test', 'ctl', 'x.tsv') == ('test', 'ctl')

def test_select_condition_pair_unknown_name():
    with pytest.raises(FlairInputDataError, match="is not a condition in"):
        select_condition_pair(['ctl', 'test'], 'nope', 'ctl', 'x.tsv')

def test_select_condition_pair_needs_both():
    with pytest.raises(FlairInputDataError, match="must both be given"):
        select_condition_pair(['ctl', 'test'], 'ctl', '', 'x.tsv')

def test_select_condition_pair_needs_exactly_two():
    with pytest.raises(FlairInputDataError, match="has 3 conditions"):
        select_condition_pair(['a', 'b', 'c'], '', '', 'x.tsv')

def test_condition_column_indexes():
    assert condition_column_indexes(['ctl', 'test', 'ctl'], 'ctl') == [0, 2]
    assert condition_column_indexes(['ctl', 'test', 'ctl'], 'test') == [1]
