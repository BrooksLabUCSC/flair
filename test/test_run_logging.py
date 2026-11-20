"""What a run reports about itself: the record format, the progress of a parallel
stage, and the level passed to the programs a run starts."""
import logging
import re
import time
import pytest
from flair.run_logging import (FlairLogFormatter, child_log_level_opt, plural,
                               report_task_done, start_task_progress)

TIME_STAMP = r'\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}'


def format_record(level, msg):
    record = logging.LogRecord('flair', level, __file__, 1, msg, (), None)
    return FlairLogFormatter().format(record)


def test_info_record_has_the_time_and_no_level():
    assert re.fullmatch(TIME_STAMP + ' loading genome', format_record(logging.INFO, 'loading genome'))


def test_warning_record_is_marked():
    assert re.fullmatch(TIME_STAMP + ' WARNING: no isoforms', format_record(logging.WARNING, 'no isoforms'))


def test_error_record_is_marked():
    assert re.fullmatch(TIME_STAMP + ' ERROR: bad', format_record(logging.ERROR, 'bad'))


def test_the_error_that_ends_a_run_names_its_own_exception():
    "the CRITICAL record of cli.ErrorHandler, whose message already names the error"
    assert re.fullmatch(TIME_STAMP + ' FlairInputDataError: no such file',
                        format_record(logging.CRITICAL, 'FlairInputDataError: no such file'))


@pytest.mark.parametrize('count, expect', [(0, '0 regions'), (1, '1 region'), (2, '2 regions')])
def test_plural_counts(count, expect):
    assert plural(count, 'region') == expect


def test_plural_with_an_irregular_noun():
    assert plural(2, 'batch', 'batches') == '2 batches'


def test_task_progress_counts_from_one(caplog):
    start_task_progress(2)
    with caplog.at_level(logging.INFO):
        report_task_done('partition', 'chr1:0-100, 7 items', time.time())
        report_task_done('partition', 'chr2:0-100, 9 items', time.time())
    assert [record.getMessage() for record in caplog.records] == [
        'partition 1 of 2: chr1:0-100, 7 items, 0.0 s',
        'partition 2 of 2: chr2:0-100, 9 items, 0.0 s']


def test_child_log_level_is_the_level_of_the_run():
    logging.getLogger().setLevel(logging.DEBUG)
    try:
        assert child_log_level_opt() == f'--log-level={logging.DEBUG}'
    finally:
        logging.getLogger().setLevel(logging.WARNING)
