import logging
from flair import SeqRange
from flair.flair_transcriptome import (report_partition_sizes, LARGE_PARTITION_ITEMS,
                                       BYTES_PER_PARTITION_ITEM)

REGIONS = [SeqRange('chr1', 0, 1000), SeqRange('chr2', 0, 2000), SeqRange('chr3', 0, 3000)]


def test_largest_partition_reported(caplog):
    with caplog.at_level(logging.INFO):
        report_partition_sizes(REGIONS, [10, 300, 20])
    assert 'largest partition chr2:0-2000, 300 reads and annotated transcripts' in caplog.text
    assert caplog.records[-1].levelno == logging.INFO


def test_large_partition_warns_with_what_it_needs(caplog):
    items = LARGE_PARTITION_ITEMS
    with caplog.at_level(logging.INFO):
        report_partition_sizes(REGIONS, [1, items, 2])
    warning = caplog.records[-1]
    assert warning.levelno == logging.WARNING
    assert f'{items * BYTES_PER_PARTITION_ITEM / 1e9:.0f} GB' in warning.message
    assert '--threads' in warning.message


def test_no_partitions_reports_nothing(caplog):
    with caplog.at_level(logging.INFO):
        report_partition_sizes([], [])
    assert caplog.records == []
