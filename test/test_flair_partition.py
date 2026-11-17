import io
import os
import pytest
from flair import FlairInputDataError, SeqRange
from flair.flair_partition import (gtf_to_file, read_ranges, convert_inputs, bam_chroms, bam_to_file,
                                   check_bam_indexed, closed_partitions, scan_partitions, RegionSet)
from flair.gtf_to_bed import gtf_to_bed
from flair.partition_runner import run_flair_partition


def _bed3(path):
    with open(path) as fh:
        return sorted(tuple(line.rstrip('\n').split('\t')[:3]) for line in fh)


def test_gtf_spans_match_gtf_to_bed(tmp_path):
    gtf = 'input/seg1.gencode_v47.gtf'
    gtf_to_bed(str(tmp_path / 'full.bed'), gtf, include_gene=True)
    gtf_to_file(gtf, str(tmp_path / 'spans.bed'))
    assert _bed3(tmp_path / 'spans.bed') == _bed3(tmp_path / 'full.bed')


def test_read_ranges_only_needs_three_columns():
    fh = io.StringIO('track name=x\n#comment\n\nchr1\t10\t20\n'
                     'chr2\t5\t8\tname\t0\t+\n')
    assert [(r.chrom, r.chromStart, r.chromEnd) for r in read_ranges(fh)] == [('chr1', 10, 20), ('chr2', 5, 8)]


def test_indexed_bam_split_by_chrom(tmp_path, make_output_dir):
    # converting by chrom gives the same alignments as converting the whole BAM
    bam = os.path.join(make_output_dir, 'basic.reads.bam')
    assert sorted(bam_chroms(bam)) == ['chr12', 'chr17', 'chr20']
    split = convert_inputs([], [bam], [], 3, str(tmp_path))
    assert len(split) == 3
    bam_to_file(bam, str(tmp_path / 'whole.bed'))
    split_lines = sorted(line for f in split for line in open(f))
    assert split_lines == sorted(open(tmp_path / 'whole.bed'))


def test_unindexed_sam_converted_whole(tmp_path):
    assert bam_chroms('input/tiny.sam') is None
    assert len(convert_inputs([], ['input/tiny.sam'], [], 3, str(tmp_path))) == 1


def test_single_thread_converts_bam_whole(tmp_path, make_output_dir):
    # one bamtobed process, rather than samtools view and bamtobed per chrom
    bam = os.path.join(make_output_dir, 'basic.reads.bam')
    assert len(convert_inputs([], [bam], [], 1, str(tmp_path))) == 1


def test_thread_slots_limit_concurrency():
    import threading
    import time
    from concurrent.futures import ThreadPoolExecutor
    from flair.flair_partition import _ThreadSlots
    slots, lock, running, peak = _ThreadSlots(4), threading.Lock(), [0], [0]

    def job(nslots):
        with lock:
            running[0] += nslots
            peak[0] = max(peak[0], running[0])
        time.sleep(0.02)
        with lock:
            running[0] -= nslots

    with ThreadPoolExecutor(max_workers=4) as executor:
        for f in [executor.submit(slots.run, n, job, n) for n in (1, 2, 2, 2, 1, 2, 1, 1)]:
            f.result()
    assert peak[0] <= 4


def _partition_spans(partitions):
    return sorted((bed.chrom, bed.chromStart, bed.chromEnd, count) for bed, count in partitions)


def test_region_set_merges_and_finds_overlaps():
    regions = RegionSet([SeqRange('chr1', 100, 200), SeqRange('chr1', 150, 300), SeqRange('chr2', 1, 5)])
    assert regions.regions == [SeqRange('chr1', 100, 300), SeqRange('chr2', 1, 5)]
    assert regions.overlaps('chr1', 250, 260)
    assert not regions.overlaps('chr1', 300, 400)
    assert not regions.overlaps('chr3', 100, 200)


def test_region_grows_over_a_chain_of_items(tmp_path):
    # the region grows over the items reaching it, then over the items reaching
    # those, so it ends as the partition of the whole input, 100-400
    chain_bed = tmp_path / 'chain.bed'
    chain_bed.write_text('chr1\t100\t200\nchr1\t190\t300\nchr1\t290\t400\nchr1\t1000\t1100\n')
    work_dir = tmp_path / 'work'
    work_dir.mkdir()
    regions = RegionSet([SeqRange('chr1', 100, 110)])
    partitions = closed_partitions([str(chain_bed)], [], [], 1, 0, regions, str(work_dir))
    assert _partition_spans(partitions) == [('chr1', 100, 400, 3)]


def test_region_partitions_are_those_of_the_whole_input(tmp_path, make_output_dir):
    # a region inside a gene gives the partition, and its item count, that
    # partitioning the whole BAM gives
    bam = os.path.join(make_output_dir, 'basic.reads.bam')
    full_dir, region_dir = tmp_path / 'full', tmp_path / 'region'
    full_dir.mkdir()
    region_dir.mkdir()
    full = _partition_spans(scan_partitions([], [bam], [], 2, 0, 0, None, str(full_dir)))
    regions = RegionSet([SeqRange('chr17', 64490000, 64490100)])
    partitions = closed_partitions([], [bam], [], 2, 0, regions, str(region_dir))
    assert _partition_spans(partitions) == [part for part in full if part[0] == 'chr17']


def test_region_with_no_items_gives_no_partitions(tmp_path, make_output_dir):
    bam = os.path.join(make_output_dir, 'basic.reads.bam')
    regions = RegionSet([SeqRange('chr12', 60000000, 60010000)])
    assert closed_partitions([], [bam], [], 2, 0, regions, str(tmp_path)) == []


def test_bam_regions_hold_the_alignments_overlapping_them(tmp_path, make_output_dir):
    bam = os.path.join(make_output_dir, 'basic.reads.bam')
    regions = RegionSet([SeqRange('chr20', 35550000, 35551000)])
    regions_bed = str(tmp_path / 'regions.bed')
    regions.write_bed(regions_bed)
    bam_to_file(bam, str(tmp_path / 'regions.out.bed'), None, regions_bed)
    bam_to_file(bam, str(tmp_path / 'whole.bed'))
    with open(tmp_path / 'whole.bed') as fh:
        wanted = sorted(tuple(line.split('\t')[:3]) for line in fh
                        if regions.overlaps(line.split('\t')[0], int(line.split('\t')[1]), int(line.split('\t')[2])))
    assert _bed3(tmp_path / 'regions.out.bed') == wanted


def test_unindexed_bam_cannot_be_read_by_region():
    with pytest.raises(FlairInputDataError, match="indexed BAM or CRAM"):
        check_bam_indexed('input/tiny.sam')


def _alignment_spans(bams, tmp_path):
    spans = []
    for num, bam in enumerate(bams):
        bed = str(tmp_path / f'{num}.bed')
        bam_to_file(bam, bed)
        with open(bed) as fh:
            spans += [(f[0], int(f[1]), int(f[2])) for f in (line.split('\t') for line in fh)]
    return spans


def _containing(regions, chrom, start, end):
    return [r for r in regions if r.name == chrom and r.start < end and start < r.end]


def test_no_alignment_crosses_a_partition(tmp_path, make_output_dir):
    # what flair spliceevents corrects and counts reads by: a read is in one partition
    bams = [os.path.join(make_output_dir, 'basic.reads.bam'), os.path.join(make_output_dir, 'seg1.reads.bam')]
    regions, counts = run_flair_partition(bams, 'input/basic.annotation.gtf', 2)
    assert len(regions) == len(counts) > 3
    for chrom, start, end in _alignment_spans(bams, tmp_path):
        holding = _containing(regions, chrom, start, end)
        assert len(holding) == 1
        assert holding[0].start <= start and end <= holding[0].end


def test_partitions_are_merged_up_to_the_item_count(tmp_path, make_output_dir):
    # the two chr20 loci of the two BAMs are separate partitions, merged into one
    bams = [os.path.join(make_output_dir, 'basic.reads.bam'), os.path.join(make_output_dir, 'seg1.reads.bam')]
    clusters = [r for r in run_flair_partition(bams, None, 2)[0] if r.name == 'chr20']
    merged = [r for r in run_flair_partition(bams, None, 2, min_partition_items=100000)[0] if r.name == 'chr20']
    assert len(clusters) > 1
    assert len(merged) == 1
    assert (merged[0].start, merged[0].end) == (clusters[0].start, clusters[-1].end)
