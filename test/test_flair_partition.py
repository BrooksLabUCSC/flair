import io
from flair.flair_partition import gtf_to_file, read_ranges, convert_inputs, bam_chroms, bam_to_file
from flair.gtf_to_bed import gtf_to_bed


def _bed3(path):
    with open(path) as fh:
        return sorted(tuple(line.rstrip('\n').split('\t')[:3]) for line in fh)


def test_gtf_spans_match_gtf_to_bed(tmp_path):
    gtf = 'input/seg1.gencodeV47.gtf'
    gtf_to_bed(str(tmp_path / 'full.bed'), gtf, include_gene=True)
    gtf_to_file(gtf, str(tmp_path / 'spans.bed'))
    assert _bed3(tmp_path / 'spans.bed') == _bed3(tmp_path / 'full.bed')


def test_read_ranges_only_needs_three_columns():
    fh = io.StringIO('track name=x\n#comment\n\nchr1\t10\t20\n'
                     'chr2\t5\t8\tname\t0\t+\n')
    assert [(r.chrom, r.chromStart, r.chromEnd) for r in read_ranges(fh)] == [('chr1', 10, 20), ('chr2', 5, 8)]


def test_indexed_bam_split_by_chrom(tmp_path):
    # converting by chrom gives the same alignments as converting the whole BAM
    bam = 'output/test-align.bam'
    assert sorted(bam_chroms(bam)) == ['chr12', 'chr17', 'chr20']
    split = convert_inputs([], [bam], [], 3, str(tmp_path))
    assert len(split) == 3
    bam_to_file(bam, str(tmp_path / 'whole.bed'))
    split_lines = sorted(line for f in split for line in open(f))
    assert split_lines == sorted(open(tmp_path / 'whole.bed'))


def test_unindexed_sam_converted_whole(tmp_path):
    assert bam_chroms('input/tiny.sam') is None
    assert len(convert_inputs([], ['input/tiny.sam'], [], 3, str(tmp_path))) == 1


def test_single_thread_converts_bam_whole(tmp_path):
    # one bamtobed process, rather than samtools view and bamtobed per chrom
    assert len(convert_inputs([], ['output/test-align.bam'], [], 1, str(tmp_path))) == 1


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
