import pysam
from flair.count_sam_transcripts import check_splicesites, process_cigar

# a transcript of two 100-base exons, its splice site at offset 100; coveredpos has
# 1 for a matched position, 0 for unmatched, and 1 + insert size before an insertion,
# indexed by 0-based transcript position, and alignments span [tstart, tend)
EXONS = [100, 100]


def test_coveredpos_is_indexed_by_transcript_position():
    # 10 aligned bases from 0-based transcript position 5, the third a mismatch
    _, coveredpos, _, blockstarts, _, tendpos = process_cigar([1, 1, 0] + [1] * 7, [(pysam.CMATCH, 10)], 5, None, None)
    assert coveredpos == [0] * 5 + [1, 1, 0] + [1] * 7
    assert (blockstarts, tendpos) == ([5], 15)


def test_exact_match_has_no_divergence():
    assert check_splicesites([1] * 200, EXONS, 0, 200, 't') == 0


def test_indel_moved_away_from_the_splice_site_still_ranks():
    # aligned to a transcript whose splice site is 15 bases off the read's, the read
    # has a 15-base insertion, here put 14 bases before the splice site: it passes
    # the splice site check, but ranks below an exact match
    coveredpos = [1] * 200
    coveredpos[85] += 15
    assert check_splicesites(coveredpos, EXONS, 0, 200, 't') == 15


def test_indel_at_the_splice_site_fails():
    coveredpos = [1] * 200
    coveredpos[99] += 15
    assert check_splicesites(coveredpos, EXONS, 0, 200, 't') is None


def test_unaligned_positions_are_not_divergence():
    # a read starting 10 bases before the splice site
    assert check_splicesites([0] * 90 + [1] * 110, EXONS, 90, 200, 't') == 0


def test_read_ending_just_past_the_splice_site_fails():
    # the alignment ends a base into the second exon, as one of the read's clipped
    # bases can match by chance; the window's other 4 bases there are unmatched
    assert check_splicesites([1] * 101, EXONS, 0, 101, 't') is None


def test_read_starting_just_before_the_splice_site_fails():
    assert check_splicesites([0] * 99 + [1] * 101, EXONS, 99, 200, 't') is None

def _passing(tname, left_clipping, right_clipping):
    # an entry of get_best_transcript's passing transcripts: ranking keys, then the name
    # and each end's (intron index, distance to it, distance to the transcript end, clipping)
    return [0, -1000, -900, -2, left_clipping + right_clipping, 1500, tname,
            (0, 100, 5, left_clipping), (1, 100, 5, right_clipping)]


def test_no_extra_clipping_rejects_alignments_clipped_at_a_transcript_end():
    from flair.count_sam_transcripts import return_best_transcript_stringent
    clipped = [_passing('t1', 0, 30)]
    assert return_best_transcript_stringent(clipped, [0, 0], 50, 'r')[0][0] == 't1'
    assert return_best_transcript_stringent(clipped, [0, 0], 50, 'r', no_extra_clipping=True) is None
    both = [_passing('t1', 0, 30), _passing('t2', 0, 0)]
    assert return_best_transcript_stringent(both, [0, 0], 50, 'r', no_extra_clipping=True)[0][0] == 't2'
