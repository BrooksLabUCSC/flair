from flair.count_sam_transcripts import check_splicesites

# a transcript of two 100-base exons, its splice site at offset 100; coveredpos has
# 1 for a matched position, 0 for unmatched, and 1 + insert size before an insertion
EXONS = [100, 100]


def test_exact_match_has_no_divergence():
    assert check_splicesites([1] * 200, EXONS, 1, 201, 't') == 0


def test_indel_moved_away_from_the_splice_site_still_ranks():
    # aligned to a transcript whose splice site is 15 bases off the read's, the read
    # has a 15-base insertion, here put 14 bases before the splice site: it passes
    # the splice site check, but ranks below an exact match
    coveredpos = [1] * 200
    coveredpos[85] += 15
    assert check_splicesites(coveredpos, EXONS, 1, 201, 't') == 15


def test_indel_at_the_splice_site_fails():
    coveredpos = [1] * 200
    coveredpos[99] += 15
    assert check_splicesites(coveredpos, EXONS, 1, 201, 't') is None


def test_unaligned_positions_are_not_divergence():
    # a read starting 10 bases before the splice site
    assert check_splicesites([0] * 90 + [1] * 110, EXONS, 91, 201, 't') == 0
