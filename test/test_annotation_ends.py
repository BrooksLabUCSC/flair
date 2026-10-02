from flair.annotation_data import AnnotData, _save_transcript_ends


def _annots(*transcripts):
    annots = AnnotData()
    for strand, start, end, tags in transcripts:
        _save_transcript_ends(annots, strand, start, end, tags)
    for strand in ('+', '-'):
        annots.transcript_ends[strand].sort()
    return annots


def test_matches_both_ends_within_window():
    annots = _annots(('+', 1000, 5000, []))
    assert annots.has_transcript_ends('+', 1050, 4950, 100)
    assert annots.has_transcript_ends('+', 900, 5100, 100)
    assert not annots.has_transcript_ends('+', 1150, 5000, 100)   # start too far
    assert not annots.has_transcript_ends('+', 1000, 4800, 100)   # end too far
    assert not annots.has_transcript_ends('-', 1000, 5000, 100)   # other strand


def test_both_ends_from_one_transcript():
    # the start of one transcript and the end of another don't match
    annots = _annots(('+', 1000, 3000, []), ('+', 2000, 5000, []))
    assert not annots.has_transcript_ends('+', 1000, 5000, 100)
    assert annots.has_transcript_ends('+', 2000, 5000, 100)


def test_unconfirmed_ends_excluded():
    annots = _annots(('+', 1000, 5000, ['mRNA_end_NF']), ('-', 1000, 5000, ['mRNA_start_NF', 'basic']),
                     ('+', 7000, 9000, ['cds_start_NF', 'basic']))
    assert not annots.has_transcript_ends('+', 1000, 5000, 100)
    assert not annots.has_transcript_ends('-', 1000, 5000, 100)
    # an incomplete CDS doesn't make the transcript's ends unconfirmed
    assert annots.has_transcript_ends('+', 7000, 9000, 100)


def test_same_start_sorts():
    annots = _annots(('+', 1000, 5000, []), ('+', 1000, 3000, []))
    assert annots.has_transcript_ends('+', 1000, 3050, 100)
