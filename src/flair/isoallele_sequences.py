"""An isoform's sequence with variants applied, its coding sequence and protein, and
how the protein differs from the isoform's own, for flair isoalleles.
"""
import difflib
import hashlib
from flair.isoform_data import COMPBASE, translate_codon
from flair.predictProductivity import STOP_CODON_SEQS, calc_ptc


def variant_label(key):
    "a variant (chrom, pos, ref, alt), pos 1-based as in a VCF, as chrom:pos:ref>alt"
    return f'{key[0]}:{key[1]}:{key[2]}>{key[3]}'


def exon_index(exons, pos, ref):
    """the index of the exon, of (start, end) 0-based half open, that wholly holds a
    variant at 1-based pos with this ref allele, or None"""
    return next((j for j, (s, e) in enumerate(exons) if s < pos and pos - 1 + len(ref) <= e), None)


def build_transcript(genome, iso, variants):
    """The isoform's sequence, 5' to 3', with the variants applied, as (sequence,
    origins, exon sizes, variants not applied).

    origins gives each base's 0-based genomic position, or None for a base an
    insertion adds, so a position in the sequence can be taken back to the genome, for
    the CDS ends.  The exon sizes, 5' to 3', include the variants' indels, for the
    premature stop rule.  A variant isn't applied if it isn't wholly in one exon,
    overlaps another, or has a ref allele the genome disagrees with."""
    exons = sorted((b.start, b.end) for b in iso.blocks)
    pieces = [[(base, s + k) for k, base in enumerate(genome.fetch(iso.chrom, s, e).upper())] for s, e in exons]
    skipped, taken = [], []
    for key in sorted(variants, key=lambda k: -k[1]):  # from the right, so earlier positions don't shift
        _, pos, ref, alt = key
        i = exon_index(exons, pos, ref)
        if i is None or any(t[1] < pos + len(ref) and pos < t[1] + len(t[2]) for t in taken):
            skipped.append(key)
            continue
        offset = pos - 1 - exons[i][0]
        if ''.join(b for b, _ in pieces[i][offset:offset + len(ref)]) != ref.upper():
            skipped.append(key)
            continue
        # the alt's bases take the ref's positions where it has them; the rest of a longer alt is inserted
        new = [(b, pos - 1 + k if k < len(ref) else None) for k, b in enumerate(alt.upper())]
        pieces[i] = pieces[i][:offset] + new + pieces[i][offset + len(ref):]
        taken.append(key)
    bases = [x for piece in pieces for x in piece]
    sizes = [len(piece) for piece in pieces]
    if iso.strand == '-':
        bases = [(COMPBASE.get(b, 'N'), p) for b, p in reversed(bases)]
        sizes = sizes[::-1]
    return ''.join(b for b, _ in bases), [p for _, p in bases], sizes, skipped


def call_cds(iso, seq, origins, sizes, nmd_exempt):
    """The coding sequence of a variant transcript, from the isoform's annotated start
    codon, the BED's thickStart and thickEnd, found in the variant sequence by its
    genomic position and read to the first in-frame stop codon, so a variant that
    changes the stop (nonsense, frameshift, stop lost) moves the CDS end.  Productivity,
    as flair transcriptome calls it:
      PRO  a stop codon, not a premature one
      PTC  a stop codon more than 55 bases before the last exon junction (calc_ptc),
           unless nmd_exempt: the annotation calls the unchanged isoform productive
      NST  no stop codon before the transcript's end
      NGO  no ORF: a non-coding isoform, or a start codon a variant deleted or changed
    Returns (productivity, protein or None, (thickStart, thickEnd), note).  The protein,
    without its stop, is given for PRO, PTC and NST: a truncated protein is also what
    the transcript encodes."""
    none = (iso.chromStart, iso.chromStart)
    if iso.thickStart is None or iso.thickEnd is None or iso.thickStart >= iso.thickEnd or iso.productivity == 'NGO':
        return 'NGO', None, none, 'non-coding isoform'
    start_g = iso.thickStart if iso.strand == '+' else iso.thickEnd - 1  # the start codon's first base
    if start_g not in origins:
        return 'NGO', None, none, 'start codon deleted'
    i = origins.index(start_g)
    if seq[i:i + 3] != 'ATG':
        return 'NGO', None, none, f'start codon {seq[i:i + 3]}, not ATG'
    aas, stop = [], None
    for j in range(i, len(seq) - 2, 3):
        codon = seq[j:j + 3]
        if codon in STOP_CODON_SEQS:
            stop = j
            break
        aas.append(translate_codon(codon))
    last = (stop + 2) if stop is not None else (i + 3 * len(aas) - 1)
    end_g = next(origins[k] for k in range(last, -1, -1) if origins[k] is not None)  # an inserted base has none
    thick = (start_g, end_g + 1) if iso.strand == '+' else (end_g, start_g + 1)
    if stop is None:
        return 'NST', ''.join(aas), thick, 'no stop codon'
    productivity = calc_ptc(sizes, stop + 3, 'isoform', {'isoform': nmd_exempt})
    return productivity, ''.join(aas), thick, ''


def _cds_bases_before(iso, key):
    "the isoform's CDS bases 5' of a variant, in transcript order"
    exons = [(b.start, b.end) for b in iso.blocks]
    if iso.strand == '+':
        lo, hi = iso.thickStart, key[1] - 1
    else:
        lo, hi = key[1] - 1 + len(key[2]), iso.thickEnd
    return sum(max(0, min(e, hi) - max(s, lo)) for s, e in exons)


def frameshift_residue(iso, variants):
    """The residue, 0-based, of the isoform's own protein at which the variants shift
    its frame: of the variants in its CDS, taken 5' to 3', the first after which their
    indels add up to a non-multiple of 3.  None if they don't shift it."""
    if iso.thickStart is None or iso.thickEnd is None or iso.thickStart >= iso.thickEnd:
        return None
    in_cds = sorted((k for k in variants if iso.thickStart < k[1] <= iso.thickEnd),
                    key=lambda k: k[1], reverse=iso.strand == '-')
    shift = 0
    for key in in_cds:
        shift += len(key[3]) - len(key[2])
        if shift % 3:
            return _cds_bases_before(iso, key) // 3
    return None


def _aligned(ref, alt, frameshift_at):
    """The alignment of two proteins, as difflib opcodes.  Proteins of one length,
    without a frameshift, differ only by substitutions, compared residue by residue.
    Otherwise their shared start is left out of the alignment, which difflib places
    badly in a run of one residue, as a substitution in poly-alanine taken for an
    insertion."""
    if len(ref) == len(alt) and frameshift_at is None:
        return [('replace', i, i + 1, i, i + 1) for i in range(len(ref)) if ref[i] != alt[i]]
    shared = next((i for i, (a, b) in enumerate(zip(ref, alt)) if a != b), min(len(ref), len(alt)))
    ops = difflib.SequenceMatcher(None, ref[shared:], alt[shared:], autojunk=False).get_opcodes()
    return [(op, i1 + shared, i2 + shared, j1 + shared, j2 + shared) for op, i1, i2, j1, j2 in ops]


def protein_change(ref, alt, frameshift_at=None):  # noqa: C901 - a case per kind of HGVS change
    """The changes from a reference protein to a variant one, in HGVS-like terms,
    ';' joined: p.G12A (missense), p.R196* (a new stop), p.E746_A750del,
    p.N771_H773dup, p.V600_K601insX, p.X10_Y12delinsZ, p.*1211ext (the stop is lost),
    p.L858fs (a frameshift: the first residue changed from frameshift_at, the residue
    frameshift_residue gives, on; nothing after it is described), or 'none'.  The proteins are aligned (difflib), so separate changes, such as a
    germline missense and a somatic duplication, are described separately.  Deletions,
    insertions and duplications are placed as far toward the C terminus as they go, as
    HGVS does.  '' for no reference protein, 'no protein' for no variant one."""
    if ref is None:
        return ''
    if alt is None:
        return 'no protein'
    if ref == alt:
        return 'none'

    def span(a, b):
        return f'{ref[a]}{a + 1}' + (f'_{ref[b]}{b + 1}' if b > a else '')

    changes = []
    for op, i1, i2, j1, j2 in _aligned(ref, alt, frameshift_at):
        if op == 'equal':
            continue
        if frameshift_at is not None and i2 > frameshift_at:
            k = next((i for i in range(max(i1, frameshift_at), len(ref)) if i >= len(alt) or ref[i] != alt[i]),
                     len(ref))
            changes.append(f'p.{ref[k]}{k + 1}fs' if k < len(ref) else f'p.*{len(ref) + 1}ext')
            break
        if op == 'replace' and i2 - i1 == j2 - j1:  # substitutions
            changes += [f'p.{ref[i]}{i + 1}{alt[j]}' for i, j in zip(range(i1, i2), range(j1, j2)) if ref[i] != alt[j]]
            continue
        if op == 'delete' and i2 == len(ref) and j1 == len(alt):  # the variant protein stops early
            changes.append(f'p.{ref[i1]}{i1 + 1}*')
        elif op == 'insert' and i1 == len(ref):  # it reads past the reference stop
            changes.append(f'p.*{len(ref) + 1}ext')
        elif op == 'delete':
            while i2 < len(ref) and ref[i1] == ref[i2]:  # shift right
                i1, i2 = i1 + 1, i2 + 1
            changes.append(f'p.{span(i1, i2 - 1)}del')
        elif op == 'insert':
            ins = alt[j1:j2]
            while i1 < len(ref) and ref[i1] == ins[0]:
                ins, i1 = ins[1:] + ref[i1], i1 + 1
            n = len(ins)
            if i1 >= n and ref[i1 - n:i1] == ins:
                changes.append(f'p.{span(i1 - n, i1 - 1)}dup')
            elif 0 < i1 < len(ref):
                changes.append(f'p.{ref[i1 - 1]}{i1}_{ref[i1]}{i1 + 1}ins{ins}')
            else:
                changes.append(f'p.ins{ins}')
        else:
            changes.append(f'p.{span(i1, i2 - 1)}delins{alt[j1:j2]}')
    return ';'.join(changes) or 'none'


def protein_id(aaseq, reference_ids):
    """A protein's id: the isoform set's own (its aaseq table) for a protein it already
    has, else FLPV and the first 10 hex digits of the protein's SHA-1, the same for the
    same protein in any run"""
    if aaseq in reference_ids:
        return reference_ids[aaseq]
    return 'FLPV' + hashlib.sha1(aaseq.encode()).hexdigest()[:10]
