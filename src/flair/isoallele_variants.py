"""The variants flair isoalleles puts in isoform sequences, and each read's allele and
haplotype at them.

Variants are keys (chrom, pos, ref, alt), pos 1-based as in a VCF.  A VCF record
passes with FILTER PASS or '.', as flair variantquant takes them.

Germline variants are the calls of the sample, or of the normal of a pair, phased
by the VCF that haplotagged the reads (WhatsHap, DeepVariant, HiPhase and the like),
which can hold them all or, when only some were phased, as only the hets, just
those:
  homozygous     in every sequence
  phased het     on the haplotype its phased genotype gives, H1 for the first allele
  unphased het   placed per isoform from its reads, or left out (place_unphased)

Somatic variants, of a pair only, are one of:
  tumor_vcf      the tumor's passing calls the phased VCF does not call, with at least
                 min_normal_depth normal reads informing the site; the rest are
                 excluded, as the normal can't show they aren't germline
  somatic_vcf    every passing call of a somatic caller's VCF, as it calls them
  --variant      only the variants given, whatever the VCFs say
"""
import re
from collections import namedtuple
import pysam
from flair import FlairInputDataError

Variant = namedtuple('Variant', ('chrom', 'pos', 'ref', 'alt'))

GermlineVariant = namedtuple('GermlineVariant', ('key', 'kind', 'haplotype', 'phase_set'))
HOMOZYGOUS, PHASED_HET, UNPHASED_HET = 'homozygous', 'phased het', 'unphased het'

# a read's haplotype without an HP tag
UNTAGGED = 'untagged'

# an indel read within this many bases of a called indel, of the same kind and
# length, carries it: long reads don't always place an indel where the caller did
INDEL_WINDOW = 5


def passes(rec):
    "whether a VCF record's FILTER is PASS or '.'"
    filters = list(rec.filter.keys())
    return not filters or filters == ['PASS']


def _genotype(rec):
    gt = rec.samples[0].get('GT') if len(rec.samples) else None
    return tuple(gt) if gt else ()


def called_alts(rec):
    """the alts of a passing record its first sample's genotype holds; every alt for a
    VCF of no samples or genotypes"""
    if not passes(rec):
        return []
    gt = _genotype(rec)
    if not gt or None in gt:
        return list(rec.alts or ()) if not gt else []
    return [alt for i, alt in enumerate(rec.alts or ()) if i + 1 in gt]


def open_vcf(path):
    "a VCF to fetch regions of: bgzipped with a tabix or csi index"
    try:
        return pysam.VariantFile(path)
    except (OSError, ValueError) as ex:
        raise FlairInputDataError(f'can\'t read VCF {path}: {ex}') from ex


def fetch(vcf, path, chrom, start, end):
    "the records of a region, none for a chromosome the VCF doesn't have"
    if chrom not in vcf.header.contigs:
        return []
    try:
        return list(vcf.fetch(chrom, start, end))
    except ValueError as ex:
        raise FlairInputDataError(f'{path} must be bgzipped and indexed (tabix -p vcf) to be read by region: '
                                  f'{ex}') from ex


def _phases(records):
    "(pos, ref, alt) -> (haplotype, phase set) of the phased genotypes of records"
    phases = {}
    for rec in records:
        gt, sample = _genotype(rec), rec.samples[0] if len(rec.samples) else None
        if sample is None or not sample.phased or not gt or None in gt:
            continue
        for i, alt in enumerate(rec.alts or ()):
            if i + 1 in gt and len(set(gt)) > 1:
                phases[(rec.pos, rec.ref, alt)] = (f'H{gt.index(i + 1) + 1}', sample.get('PS'))
    return phases


def germline_variants(records, phased_records=None):
    """The germline variants of the calls of a region, records, as GermlineVariant,
    phased by phased_records, the phased VCF's records of the region, or by their own
    genotypes when the calls are the phased VCF's.  A het genotype of two alts,
    phased, puts one on each haplotype."""
    records = list(records)
    phases = _phases(records if phased_records is None else phased_records)
    out = []
    for rec in records:
        gt = _genotype(rec)
        if not passes(rec) or not gt or None in gt:
            continue
        for i, alt in enumerate(rec.alts or ()):
            if i + 1 not in gt:
                continue
            key = Variant(rec.chrom, rec.pos, rec.ref, alt)
            phase = phases.get((rec.pos, rec.ref, alt))
            if len(set(gt)) == 1:
                out.append(GermlineVariant(key, HOMOZYGOUS, None, None))
            elif phase:
                out.append(GermlineVariant(key, PHASED_HET, *phase))
            else:
                out.append(GermlineVariant(key, UNPHASED_HET, None, None))
    return out


def called_variants(records):
    "the variants passing records call"
    return {Variant(rec.chrom, rec.pos, rec.ref, alt) for rec in records for alt in called_alts(rec)}


_VARIANT_RE = re.compile(r'(\S+):(\d+)\s*(?::?([ACGTNacgtn]+)>([ACGTNacgtn]+))?')


def parse_variant(text):
    """A --variant: 'chr7:140753336' -> ('chr7', 140753336, None, None), and
    'chr7:140753336A>T' or 'chr7:140753336:A>T' -> ('chr7', 140753336, 'A', 'T')"""
    m = _VARIANT_RE.fullmatch(text.strip())
    if not m:
        raise FlairInputDataError(f'--variant {text}: expected CHROM:POS or CHROM:POS:REF>ALT')
    chrom, pos, ref, alt = m.groups()
    return chrom, int(pos), ref and ref.upper(), alt and alt.upper()


def variant_keys(texts, vcf_paths):
    """The --variant variants as Variant keys.  CHROM:POS:REF>ALT is taken as given;
    CHROM:POS takes every alt that a record of one of the VCFs has at that position,
    and is an error if none has."""
    keys = []
    vcfs = [(path, open_vcf(path)) for path in vcf_paths if path]
    for text in texts:
        chrom, pos, ref, alt = parse_variant(text)
        if ref:
            keys.append(Variant(chrom, pos, ref, alt))
            continue
        found = {Variant(rec.chrom, rec.pos, rec.ref, a) for path, vcf in vcfs
                 for rec in fetch(vcf, path, chrom, pos - 1, pos) if rec.pos == pos for a in rec.alts or ()}
        if not found:
            raise FlairInputDataError(f'--variant {text}: no VCF has a record at {chrom}:{pos}; give it as '
                                      'CHROM:POS:REF>ALT')
        keys.extend(sorted(found))
    for vcf in vcfs:
        vcf[1].close()
    return sorted(set(keys))


def haplotype(read):
    "'H1', 'H2' and so on from the read's HP tag, or UNTAGGED if the phasing couldn't place it"
    return f'H{read.get_tag("HP")}' if read.has_tag('HP') else UNTAGGED


class ReadAlleles:
    """A read's allele at variants: 'alt', 'ref' or None, when the read doesn't inform
    the variant (it doesn't reach it, has an intron or deletion over it, or a third
    base there).  For an SNV it's the read's base; for an indel, whether the read has
    an insertion or deletion of the same length within INDEL_WINDOW of it, counted
    only for a read aligned at least INDEL_WINDOW bases past the site both sides.
    The read's aligned bases are looked up once, however many variants."""

    def __init__(self, read):
        self.read = read
        self._ref_to_query = None

    def _bases(self):
        if self._ref_to_query is None:
            self._ref_to_query = {rpos: qpos for qpos, rpos in self.read.get_aligned_pairs(matches_only=True)}
        return self._ref_to_query

    def allele(self, key):
        read = self.read
        _, pos, ref, alt = key
        if len(ref) == len(alt) == 1:
            qpos = self._bases().get(pos - 1)
            if qpos is None:
                return None
            base = read.query_sequence[qpos]
            return 'alt' if base == alt else 'ref' if base == ref else None
        kind, length = (1, len(alt) - len(ref)) if len(alt) > len(ref) else (2, len(ref) - len(alt))
        if read.reference_start > pos - 1 - INDEL_WINDOW or read.reference_end < pos + len(ref) + INDEL_WINDOW:
            return None
        rpos, covered = read.reference_start, False
        for op, oplen in read.cigartuples:
            if op in (0, 7, 8) and rpos <= pos < rpos + oplen:
                covered = True
            if op == kind and oplen == length and abs(rpos - pos) <= INDEL_WINDOW:
                return 'alt'
            if op in (0, 2, 3, 7, 8):
                rpos += oplen
        return 'ref' if covered else None
