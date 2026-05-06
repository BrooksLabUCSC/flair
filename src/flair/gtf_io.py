"""
Simplistic, non-validating GTF parser.
"""
import re
from enum import Enum
from typing import Optional
from collections import defaultdict
from flair.pycbio.sys import fileOps
from flair import SeqRange
from flair.interval_index import IntervalIndex

StrNone = Optional[str]
StrSetNone = Optional[set[str]]
Attrs = dict[str, str]


##
# Feature names. Use sets for filters
##
TRANSCRIPT_FEATURES = frozenset([
    'transcript', 'mRNA', 'lncRNA', 'miRNA',
    'ncRNA', 'rRNA', 'tRNA', 'snRNA', 'snoRNA',
    'processed_transcript', 'pseudogenic_transcript'
])

EXON_FEATURE = "exon"
EXON_FEATURES = frozenset([EXON_FEATURE])

CDS_FEATURE = "CDS"
CDS_FEATURES = frozenset([CDS_FEATURE])

TRANSCRIPT_EXON_FEATURES = TRANSCRIPT_FEATURES | EXON_FEATURES

FLAIR_ATTRS = ('gene_id', 'gene_name', 'transcript_id')

_ALL_ATTR_RE = re.compile(r'(\w+)\s+(?:"([^"]*)"|([^;\s]+))')

class GtfAttrsSet(Enum):
    """Selects which GTF attributes to parse.
    ALL parses every attribute; FLAIR parses only the attributes in FLAIR_ATTRS,
    which is faster for flair's use."""
    ALL = 'all'
    FLAIR = 'flair'


class GtfParseError(Exception):
    """Error parsing GTF record."""
    pass

class GtfIdError(KeyError):
    """GTF id lookup error"""
    pass

class GtfRecord:
    """Base GTF record.  Text columns of '.' are set to None.  Repeated attributes
    are converted to list of values.  If attrs is specified, ownership is passed
    to this object.
    """
    def __init__(self, chrom: str, source: str, feature: str,
                 start: int, end: int, score: str, strand: str, frame: str, *,
                 attrs: Attrs = None,
                 gene_id: str = None, gene_name: str = None, transcript_id: str = None,
                 exon_number: int = None):
        self.chrom = chrom
        self.source = source
        self.feature = feature
        self.start = start  # 0-based
        self.end = end
        self.score = score
        self.strand = strand
        self.frame = frame

        # Make a copy of attrs and add optional parameters
        self.attrs = attrs if attrs is not None else {}
        if gene_id is not None:
            if re.search(r'\s', gene_id):
                raise GtfParseError(f"white space not allowed  in gene_id: {gene_id}")
            self.attrs['gene_id'] = gene_id
        if gene_name is not None:
            self.attrs['gene_name'] = gene_name
        if transcript_id is not None:
            if re.search(r'\s', transcript_id):
                raise GtfParseError(f"white space not allowed  in transcript_id: {transcript_id}")
            self.attrs['transcript_id'] = transcript_id
        if exon_number is not None:
            self.attrs['exon_number'] = exon_number

    @property
    def gene_id(self):
        return self.attrs.get("gene_id")

    @gene_id.setter
    def gene_id(self, value):
        self.attrs['gene_id'] = value

    @property
    def gene_name(self):
        return self.attrs.get("gene_name")

    @property
    def transcript_id(self):
        return self.attrs.get("transcript_id")

    @property
    def exon_number(self):
        return self.attrs.get("exon_number")

    def __len__(self) -> int:
        return self.end - self.start

    @property
    def coords(self) -> SeqRange:
        """coordinates of sequence"""
        return SeqRange(self.chrom, self.start, self.end, self.strand)

    @property
    def coords_no_strand(self) -> SeqRange:
        """coordinates of sequence, without strand"""
        return SeqRange(self.chrom, self.start, self.end)

    def __str__(self) -> str:
        """Return GTF format line."""
        return _format_gtf_columns(self.chrom, self.source, self.feature, self.start, self.end,
                                   self.score, self.strand, self.frame, self.attrs)

def gtf_record_sort_key(rec):
    return (rec.chrom, rec.start, rec.end)

class GtfExon(GtfRecord):
    """GTF exon."""

    def __init__(self, chrom: str, source: str, feature: str, start: int, end: int,
                 score: str, strand: str, frame: str, *,
                 attrs: Attrs = None,
                 gene_id: str = None, gene_name: str = None, transcript_id: str = None,
                 exon_number: int = None):
        super().__init__(chrom, source, feature, start, end, score, strand, frame,
                         attrs=attrs, gene_id=gene_id, gene_name=gene_name, transcript_id=transcript_id,
                         exon_number=exon_number)

class GtfCDS(GtfRecord):
    """GTF CDS (coding sequence)."""

    def __init__(self, chrom: str, source: str, feature: str, start: int, end: int,
                 score: str, strand: str, frame: str, *,
                 attrs: Attrs = None,
                 gene_id: str = None, gene_name: str = None, transcript_id: str = None,
                 exon_number: int = None):
        super().__init__(chrom, source, feature, start, end, score, strand, frame,
                         attrs=attrs, gene_id=gene_id, gene_name=gene_name, transcript_id=transcript_id,
                         exon_number=exon_number)

class GtfTranscript(GtfRecord):
    """GTF transcript with exons."""
    def __init__(self, chrom: str, source: str, feature: str, start: int, end: int,
                 score: str, strand: str, frame: str, *,
                 attrs: Attrs = None,
                 gene_id: str = None, gene_name: str = None, transcript_id: str = None,
                 exon_number: int = None):
        super().__init__(chrom, source, feature, start, end, score, strand, frame,
                         attrs=attrs, gene_id=gene_id, gene_name=gene_name, transcript_id=transcript_id,
                         exon_number=exon_number)

        self.exons: list[GtfExon] = []
        self.cds_recs: list[GtfCDS] = []

    def add_exon(self, exon: GtfExon) -> None:
        """Add exon to transcript."""
        self.exons.append(exon)

    def add_cds(self, cds: GtfCDS) -> None:
        """Add CDS to transcript."""
        self.cds_recs.append(cds)

    def sort_children(self) -> None:
        self.exons.sort(key=gtf_record_sort_key)
        self.cds_recs.sort(key=gtf_record_sort_key)

class GtfData:
    """Data from a GTF file.  This is a container of GtfTranscript objects
    and there children.  There is no gene objects, as they are not required records
    in GTF."""
    def __init__(self, gtf_file=None):
        self.gtf_file = gtf_file  # saved for error messages
        self.transcripts = []
        self.transcripts_by_id: dict[str, GtfTranscript] = {}
        # transcripts by chrom then range overlap.
        self.transcripts_by_range = defaultdict(IntervalIndex)

    def add_transcript(self, transcript: GtfTranscript):
        if transcript.transcript_id in self.transcripts_by_id:
            raise GtfParseError(f"adding duplicate transcript id: `{transcript.transcript_id}'")
        self.transcripts.append(transcript)
        self.transcripts_by_id[transcript.transcript_id] = transcript
        self.transcripts_by_range[transcript.chrom].add(transcript.start, transcript.end, transcript)

    def get_transcript(self, transcript_id):
        """return transcript for id or None if not found"""
        return self.transcripts_by_id.get(transcript_id)

    def fetch_transcript(self, transcript_id):
        """return transcript for id or error if not found"""
        try:
            return self.transcripts_by_id[transcript_id]
        except KeyError:
            raise GtfIdError(f"unknown transcript id `{transcript_id}'")

    def iter_transcript_ids(self):
        return self.transcripts_by_id.keys()

    def get_chroms(self):
        """get a sorted list of the chrom names"""
        return sorted(self.transcripts_by_range.keys())

    def iter_overlap_transcripts(self, chrom, start, end, *, strand=None):
        """Generator overlapping transcripts, optionally filtering for strand"""
        # defaultdict will handle chrom not in GTF
        for transcript in self.transcripts_by_range[chrom].overlap(start, end):
            if (strand is None) or (transcript.strand == strand):
                yield transcript

    def iter_overlap_transcripts_sr(self, seq_range):
        """Generator overlapping transcripts given a SeqRange object,
        optionally filtering for strand."""
        yield from self.iter_overlap_transcripts(seq_range.name, seq_range.start, seq_range.end,
                                                 strand=seq_range.strand)

    def subset_for_region(self, chrom, start, end):
        """Return a new GtfData with transcripts overlapping [start, end) on chrom.
        Transcript objects are shared, not copied."""
        sub = GtfData(self.gtf_file)
        for transcript in self.iter_overlap_transcripts(chrom, start, end):
            sub.add_transcript(transcript)
        return sub


def _parse_attribute_match(match: re.Match) -> tuple[str, str | int | float]:
    """Parse a single attribute match into key-value pair, converting
    unquote values to int or float if possible"""

    key, quoted_value, unquoted_value = match.groups()
    if quoted_value is not None:
        return key, quoted_value
    else:
        # It was unquoted - try to convert to number
        try:
            # Try int first
            if '.' not in unquoted_value:
                return key, int(unquoted_value)
            else:
                return key, float(unquoted_value)
        except ValueError:
            # If conversion fails, keep as string
            return key, unquoted_value

def _parse_all_attributes(attrs_str: str, attr_re=_ALL_ATTR_RE) -> Attrs:
    """Parse GTF attributes string into dict."""
    attrs = {}
    for attr_str in attr_re.finditer(attrs_str):
        key, value = _parse_attribute_match(attr_str)
        if key in attrs:
            if not isinstance(attrs[key], list):
                attrs[key] = [attrs[key]]  # convert to list
            attrs[key].append(value)
        else:
            attrs[key] = value

    return attrs

def _find_flair_attr_value(attrs_str: str, key: str):
    """Find the quoted value of a single GTF attribute by key using str.find().
    Returns the value as a str or None if not found."""
    alen = len(attrs_str)
    idx = attrs_str.find(key)
    if idx < 0:
        return None
    # skip whitespace to opening quote
    val_idx = idx + len(key)
    while val_idx < alen and attrs_str[val_idx] == ' ':
        val_idx += 1
    if val_idx >= alen or attrs_str[val_idx] != '"':
        return None
    val_start = val_idx + 1
    val_end = attrs_str.find('"', val_start)
    return attrs_str[val_start:val_end] if val_end >= 0 else None

def _parse_flair_attributes(attrs_str: str) -> Attrs:
    """Fast-path parser for FLAIR_ATTRS using str.find() instead of regex."""
    attrs = {}
    for key in FLAIR_ATTRS:
        value = _find_flair_attr_value(attrs_str, key)
        if value is not None:
            attrs[key] = value
    return attrs


_ATTRS_SET_TO_PARSER = {
    GtfAttrsSet.ALL: _parse_all_attributes,
    GtfAttrsSet.FLAIR: _parse_flair_attributes,
}


def _parse_coordinates(start_str: str, end_str: str) -> tuple[int, int]:
    try:
        start = int(start_str) - 1  # Convert to 0-based
        end = int(end_str)
    except ValueError as exc:
        raise GtfParseError(f"Invalid coordinates `{start_str}' `{end_str}'") from exc
    if start >= end:
        raise GtfParseError(f"Coordinates `{start_str}' > `{end_str}'")
    return start, end

def _parse_strand(strand: str) -> str:
    if strand not in {'+', '-', '.'}:
        raise GtfParseError(f"Invalid strand '{strand}'")
    return strand

def _parse_score(score_str: str) -> str:
    try:
        if score_str == '.':
            return None
        else:
            return float(score_str)
    except Exception as exc:
        raise GtfParseError(f"Invalid score '{score_str}'") from exc

def _parse_frame(frame_str: str) -> str:
    try:
        if frame_str == '.':
            return None
        else:
            frame = int(frame_str)
    except Exception as exc:
        raise GtfParseError(f"Invalid frame '{frame_str}'") from exc

    if not (0 <= frame <= 2):
        raise GtfParseError(f"Frame must be in the range 0..2 or `.', got'{frame}'")
    return frame


def _gtf_record_class(feature):
    if feature in TRANSCRIPT_FEATURES:
        return GtfTranscript
    elif feature == EXON_FEATURE:
        return GtfExon
    elif feature == CDS_FEATURE:
        return GtfCDS
    else:
        return GtfRecord

def _check_id_whitespace(attrs):
    for key in ('gene_id', 'transcript_id'):
        val = attrs.get(key)
        if val is not None and re.search(r'\s', val):
            raise GtfParseError(f"white space not allowed in {key}: {val!r}")


def _parse_gtf_line(line: str, include_features: StrSetNone, attrs_parser=_parse_flair_attributes) -> GtfRecord:
    """Parse a single GTF line into a GtfRecord or derived class."""
    # skip empty and comments
    line = line.rstrip()
    if (len(line) == 0) or line.startswith('#'):
        return None

    fields = line.split('\t')
    if len(fields) != 9:
        raise GtfParseError(f"Expected 9 fields, got {len(fields)}")
    if (include_features is not None) and (fields[2] not in include_features):
        return None

    start, end = _parse_coordinates(fields[3], fields[4])
    attrs = attrs_parser(fields[8])
    _check_id_whitespace(attrs)

    cls = _gtf_record_class(fields[2])
    return cls(chrom=fields[0],
               source=fields[1],
               feature=fields[2],
               start=start,
               end=end,
               score=_parse_score(fields[5]),
               strand=_parse_strand(fields[6]),
               frame=_parse_frame(fields[7]),
               attrs=attrs)

def _format_attr(key, value):
    if isinstance(value, (int, float)):
        return f'{key} {value};'
    else:
        return f'{key} "{value}";'

def _format_attr_key_val(key, value):
    attr_strs = []
    if isinstance(value, list):
        for value_n in value:
            attr_strs.append(_format_attr(key, value_n))
    else:
        attr_strs = [_format_attr(key, value)]
    return attr_strs

def _format_attrs(attrs):
    """Build list of key=val, handling duplicate keys"""
    attrs_strs = []
    for key, value in attrs.items():
        attrs_strs += _format_attr_key_val(key, value)
    return ' '.join(attrs_strs)

def _str_or_dot(val):
    return str(val) if val is not None else '.'

def _format_gtf_columns(chrom, source, feature, start, end, score, strand, frame,
                        attrs) -> str:
    """Return GTF format line."""

    # Join all fields with tabs
    fields = [chrom, source, feature, str(start + 1), str(end),
              _str_or_dot(score), _str_or_dot(strand), _str_or_dot(frame),
              _format_attrs(attrs)]
    return '\t'.join(fields)

def _gtf_record_iter(gtf_file: str, include_features, attrs_parser):
    """Internal: iterate parsed GTF records using an attribute parser callable."""
    with fileOps.opengz(gtf_file) as fh:
        for line_num, line in enumerate(fh, start=1):
            try:
                rec = _parse_gtf_line(line, include_features, attrs_parser)
                if rec is not None:
                    yield rec
            except GtfParseError as exc:
                raise GtfParseError(f"{gtf_file}:{line_num}: invalid GTF record") from exc

def gtf_record_parser(gtf_file: str, *, include_features: StrSetNone = None, attrs: GtfAttrsSet = GtfAttrsSet.FLAIR):
    """Parse GTF file, yields GtfRecord GtfTranscript, GtfExon, or GtfCDS objects.
    File maybe compressed"""
    yield from _gtf_record_iter(gtf_file, include_features, _ATTRS_SET_TO_PARSER[attrs])

def _load_gtf_records(gtf_file, gtf_data, transcript_id_to_exons, transcript_id_to_cds_recs, include_features, attrs_parser):
    for rec in _gtf_record_iter(gtf_file, include_features, attrs_parser):
        if isinstance(rec, GtfTranscript):
            gtf_data.add_transcript(rec)
        elif isinstance(rec, GtfExon):
            transcript_id_to_exons[rec.transcript_id].append(rec)
        elif isinstance(rec, GtfCDS):
            transcript_id_to_cds_recs[rec.transcript_id].append(rec)

def _add_children(transcript, exons, cds_recs):
    if exons is not None:
        for exon in exons:
            transcript.add_exon(exon)
    if cds_recs is not None:
        for cds_rec in cds_recs:
            transcript.add_cds(cds_rec)
    transcript.sort_children()

def _resolve_gtf_records(gtf_data, transcript_id_to_exons, transcript_id_to_cds_recs):
    """add exon and CDS records to transcripts and sort"""
    for transcript in gtf_data.transcripts:
        _add_children(transcript,
                      transcript_id_to_exons.get(transcript.transcript_id),
                      transcript_id_to_cds_recs.get(transcript.transcript_id))

def gtf_data_parser(gtf_file, *, include_features: StrSetNone = None, attrs: GtfAttrsSet = GtfAttrsSet.FLAIR):
    """parse a GTF file into a GtfData object.  Use attrs=GtfAttrsSet.FLAIR
    to only parse the four attributes used by flair for faster loading."""
    # must save up exons and CDS, as sorting of GTF files is not required
    gtf_data = GtfData()
    transcript_id_to_exons = defaultdict(list)
    transcript_id_to_cds_recs = defaultdict(list)
    _load_gtf_records(gtf_file, gtf_data, transcript_id_to_exons, transcript_id_to_cds_recs, include_features, _ATTRS_SET_TO_PARSER[attrs])
    _resolve_gtf_records(gtf_data, transcript_id_to_exons, transcript_id_to_cds_recs)
    return gtf_data

def gtf_write_row(gtf_fh, chrom, source, feature, start, end, score, strand, frame, *,
                  gene_id=None, gene_name=None, transcript_id=None, attrs=None):
    """
    Write a row of a GTF without creating a GtfRecord object.  Handle score, strand, frame being
    None.
    """
    # Individual attributes override attrs
    merged = attrs.copy() if attrs is not None else {}
    if gene_id is not None:
        merged['gene_id'] = gene_id
    if gene_name is not None:
        merged['gene_name'] = gene_name
    if transcript_id is not None:
        merged['transcript_id'] = transcript_id

    print(_format_gtf_columns(chrom, source, feature, start, end, score, strand, frame, merged),
          file=gtf_fh)
