"""The --region option: a genome range written on the command line, in either the
one-based, closed form of a genome browser or the zero-based, half-open form of a
BED line.  Combined with the regions of a --region_bed file."""

import re
from flair import SeqRange, FlairInputDataError, combine_regions
from flair.pycbio.hgdata.bed import BedReader

# appended to the help of each program's --region option
REGION_FORMAT_HELP = ("A region is either one-based and closed, chr1:1,000-2,000, or "
                      "zero-based and half-open, chr1 999 2000, as a BED line.  May be "
                      "repeated, and may be used with --region_bed, whose regions are "
                      "added to it")

# `chr1:1,000-2,000'; the sequence name may contain a colon, so the last one starts the range
_ONE_BASED_RE = re.compile(r"^(\S+):([\d,]+)-([\d,]+)$")

# `chr1 999 2000', a BED line, so also tab-separated
_ZERO_BASED_RE = re.compile(r"^(\S+)\s+([\d,]+)\s+([\d,]+)$")


def _coord(value):
    return int(value.replace(',', ''))


def _bad_region(spec, problem):
    return FlairInputDataError(f"invalid --region `{spec}': {problem}; write it as one-based, "
                               f"closed `chr1:1,000-2,000' or zero-based, half-open `chr1 999 2000'")


def _make_region(spec, match, base):
    "SeqRange of a matched region, base 1 for a browser range, 0 for a BED one"
    start, end = _coord(match.group(2)) - base, _coord(match.group(3))
    if start < 0:
        raise _bad_region(spec, f"the start must be at least {base}")
    if end <= start:
        raise _bad_region(spec, "the end must be after the start")
    return SeqRange(match.group(1), start, end)


def parse_region(spec):
    "the SeqRange of one --region"
    spec = spec.strip()
    match = _ONE_BASED_RE.match(spec)
    if match is not None:
        return _make_region(spec, match, 1)
    match = _ZERO_BASED_RE.match(spec)
    if match is not None:
        return _make_region(spec, match, 0)
    raise _bad_region(spec, "not a sequence name with a range")


def read_region_bed(region_bed):
    "the SeqRanges of a region BED, unmerged"
    return [SeqRange(bed.chrom, bed.chromStart, bed.chromEnd) for bed in BedReader(region_bed)]


def parse_regions(region_specs, region_bed=None):
    """The regions of the --region options and the --region_bed file, merged, or None
    if neither was given"""
    if not region_specs and (region_bed is None):
        return None
    regions = [parse_region(spec) for spec in region_specs]
    if region_bed is not None:
        regions += read_region_bed(region_bed)
    if len(regions) == 0:
        raise FlairInputDataError(f"no regions in --region_bed {region_bed}; it must have at least "
                                  f"one BED record, or use --region instead")
    return combine_regions(regions)


def write_regions(regions, fh):
    "the regions as BED3 lines on an open file"
    for region in regions:
        fh.write(f"{region.name}\t{region.start}\t{region.end}\n")


def write_regions_bed(regions, bed_file):
    "the regions as a BED3, for a program that takes a region BED; returns the path"
    with open(bed_file, 'w') as fh:
        write_regions(regions, fh)
    return bed_file
