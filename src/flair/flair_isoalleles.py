"""flair isoalleles: each isoform's reads by haplotype, and in a tumor and normal pair
by somatic variants, and the sequence and protein of each, from reads haplotagged by
an upstream phasing tool and the read-to-isoform map of flair quantify.

See isoallele_genes for how reads are grouped, isoallele_variants for the variants,
and isoalleles_tsv for the outputs.  The genes are processed in batches of loci, as
flair quantify counts them (make_batches), by a pool of workers each with its own
open files (WorkerContext).
"""
import inspect
import logging
import multiprocessing as mp
import os
import shutil
from contextlib import ExitStack
import pysam
from flair import FlairInputDataError
from flair.aaseq_tsv import write_aaseqs
from flair.counts_matrix_tsv import CountsRow, counts_matrix_writer
from flair.flair_quantify import load_isoform_data, make_batches
from flair.io_utils import make_run_temp_dir, open_sorted_bed
from flair.isoallele_genes import GeneOutput, process_gene
from flair.isoallele_variants import open_vcf, variant_keys
from flair.isoalleles_tsv import TABLES, isoalleles_table_writer, table_path
from flair.isoform_data import make_big_bed
from flair.output_compress import compressed_path, open_fasta_output
from flair.read_map_index import ReadMapIndex
from flair.read_map_tsv import ReadMapWriter

# batches of genes per thread, so threads finishing early take more
BATCHES_PER_THREAD = 4
# the tag of a labelled read: its isoform and group
LABEL_TAG = 'IA'

SOMATIC_MODES = ('tumor_vcf', 'somatic_vcf', 'variant')


def add_subparser(subparsers):
    desc = ("Each isoform's reads by haplotype, and in a tumor and normal pair by somatic variants, with the "
            'sequence and protein of each, from reads haplotagged by a phasing tool (such as WhatsHap) and the '
            'read map of flair quantify.  One sample: --bam, --read_map and --phased_vcf.  A pair adds '
            '--normal_bam and --normal_read_map, --phased_vcf being the normal\'s, and one way of choosing the '
            'somatic variants: --tumor_vcf, --somatic_vcf or --variant')
    parser = subparsers.add_parser('isoalleles', help="Isoform usage by haplotype and somatic variants",
                                   description=desc)
    required = parser.add_argument_group('required named arguments')
    required.add_argument('--isoform_bed', required=True,
                          help='the isoforms, the flair isoform bed the read maps were made with by flair quantify; '
                               'its CDS gives each isoform\'s start codon.  Fusion isoforms are left out')
    required.add_argument('-g', '--genome', required=True, help='FASTA of the genome, indexed with samtools faidx')
    required.add_argument('--bam', required=True,
                          help='the sample\'s (in a pair, the tumor\'s) reads aligned to the genome, sorted and '
                               'indexed, haplotagged (HP and PS tags) by the phasing of --phased_vcf')
    required.add_argument('--read_map', required=True,
                          help='the read map flair quantify --generate_map wrote for --bam\'s sample, '
                               '<output>.<sample>.read.map.tsv, plain, bgzip or gzip')
    required.add_argument('--phased_vcf', required=True,
                          help='the germline variants, phased: the VCF the phasing tool wrote and haplotagged '
                               '--bam with, the sample\'s own, or in a pair the normal\'s.  bgzipped and indexed')
    parser.add_argument('--germline_vcf',
                        help='the germline variant calls, of the sample, or in a pair the normal, when --phased_vcf '
                             'holds only some of them, as when only hets were phased: every passing call is a '
                             'germline variant, phased by --phased_vcf.  bgzipped and indexed (default: --phased_vcf)')
    parser.add_argument('--sample_name', help='the name of --bam\'s sample in the outputs (default: its file name)')
    pair = parser.add_argument_group('tumor and normal pair')
    pair.add_argument('--normal_bam',
                      help='the normal\'s reads, haplotagged by the same phasing as --bam; makes the run a pair, '
                           'whose tumor reads, of --bam, are grouped by somatic variants too')
    pair.add_argument('--normal_read_map', help='the read map of --normal_bam\'s sample, from the same flair quantify')
    pair.add_argument('--normal_name', help='the name of the normal in the outputs (default: its file name)')
    somatic = pair.add_mutually_exclusive_group()
    somatic.add_argument('--tumor_vcf',
                         help='the tumor\'s variant calls: its passing calls that --phased_vcf, the normal\'s, '
                              'does not call are somatic if at least --min_normal_depth normal reads inform the '
                              'site, and are otherwise left out, in excluded.tsv.  bgzipped and indexed')
    somatic.add_argument('--somatic_vcf',
                         help='the somatic variants, every passing call of a somatic caller\'s VCF (such as '
                              'ClairS or DeepSomatic), as it calls them.  bgzipped and indexed')
    somatic.add_argument('--variant', action='append',
                         help='a somatic variant, CHROM:POS:REF>ALT, or CHROM:POS for every alt the VCFs have there; '
                              'may be given more than once.  Only the variants given are somatic, whatever the '
                              'VCFs call, as for a known driver the normal calls because tumor cells are in it')
    pair.add_argument('--min_normal_depth', type=int, default=10,
                      help='with --tumor_vcf, the normal reads informing a site for a tumor-only call there to be '
                           'somatic (default: %(default)s)')
    parser.add_argument('--carrier_min', type=float, default=0.9,
                        help='an unphased germline het is put on the haplotype holding at least this share of an '
                             'isoform\'s tagged reads carrying it (default: %(default)s)')
    parser.add_argument('--min_tagged', type=int, default=10,
                        help='the tagged reads carrying an unphased germline het needed to place it (default: '
                             '%(default)s)')
    parser.add_argument('-o', '--output', default='flair', help='output file name base (default: %(default)s)')
    parser.add_argument('-t', '--threads', type=int, default=4,
                        help='batches of genes processed at once (default: %(default)s)')
    parser.add_argument('--generate_map', action='store_true',
                        help='write each sample\'s reads by isoform and group, '
                             '<output>.isoalleles.<sample>.read.map.tsv')
    parser.add_argument('--label_bams', action='store_true',
                        help=f'write each sample\'s assigned reads with their isoform and group in the {LABEL_TAG} '
                             f'tag, <output>.isoalleles.<sample>.bam (in IGV, group by tag {LABEL_TAG})')
    parser.add_argument('--temp_dir',
                        help='directory for temporary files, in a directory of the run\'s own (default: $TMPDIR)')
    parser.set_defaults(entry=isoalleles_cmd)
    return parser


def _sample_name(name, bam):
    base = os.path.basename(bam)
    return name or (base[:-len('.bam')] if base.endswith('.bam') else base)


def check_options(*, bam, sample_name, normal_bam, normal_read_map, normal_name, somatic, paths):
    """the options of a run of one sample or of a pair, consistent: somatic gives each
    of SOMATIC_MODES's option, paths every input file.  Returns the samples' names,
    the tumor first, and the somatic mode"""
    pair = normal_bam is not None
    given = [mode for mode in SOMATIC_MODES if somatic[mode]]
    if pair and not normal_read_map:
        raise FlairInputDataError('--normal_bam needs --normal_read_map, the normal\'s flair quantify read map')
    if not pair:
        pair_options = [('normal_read_map', normal_read_map), ('normal_name', normal_name), *somatic.items()]
        extra = [f'--{opt}' for opt, value in pair_options if value]
        if extra:
            raise FlairInputDataError(f'{", ".join(extra)} are of a tumor and normal pair: give --normal_bam too, or '
                                      'leave them out for one sample')
    elif not given:
        raise FlairInputDataError('a pair needs the somatic variants: --tumor_vcf (the tumor\'s calls, those the '
                                  'normal doesn\'t call), --somatic_vcf (a somatic caller\'s) or --variant')
    names = [_sample_name(sample_name, bam)] + ([_sample_name(normal_name, normal_bam)] if pair else [])
    if len(set(names)) != len(names):
        raise FlairInputDataError(f'the tumor and normal are both named {names[0]}: give --sample_name and '
                                  '--normal_name')
    for path in paths:
        if path and not os.path.exists(path):
            raise FlairInputDataError(f'{path} does not exist')
    return names, (given[0] if given else None)


class WorkerContext:
    """One process's open files and the settings process_gene reads: samples in the
    order their reads are scanned, the normal first; output_samples in the order of
    the outputs' columns, the tumor first"""

    def __init__(self, settings):
        self.__dict__.update(settings)
        self.genome = pysam.FastaFile(self.genome_path)
        self.bams = {name: pysam.AlignmentFile(bam) for name, bam in self.bam_paths.items()}
        self.map_readers = {name: index.open_reader() for name, index in self.map_indexes.items()}
        self.phased_vcf = open_vcf(self.phased_vcf_path)
        self.germline_vcf = open_vcf(self.germline_vcf_path) if self.germline_vcf_path else None
        self.tumor_vcf = open_vcf(self.tumor_vcf_path) if self.tumor_vcf_path else None
        self.somatic_vcf = open_vcf(self.somatic_vcf_path) if self.somatic_vcf_path else None


_context = None


def _init_worker(settings):
    global _context
    _context = WorkerContext(settings)


def _process_batch(packed):
    "a batch's genes, as one GeneOutput, and the paths of its labelled reads"
    batch_genes, batch_dir = packed
    out = GeneOutput()
    for gene_id, isos in batch_genes:
        process_gene(_context, gene_id, isos, out)
    label_bams = _write_batch_labels(_context, batch_genes, out, batch_dir) if _context.label_bams else {}
    if not _context.generate_map:
        out.labels = {}
    return out, label_bams


def _write_batch_labels(ctx, batch_genes, out, batch_dir):
    """each sample's reads of the batch's genes, with LABEL_TAG, in a sorted bam of
    the batch in batch_dir"""
    paths = {}
    for sample, labelled in out.labels.items():
        labels = dict(labelled)
        if not labels:
            continue
        unsorted = os.path.join(batch_dir, f'{sample}.unsorted.bam')
        written = set()
        with pysam.AlignmentFile(unsorted, 'wb', template=ctx.bams[sample]) as fh:
            for _, isos in batch_genes:
                chrom, start, end = isos[0].chrom, min(i.chromStart for i in isos), max(i.chromEnd for i in isos)
                for read in ctx.bams[sample].fetch(chrom, start, end):
                    label = labels.get(read.query_name)
                    if label and read.query_name not in written and not (read.is_secondary or read.is_supplementary):
                        read.set_tag(LABEL_TAG, label)
                        fh.write(read)
                        written.add(read.query_name)
        paths[sample] = os.path.join(batch_dir, f'{sample}.bam')
        pysam.sort('-o', paths[sample], unsorted)
        os.remove(unsorted)
    return paths


class Outputs:
    """the outputs, open in stack, written batch by batch in the order of the batches,
    the proteins once all are"""

    def __init__(self, stack, output, samples, generate_map):
        self.writers = {table: stack.enter_context(isoalleles_table_writer(compressed_path(table_path(output, table)),
                                                                           table))
                        for table in TABLES}
        self.counts = stack.enter_context(counts_matrix_writer(compressed_path(f'{output}.isoalleles.counts.tsv'),
                                                               samples))
        self.fasta = stack.enter_context(open_fasta_output(f'{output}.isoalleles.fa'))
        self.bed = stack.enter_context(open_sorted_bed(compressed_path(f'{output}.isoalleles.bed')))
        self.proteins = {}
        self.read_maps = ({s: stack.enter_context(ReadMapWriter(compressed_path(f'{output}.isoalleles.{s}.read.map.tsv')))
                           for s in samples} if generate_map else {})

    def add(self, out):
        for table in TABLES:
            for row in getattr(out, table):
                self.writers[table].writeRow(row)
        for gene_id, isoform_id, counts in out.counts:
            self.counts.writeCountsRow(CountsRow(gene_id, isoform_id, counts))
        self.fasta.writelines(out.fasta)
        self.bed.writelines(out.bed)
        for pid, (aaseq, _) in out.proteins.items():
            self.proteins[aaseq] = pid
        for sample, writer in self.read_maps.items():
            by_label = {}
            for read, label in out.labels.get(sample, ()):
                by_label.setdefault(label, []).append(read)
            for label in sorted(by_label):
                writer.writeReads(label, sorted(by_label[label]))


def isoalleles_cmd(args):
    "the subcommand: its own options, not flair's global ones, to isoalleles"
    params = inspect.signature(isoalleles).parameters
    isoalleles(**{k: v for k, v in vars(args).items() if k in params})


def isoalleles(*, isoform_bed, genome, bam, read_map, phased_vcf, germline_vcf=None, sample_name=None, normal_bam=None,
               normal_read_map=None, normal_name=None, tumor_vcf=None, somatic_vcf=None, variant=None,
               min_normal_depth=10, carrier_min=0.9, min_tagged=10, output='flair', threads=4,
               generate_map=False, label_bams=False, temp_dir=None):
    names, somatic_mode = check_options(
        bam=bam, sample_name=sample_name, normal_bam=normal_bam, normal_read_map=normal_read_map,
        normal_name=normal_name, somatic={'tumor_vcf': tumor_vcf, 'somatic_vcf': somatic_vcf, 'variant': variant},
        paths=[bam, read_map, phased_vcf, germline_vcf, isoform_bed, genome, normal_bam, normal_read_map, tumor_vcf,
               somatic_vcf])
    gene_data, fusions = load_isoform_data(isoform_bed)
    if fusions:
        logging.info(f'{len(fusions)} fusion isoforms left out: isoalleles takes the isoforms of genes')
    run_temp_dir = make_run_temp_dir(output, temp_dir)
    try:
        tumor, normal = names[0], (names[1] if len(names) > 1 else None)
        bam_paths = {tumor: bam, **({normal: normal_bam} if normal else {})}
        map_paths = {tumor: read_map, **({normal: normal_read_map} if normal else {})}
        logging.info('indexing the read maps')
        map_indexes = {name: ReadMapIndex(path, run_temp_dir) for name, path in map_paths.items()}
        _check_read_maps(map_indexes, map_paths, gene_data)
        keys = variant_keys(variant, [tumor_vcf, somatic_vcf, germline_vcf, phased_vcf]) if variant else []
        settings = dict(genome_path=genome, bam_paths=bam_paths, map_indexes=map_indexes,
                        phased_vcf_path=phased_vcf, germline_vcf_path=germline_vcf, tumor_vcf_path=tumor_vcf,
                        somatic_vcf_path=somatic_vcf,
                        somatic_mode=somatic_mode, variant_keys=keys, min_normal_depth=min_normal_depth,
                        carrier_min=carrier_min, min_tagged=min_tagged, normal=normal,
                        samples=([normal, tumor] if normal else [tumor]), output_samples=names,
                        keep_labels=generate_map or label_bams, generate_map=generate_map, label_bams=label_bams)
        batches = make_batches(gene_data, threads * BATCHES_PER_THREAD)
        packed = []
        for n, batch in enumerate(batches):
            batch_dir = os.path.join(run_temp_dir, f'batch{n}')
            os.makedirs(batch_dir)
            packed.append(([(gene_id, gene_data[gene_id].isoform_beds) for gene_id in batch['gene_ids']], batch_dir))
        logging.info(f'{len(gene_data)} genes in {len(batches)} batches')
        _run(settings, packed, threads, output, names, generate_map, label_bams, genome)
    finally:
        shutil.rmtree(run_temp_dir, ignore_errors=True)


def _check_read_maps(map_indexes, map_paths, gene_data):
    "every isoform a read map names is one of the bed's"
    names = {bed.name for gene_info in gene_data.values() for bed in gene_info.isoform_beds}
    for sample, index in map_indexes.items():
        unknown = sorted(set(index.offsets) - names)
        if unknown:
            raise FlairInputDataError(f'{map_paths[sample]} names isoforms not in --isoform_bed, such as '
                                      f'{", ".join(unknown[:3])}: give the read map flair quantify made with this '
                                      'isoform bed')


def _run(settings, packed, threads, output, names, generate_map, label_bams, genome):
    label_parts = {name: [] for name in names}
    with ExitStack() as stack:
        outputs = Outputs(stack, output, names, generate_map)
        if threads == 1:
            _init_worker(settings)
            results = map(_process_batch, packed)
        else:
            mp.set_start_method('fork', force=True)
            pool = stack.enter_context(mp.Pool(threads, initializer=_init_worker, initargs=(settings,)))
            results = pool.imap(_process_batch, packed)
        for out, label_bams_of_batch in results:
            outputs.add(out)
            for sample, path in label_bams_of_batch.items():
                label_parts[sample].append(path)
    write_aaseqs(compressed_path(f'{output}.isoalleles.aaseq.tsv'), outputs.proteins)
    with pysam.FastaFile(genome) as genome_fa:
        make_big_bed(genome_fa, f'{output}.isoalleles.chrom.sizes', f'{output}.isoalleles')
    if label_bams:
        for sample, parts in label_parts.items():
            path = f'{output}.isoalleles.{sample}.bam'
            if parts:
                pysam.merge('-f', '-o', path, *parts)
                pysam.index(path)
