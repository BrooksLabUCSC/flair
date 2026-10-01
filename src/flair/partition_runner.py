"""
Genome partitioning and parallel execution.

Manages per-region Partition objects, each with a temporary working directory
and the annotation data for its region, and runs a user-supplied function on
each partition.
"""
import gc
import os
import re
import shutil
import multiprocessing as mp
import pipettor
import pysam
from flair import SeqRange, FlairInputDataError
from flair import thread_share
from flair.pycbio.hgdata.bed import BedReader

# the full annotation data of each PartitionRunner, keyed by its work_dir.  The
# pool workers are forked, so they inherit this and subset their own region,
# rather than the main process subsetting and pickling every region up front,
# which took minutes, or pickling the full data into every task
_shared_data = {}


def combine_temp_files_by_suffix(output, temp_prefixes, suffixes):
    for filesuffix in suffixes:
        with open(output + filesuffix, 'wb') as combined_fh:
            for temp_prefix in temp_prefixes:
                with open(temp_prefix + filesuffix, 'rb') as in_fh:
                    shutil.copyfileobj(in_fh, combined_fh, 1024 * 1024 * 10)


def parallel_mode_parse(parallel_mode):
    """Parse --parallel_mode option string into a tuple.

    Valid values: auto:10GB, bychrom, byregion
    Returns ('auto', size_gb), ('bychrom', None), or ('byregion', None).
    """
    match = re.match(r'^(auto):(\d+)GB$|^(bychrom|byregion)$', parallel_mode)
    if match is None:
        raise FlairInputDataError(f"Invalid value for --parallel_mode: '{parallel_mode}', expected auto:10GB, bychrom, or byregion")
    if match.group(1) is not None:
        size = int(match.group(2))
        if size < 1:
            raise FlairInputDataError("auto parallel_mode must have a size greater than zero")
        return (match.group(1), size)
    else:
        return (match.group(3), None)


class Partition:
    """A single genome region with a temporary directory.

    Attributes:
        region:    SeqRange for this partition
        temp_dir:  per-partition working directory (created on construction)

    The region's GtfData and JunctionCorrector are subset from the runner's full
    data on demand, in the process that runs the partition.
    """
    def __init__(self, region, temp_dir, data_key):
        self.region = region
        self.temp_dir = temp_dir
        self._data_key = data_key
        os.makedirs(temp_dir, exist_ok=True)

    def _subset(self, name):
        data = _shared_data.get(self._data_key, {}).get(name)
        if data is None:
            return None
        return data.subset_for_region(self.region.name, self.region.start, self.region.end)

    def load_gtf_data(self):
        """Return the GtfData for this region, or None if the runner has none."""
        return self._subset('gtf_data')

    def load_junction_corrector(self):
        """Return the JunctionCorrector for this region, or None if the runner has none."""
        return self._subset('junction_corrector')

    def temp_path(self, suffix):
        """Return a path inside this partition's temp_dir with the given suffix."""
        return os.path.join(self.temp_dir, suffix)

    @property
    def file_prefix(self):
        """Base path prefix for output files inside temp_dir."""
        r = self.region
        return self.temp_path(f"{r.name}-{r.start}-{r.end}")

    def output_path(self, name):
        """Return the path for a named output file inside temp_dir.

        E.g. partition.output_path('reads.fasta')
             -> temp_dir/chr12-0-133275309.reads.fasta
        """
        return self.file_prefix + '.' + name

    def __repr__(self):
        return f"Partition({self.region.name}:{self.region.start}-{self.region.end}, temp_dir={self.temp_dir!r})"


def _call_partition_func(packed):
    partition, weight, func, func_kwargs = packed
    thread_share.task_start(weight)
    try:
        func(partition=partition,
             gtf_data=partition.load_gtf_data(),
             junction_corrector=partition.load_junction_corrector(),
             **func_kwargs)
    finally:
        thread_share.task_done()


def _run_flair_partition(genome_aligned_bam, annot_gtf, threads):
    """Run flair_partition and return a list of SeqRanges and a list of the number
    of input items in each."""
    cmd = ['flair_partition',
           '--min_partition_items=1000',
           f'--threads={threads}',
           f'--bam={genome_aligned_bam}']
    if annot_gtf is not None:
        cmd += [f'--gtf={annot_gtf}']
    cmd += ['/dev/stdout']
    regions, item_counts = [], []
    with pipettor.Popen(cmd) as fh:
        for bed in BedReader(fh, numStdCols=4):
            regions.append(SeqRange(bed.chrom, bed.chromStart, bed.chromEnd))
            item_counts.append(int(bed.extraCols[0]))
    return regions, item_counts


def _decide_parallel_mode(parallel_mode, genome_aligned_bam):
    # FIXME: remove by-chrom, decide based on really just need a size
    # if size exceeds chromosome, it still works.
    if parallel_mode[0] in ('bychrom', 'byregion'):
        return parallel_mode[0]
    # auto: choose based on BAM file size
    file_size_gb = os.path.getsize(genome_aligned_bam) / 1e9
    return 'byregion' if file_size_gb > parallel_mode[1] else 'bychrom'


class PartitionRunner:
    """A set of genome partitions, each with a temporary directory.

    Construction creates the per-partition temp directories under work_dir and
    keeps the annotation data, which each partition subsets for its region when it
    runs.  Call run() to apply a function to each partition.
    """
    def __init__(self, regions, work_dir, *, gtf_data=None, junction_corrector=None, threads=1, weights=None):
        """
        Args:
            regions: iterable of SeqRange objects
            work_dir: root directory; per-partition subdirectories are created here
            gtf_data: GtfData to subset per region, or None
            junction_corrector:  JunctionCorrector to subset per region, or None
            threads: number of parallel workers used by run()
            weights: optional size of each region, such as its number of reads;
                larger regions run first and get more of the idle threads
        """
        os.makedirs(work_dir, exist_ok=True)
        self.work_dir = work_dir
        self.threads = threads
        regions = list(regions)
        self.weights = list(weights) if weights is not None else [1] * len(regions)
        # indexes are built here, in the parent, so forked workers share them
        # rather than each building, and so copying, its own
        for data in (gtf_data, junction_corrector):
            if data is not None:
                data.build_indexes()
        _shared_data[work_dir] = {'gtf_data': gtf_data, 'junction_corrector': junction_corrector}
        self.partitions = [Partition(region, _region_temp_dir(work_dir, region), work_dir)
                           for region in regions]

    def __iter__(self):
        return iter(self.partitions)

    def __len__(self):
        return len(self.partitions)

    def run(self, func, **kwargs):
        """Run func for each partition, in parallel if threads > 1 (set at construction).

        func is called with the following keyword arguments:
            partition      -- Partition for this region; provides region,
                              temp_dir, output_path(), and file_prefix
            gtf_data       -- GtfData subset for the region, or None
            junction_corrector -- JunctionCorrector subset for the region, or None
            **kwargs       -- any additional keyword arguments passed to run()

        Side effects (e.g. writing files to partition.temp_dir) are the
        expected pattern; return values from func are discarded.
        """
        # largest partitions first, handed out one at a time, so the big ones don't
        # end up running alone at the end; programs run by the partitions, such as
        # minimap2, borrow threads left idle, more for larger partitions
        packed = sorted(((p, w, func, kwargs) for p, w in zip(self.partitions, self.weights)),
                        key=lambda x: x[1], reverse=True)
        thread_share.init(self.threads, max(self.weights, default=1))

        if self.threads == 1:
            for p in packed:
                _call_partition_func(p)
        else:
            mp.set_start_method('fork', force=True)
            # objects existing at the fork go to the collector's permanent
            # generation, so collections in the workers don't visit them; visiting
            # writes to each object, which copied the whole annotation into every
            # worker
            gc.freeze()
            try:
                with mp.Pool(self.threads) as pool:
                    for _ in pool.imap_unordered(_call_partition_func, packed, chunksize=1):
                        pass
            finally:
                gc.unfreeze()


def _count_chrom_reads(genome_aligned_bam, regions):
    "number of mapped reads on each whole-chromosome region, from the BAM index"
    with pysam.AlignmentFile(genome_aligned_bam, 'rb') as bam:
        mapped = {s.contig: s.mapped for s in bam.get_index_statistics()}
    return [mapped.get(r.name, 0) for r in regions]


def _region_temp_dir(work_dir, region):
    return os.path.join(work_dir, f"{region.name}-{region.start}-{region.end}")


def partition_regions(parallel_mode, genome, genome_aligned_bam, annot_gtf, threads):
    """Divide the genome into regions to run in parallel, choosing bychrom or byregion
    based on parallel_mode, a tuple as returned by parallel_mode_parse:
        ('bychrom', None) | ('byregion', None) | ('auto', size_gb)
    Returns a list of SeqRanges and a list of their sizes, for PartitionRunner."""
    if _decide_parallel_mode(parallel_mode, genome_aligned_bam) == 'bychrom':
        regions = [SeqRange(chrom, 0, genome.get_reference_length(chrom))
                   for chrom in genome.references]
        return regions, _count_chrom_reads(genome_aligned_bam, regions)
    else:
        return _run_flair_partition(genome_aligned_bam, annot_gtf, threads)
