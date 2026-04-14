"""
Genome partitioning and parallel execution.

Manages per-region Partition objects, each with a temporary working directory
and subsets of annotation data pickled to disk, and runs a user-supplied
function on each partition.
"""
import os
import re
import pickle
import multiprocessing as mp
import pipettor
from flair import SeqRange
from flair.pycbio.hgdata.bed import BedReader

_GTF_DATA_PKL = 'gtf_data.pkl'
_JUNCTION_CORRECTOR_PKL = 'junction_corrector.pkl'


def parallel_mode_parse(parser, parallel_mode):
    """Parse --parallel_mode option string into a tuple.

    Valid values: auto:10GB, bychrom, byregion
    Returns ('auto', size_gb), ('bychrom', None), or ('byregion', None).
    """
    match = re.match(r'^(auto):(\d+)GB$|^(bychrom|byregion)$', parallel_mode)
    if match is None:
        parser.error(f"Invalid value for --parallel_mode: '{parallel_mode}', expected auto:10GB, bychrom, or byregion")
    if match.group(1) is not None:
        size = int(match.group(2))
        if size < 1:
            parser.error("auto parallel_mode must have a size greater than zero")
        return (match.group(1), size)
    else:
        return (match.group(3), None)


class Partition:
    """A single genome region with annotation data pickled to its temporary directory.

    Attributes:
        region:    SeqRange for this partition
        temp_dir:  per-partition working directory (created on construction)

    GtfData and IntronSupport are pickled to temp_dir at construction and
    loaded on demand by the partition runner wrapper.
    """
    def __init__(self, region, temp_dir, *, gtf_data=None, junction_corrector=None):
        self.region = region
        self.temp_dir = temp_dir
        os.makedirs(temp_dir, exist_ok=True)
        self._pickle(gtf_data, _GTF_DATA_PKL)
        self._pickle(junction_corrector, _JUNCTION_CORRECTOR_PKL)

    def _pickle(self, obj, name):
        if obj is not None:
            with open(self.temp_path(name), 'wb') as fh:
                pickle.dump(obj, fh)

    def _unpickle(self, name):
        path = self.temp_path(name)
        if os.path.exists(path):
            with open(path, 'rb') as fh:
                return pickle.load(fh)
        return None

    def load_gtf_data(self):
        """Load and return the pickled GtfData, or None if not present."""
        return self._unpickle(_GTF_DATA_PKL)

    def load_junction_corrector(self):
        """Load and return the pickled IntronSupport, or None if not present."""
        return self._unpickle(_JUNCTION_CORRECTOR_PKL)

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
    partition, func, func_kwargs = packed
    func(partition=partition,
         gtf_data=partition.load_gtf_data(),
         junction_corrector=partition.load_junction_corrector(),
         **func_kwargs)


def _run_flair_partition(genome_aligned_bam, annot_gtf, threads):
    """Run flair_partition and return a list of SeqRanges."""
    cmd = ['flair_partition',
           '--min_partition_items=1000',
           f'--threads={threads}',
           f'--bam={genome_aligned_bam}']
    if annot_gtf is not None:
        cmd += [f'--gtf={annot_gtf}']
    cmd += ['/dev/stdout']
    with pipettor.Popen(cmd) as fh:
        return [SeqRange(bed.chrom, bed.chromStart, bed.chromEnd)
                for bed in BedReader(fh)]


def _decide_parallel_mode(parallel_mode, genome_aligned_bam):
    # FIXME: remove by-chrom, decide based on really just need a size
    # if size exceeds chromosome, it still works.
    if parallel_mode[0] in ('bychrom', 'byregion'):
        return parallel_mode[0]
    # auto: choose based on BAM file size
    file_size_gb = os.path.getsize(genome_aligned_bam) / 1e9
    return 'byregion' if file_size_gb > parallel_mode[1] else 'bychrom'


class PartitionRunner:
    """A set of genome partitions, each with a temporary directory and pickled annotation data.

    Construction subsets and pickles annotation data per region and creates the
    per-partition temp directories under work_dir.  Call run() to apply a function
    to each partition.
    """
    def __init__(self, regions, work_dir, *, gtf_data=None, junction_corrector=None, threads=1):
        """
        Args:
            regions: iterable of SeqRange objects
            work_dir: root directory; per-partition subdirectories are created here
            gtf_data: GtfData to subset and pickle per region, or None
            junction_corrector:  JunctionCorrector to subset and pickle per region, or None
            threads: number of parallel workers used by run()
        """
        os.makedirs(work_dir, exist_ok=True)
        self.work_dir = work_dir
        self.threads = threads
        self.partitions = []
        region_gtf = region_jc = None
        for region in regions:
            if gtf_data is not None:
                region_gtf = gtf_data.subset_for_region(region.name, region.start, region.end)
            if junction_corrector is not None:
                region_jc = junction_corrector.subset_for_region(region.name, region.start, region.end)
            self.partitions.append(Partition(region, _region_temp_dir(work_dir, region),
                                             gtf_data=region_gtf, junction_corrector=region_jc))

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
        packed = [(p, func, kwargs) for p in self.partitions]

        if self.threads == 1:
            for p in packed:
                _call_partition_func(p)
        else:
            mp.set_start_method('fork', force=True)
            with mp.Pool(self.threads) as pool:
                pool.map(_call_partition_func, packed)


def _region_temp_dir(work_dir, region):
    return os.path.join(work_dir, f"{region.name}-{region.start}-{region.end}")


def partition_runner_factory(parallel_mode, genome, genome_aligned_bam, work_dir, annot_gtf, threads, *,
                             gtf_data=None, junction_corrector=None):
    """Create a PartitionRunner, choosing bychrom or byregion based on parallel_mode.

    parallel_mode is a tuple as returned by parallel_mode_parse:
        ('bychrom', None) | ('byregion', None) | ('auto', size_gb)
    """
    if _decide_parallel_mode(parallel_mode, genome_aligned_bam) == 'bychrom':
        regions = [SeqRange(chrom, 0, genome.get_reference_length(chrom))
                   for chrom in genome.references]
    else:
        regions = _run_flair_partition(genome_aligned_bam, annot_gtf, threads)
    return PartitionRunner(regions, work_dir, gtf_data=gtf_data, junction_corrector=junction_corrector, threads=threads)
