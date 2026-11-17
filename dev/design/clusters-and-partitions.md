# Clusters and partitions

Two different things, often confused because both come out of `flair_partition` as
ranges on a chromosome.

## Cluster

A cluster is a maximal set of input items that overlap, transitively: reads of every
sample, and the spans of annotated transcripts when a GTF is given.  `flair_partition`
builds them by sorting the items and breaking wherever the next item starts after the
end of everything before it (`partition_reader`, `partition_build`).

The property that matters:

- no read and no annotated transcript crosses a cluster boundary
- so every transcript a read overlaps is in the read's own cluster

A cluster is a unit of meaning.  It is the largest region a read can influence and the
smallest region that contains everything that influences the read.

## Partition

A partition is one or more adjacent clusters run as one task.  Merging is for
scheduling: a task has a fork, an annotation subset and at least one minimap2 run
behind it, so very small tasks waste those, and very large ones leave threads idle.

- flair spliceevents merges clusters itself, up to `--partition_items` input items
  (`_work_units` in `flair_spliceevents.py`)
- flair transcriptome runs one cluster per partition

A partition is a unit of work.  Nothing about the result may depend on it.

## The rule

Anything that decides what a read becomes must be scoped to the cluster, never to the
partition:

- the annotation the reads are aligned to
- the junctions used to correct them
- the reads fetched
- any index built from them, such as the firstpass isoforms reads are realigned to

Scope one of them to the partition and results move when the partition size or the
parallel mode changes, or when a run is restricted to a region.  This is not
theoretical; both cases below were measured on A498 over hg19.

## How this went wrong before

`annot_data_from_gtf(gtf_data, region)` keeps every transcript of the `GtfData` it is
given and uses `region` only as a dict key.  Passing a cluster as `region` therefore
changes nothing; the `GtfData` itself has to be cut:

```python
annots = annot_data_from_gtf(gtf_data.subset_for_region(cluster.name, cluster.start, cluster.end), cluster)
```

Until that was done, flair spliceevents reported different events at different
`--partition_items`: on chr19, four rows of one gene appeared only when the partition
was large enough to pull a distant transcript into the alignment index.  With the
annotation and the junction corrector cut to the cluster, 1000 items and 20000 items
give identical tables, as does one task per cluster, over 91287 event rows genome-wide.

flair transcriptome had the same exposure through `--parallel_mode`: `bychrom` aligned
against a whole chromosome's transcripts and `byregion` against a partition's, and the
two gave different isoforms on chr19 (one isoform of 2559, six read counts).  It now
partitions by cluster always.

## Where the code is

| what                            | file                     | function                  |
|---------------------------------|--------------------------|---------------------------|
| clusters built                  | `flair_partition.py`     | `partition_reader`        |
| clusters expanded over a region | `flair_partition.py`     | `closed_partitions`       |
| clusters fetched by flair       | `partition_runner.py`    | `run_flair_partition`     |
| clusters run as tasks           | `partition_runner.py`    | `PartitionRunner`         |
| clusters merged into tasks      | `flair_spliceevents.py`  | `_work_units`             |
| one cluster corrected           | `flair_spliceevents.py`  | `_correct_cluster_sample` |
| one cluster's transcripts built | `flair_transcriptome.py` | `_run_region_reads`       |

## Regions

`--region` and `--region_bed` are a third thing: what the user asked about.  A region
is expanded over the items that overlap it until nothing more does
(`closed_partitions`), which yields clusters again.  Events and transcripts are then
reported against the region the user gave, not against the cluster.
