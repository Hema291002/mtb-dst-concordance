# Nextflow implementation validated against the original scripts

## What was compared

The bash and Python scripts in `bin/` produced the results reported in this
repository. The Nextflow pipeline in `main.nf` and `modules/local/` is a
separate implementation of the same per-sample path: trim, align, mark
duplicates, coverage, call, filter.

Both were run on the same 30 isolates from the same raw FASTQ files.

## Result

Filtered variant counts are identical for all 30 samples.

Position-level comparison for ERR181948 (chromosome, position, reference,
alternate) returns no differences across all 768 variants.

## Version difference

The two implementations do not use the same tool versions. The bash
environment has samtools 1.24 and bcftools 1.24; the Nextflow modules pin
1.21. The calls are identical regardless, indicating the results do not sit
on version-sensitive boundaries.

## Runtime

Nextflow: 3h 49m wall clock, 180 tasks, 14.4 CPU hours on 4 logical cores.
Bash: approximately 3h, sequential.

The Nextflow run is slower in wall clock on this hardware. The parallelism it
introduces has no room to help on a 2-core machine: 14.4 CPU hours in 3.8 wall
hours means the cores were saturated throughout, and the extra concurrency
produced contention rather than speedup. On a cluster the same dependency
graph would finish in a fraction of the time. Nextflow's benefit here is
portability, containerisation and content-addressed caching, not speed on this
laptop.

## What the pipeline does not cover

Metadata handling, sample selection, WHO catalogue interpretation, phenotype
resolution and lineage assignment remain as scripts in `bin/`. These are
one-off cohort-level analyses rather than per-sample steps, and forcing them
into processes would add complexity without adding capability. nf-core
pipelines make the same split.
