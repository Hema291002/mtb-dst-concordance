# Analysis status

Honest record of what has actually been run versus what is planned.
Updated as work proceeds. Nothing is marked complete until its output
has been inspected.

| Stage | Status | Notes |
|---|---|---|
| 1. Environment setup | complete | WSL2, ext4, conda env `mtb-meta` |
| 2. Metadata acquisition | complete | CRyPTIC v3.4.0, Zenodo 15680920 |
| 3. Sample selection | complete | 30 isolates, seed 42, balanced 8/7/7/8 |
| 4. Reference genome | complete | H37Rv GCF_000195955.2 (NC_000962.3), 4,411,532 bp; FASTA and GFF contig names verified to match |
| 5. FASTQ download | complete | 60 files, 9.2 GB, all MD5-verified. 2 files failed checksum on first download and were re-fetched. |
| 6. Raw read QC | complete | FastQC + MultiQC on 60 files. All pairs read-count matched. GC 63-66% across all samples, consistent with M. tuberculosis, no contamination signal. Estimated raw coverage 65-166x (median 110x). Read lengths heterogeneous: 80/100/150/151 bp. |
| 7. Trimming | complete | fastp 1.3.6, adapter auto-detect (Nextera found, not TruSeq), tail-only Q20 trimming, min length 50. All 30 retained, pass rate 91.9-99.7%. Post-trim coverage 64-144x. |
| 8. Alignment | complete | bwa-mem2 2.2.1 to H37Rv, coordinate-sorted, indexed, per-sample read groups. Mapping 96.5-99.9%, properly paired 89.5-99.4%. |
| 9. Post-alignment QC | complete | |
| 10. Variant calling | not started | |
| 11. Variant filtering | not started | |
| 12. Annotation | not started | |
| 13. WHO catalogue interpretation | not started | |
| 14. Phenotype comparison | not started | |
| 15. Lineage assignment | not started | |
| 16. Nextflow pipeline | not started | |
| 17. R analysis and figures | not started | |

## Held out until stage 14
- `data/meta/held_out/PREDICTIONS.parquet`
- `data/meta/held_out/EFFECTS.parquet`
- `GENOMES` columns: LINEAGE, SUBLINEAGE, ANTIBIOGRAM

## Known limitations, recorded as identified
- Resistance enriched by design in CRyPTIC; concordance estimates will
  not generalise to routine diagnostic case mix.
- Balanced design supports sensitivity and specificity only. PPV and NPV
  are prevalence-dependent and will not be reported as if generalisable.
- n=15 resistant per drug gives wide confidence intervals; the study
  cannot distinguish e.g. 90% from 98% sensitivity.
- Resistance category partially confounded with collecting site, and
  therefore with strain lineage.
- Samples span multiple processing pipeline builds (`PIPELINE_BUILD`).
- Two isolates (SRR2100603, SRR6045767) lack `center_name`; both are
  direct ENA submissions rather than CRyPTIC site collections.
- Read length varies across isolates (80-151 bp); shorter reads map less
  uniquely in repetitive regions, so effective resolution is not uniform.
- Some submitted FASTQs were already trimmed before deposition (mean read
  length below median in several runs), so preprocessing is not uniform
  across the cohort.
- FastQC "Per Base Sequence Content" fails in 46/60 files. This is expected
  for a 65% GC genome with standard library prep and was deliberately not
  acted upon.
- Library preparation differs systematically by era: 80 bp reads with ~80 bp
  inserts, 100 bp reads with ~165 bp inserts, 150 bp reads with ~250 bp
  inserts. Longer inserts resolve repetitive regions better, so mapping
  resolution is not uniform across the cohort.
- Adapter contamination ranged 0-20% across libraries. Nextera adapters were
  detected, not TruSeq; a hard-coded TruSeq sequence would have removed none.
- Mapping rates 96.5-99.9%. Three isolates (ERR13273273, ERR4818219,
  ERR4821953) map ~1-3% lower than the rest, consistent with reference bias
  against strains divergent from H37Rv (lineage 4). Accessory sequence absent
  from H37Rv cannot be assessed by reference-based calling.
- ERR8975061 has zero coverage across ethA (4326003-4327673) despite 103x
  genome-wide depth, consistent with a whole-gene deletion. Reference-based
  variant calling cannot represent deletions of this scale as variant calls;
  they are detectable only through coverage analysis. Any gene not examined
  for coverage could carry an undetected deletion.
- rrs coverage at 30x is incomplete in five isolates (43.5-96.3%), all short-read
  libraries. Aminoglycoside resistance calls from rrs would be unreliable in
  these samples. INH/RIF target genes are unaffected.
- All 30 isolates have 100% coverage at 30x across katG, rpoB, fabG1, inhA and
  ahpC including 200 bp promoter flanks, so wild-type calls for the two target
  drugs are supported by data rather than absence of data.
- Gene-scale deletions are invisible to variant calling. One was found by
  coverage analysis (see docs/findings/ethA_deletion_ERR8975061.md). Genes
  not examined for coverage could carry undetected deletions.
