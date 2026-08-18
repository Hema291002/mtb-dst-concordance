# Analysis status

Honest record of what has actually been run versus what is planned.
Updated as work proceeds. Nothing is marked complete until its output
has been inspected.

| Stage | Status | Notes |
|---|---|---|
| 1. Environment setup | complete | WSL2, ext4, conda env `mtb-meta` |
| 2. Metadata acquisition | complete | CRyPTIC v3.4.0, Zenodo 15680920 |
| 3. Sample selection | complete | 30 isolates, seed 42, balanced 8/7/7/8 |
| 4. Reference genome | not started | |
| 5. FASTQ download | not started | |
| 6. Raw read QC | not started | |
| 7. Trimming | not started | |
| 8. Alignment | not started | |
| 9. Post-alignment QC | not started | |
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
