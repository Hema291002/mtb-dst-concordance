# Variant filter thresholds and how they were chosen

Thresholds were chosen by inspecting the distributions in ERR181948 before
any were applied, not imported from another study.

| Filter | Threshold | Reason |
|---|---|---|
| QUAL | >= 30 | Scores saturate at 225; 914/1065 sit above 200. No natural break exists higher up, so this cuts only the indefensible tail. |
| FORMAT/DP | >= 10 | Raw minimum depth was 1. Below ~10 reads, allele fraction cannot be meaningfully estimated. |
| MQ | >= 40 | Read-level counterpart to the region-level mask; excludes ambiguously placed reads. |
| Allele fraction | >= 0.9 | M. tuberculosis is haploid and these are pure cultures. Over 75% of raw variants sit at exactly 1.0. Only 5.5% fall below 0.9. |
| Repeat mask | excluded | Reads from repetitive families map ambiguously and generate false calls. |

## Independent corroboration of the mask

Allele fraction and the repeat mask were derived independently: one from the
read data, one from the genome annotation. In ERR181948, 31 of 62 variants
with allele fraction below 0.9 fall inside masked regions. The mask covers
8.26% of the genome, so ~5 would be expected by chance. This six-fold
enrichment indicates the low-allele-fraction calls are largely caused by
ambiguous read placement, which is exactly what the mask predicts.

## Effect of each filter (ERR181948, 1065 raw variants)

| Step | Remaining | Removed |
|---|---|---|
| raw | 1065 | - |
| QUAL >= 30 | 1021 | 44 |
| + DP >= 10 | 962 | 59 |
| + MQ >= 40 | 943 | 19 |
| + AF >= 0.9 | 898 | 45 |
| + outside mask | 768 | 130 |

768 filtered SNPs is consistent with published divergence of a lineage 4
isolate from H37Rv (roughly 700-1200 SNPs), which is an external check rather
than an internal consistency check.

## Note on filtering and rounding

Filters operate on unrounded values. An exploratory command that rounded
allele fraction to three decimals before comparison produced a slightly
different count. Values are formatted for display only, never before a
comparison.
