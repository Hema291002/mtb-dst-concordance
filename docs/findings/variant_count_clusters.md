# Raw variant counts cluster into four groups

**Status:** observed before lineage assignment. Prediction registered.
**Date:** 2026-08-21

Raw (unfiltered) variant counts against H37Rv, all 30 isolates:

| Group | Range | n |
|---|---|---|
| A | 724-783 | 4 |
| B | 1001-1246 | 13 |
| C | 1705-1849 | 11 |
| D | 2499-2517 | 2 |

Gaps between groups (218, 459, 650) are large relative to within-group spread.

## Registered prediction

These groups correspond to distinct M. tuberculosis lineages, ordered by
evolutionary distance from H37Rv (lineage 4). Specifically:
- Group A is expected to be lineage 4, the same lineage as the reference.
- Groups B, C and D are expected to be progressively more distant lineages.
- Isolates within a group are expected to share a lineage assignment.

To be tested at the lineage assignment stage using the Coll 62-SNP barcode,
computed from these same VCFs. The GENOMES table's LINEAGE column remains
held out until then.

## What this observation cannot establish

- Which lineage each group is. Distance from H37Rv orders them but does not
  name them.
- Whether group D represents a distinct lineage or something else entirely.
- Anything causal. Variant count is a consequence of evolutionary history,
  not of resistance.

## Caveat

These are raw, unfiltered counts and include false positives from repetitive
regions. The clustering is if anything more striking for surviving that noise,
but the numbers will change after filtering.
