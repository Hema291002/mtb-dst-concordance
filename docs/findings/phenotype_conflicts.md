# Conflicting phenotype records and how they were resolved

## The problem

4 of 60 isolate/drug pairs had two DST records that disagreed. The original
selection used `drop_duplicates(keep="first")`, which selects by file order.
That is arbitrary, and for these four pairs it determined the ground truth.

All four conflicts follow the same pattern:

| | Record A | Record B |
|---|---|---|
| Source | CLIRES (2020/2023) | CRyPTIC |
| Quality | MEDIUM | HIGH |
| Method | mostly solid media | liquid media |
| MIC recorded | no | yes |

All three affected isolates are from site.05, consistent with a method effect
rather than random disagreement. Solid-media and liquid-media DST are known to
disagree near the critical concentration.

## Rule adopted

1. Prefer PHENOTYPE_QUALITY == HIGH over MEDIUM
2. Among equals, prefer records carrying an MIC value

Both clauses are stated on quality grounds and were fixed before the effect on
results was computed.

## Independent corroboration

The UKMYC 96-well plate assay is a separate CRyPTIC measurement on a different
platform. It was not used to construct the rule. It agrees with every record
the rule selects:

| Isolate | UKMYC INH | UKMYC RIF |
|---|---|---|
| ERR4813445 | R, MIC 0.2 | S, MIC 0.06 |
| ERR4813681 | S, MIC 0.05 | S, MIC 0.06 |
| ERR3287554 | R, MIC 0.2 | R, MIC >4 |

## Sample selection was not repeated

The same 30 isolates are retained. Re-selecting under a different phenotype
rule would mean choosing the sample set after seeing results. Only the
phenotype label changed, for 4 of 60 pairs.

## Sensitivity analysis

Results under the original `keep="first"` labels are reported alongside the
primary analysis so the dependence on this choice is visible.

## Unresolved

ERR4813445 and ERR4813681 both carry rpoB S450L at allele fraction 1.000 with
70 and 80 supporting reads, the same call found in 11 concordant isolates. Two
independent measurements give a rifampicin MIC of 0.06, eightfold below the
0.5 critical concentration. rpoB S450L confers high-level rifampicin
resistance, so this combination should not occur.

Variant calling error is ruled out. The remaining possibilities are a genuine
and unusual biological phenomenon, or an error linking genome to phenotype.
The CRyPTIC release notes document that some UNIQUEIDs map to multiple ENA run
accessions with no way to determine which was used, which makes the second
plausible. This analysis cannot distinguish between them.
