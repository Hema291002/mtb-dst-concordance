# Predictions registered before consulting the WHO catalogue

Date: 2026-08-22. Recorded before the catalogue was downloaded or opened.

Variants observed in our 30 isolates, with predicted WHO grading:

| Variant | Isolates | Predicted grade |
|---|---|---|
| katG S315T | 11 | 1, Associated with resistance |
| rpoB S450L | 13 | 1, Associated with resistance |
| rpoB H445Y | 2 | 1 or 2, Associated with resistance |
| fabG1 c.-15C>T | 3 | 1, Associated with resistance |
| katG R463L | 13 | 4 or 5, NOT associated with resistance (phylogenetic marker) |
| ahpC c.-88G>A | 5 | 4 or 5, expected to co-segregate with lineage |
| rpoB c.-61C>T | 5 | 4 or 5, expected to co-segregate with lineage |

Rationale for katG R463L: it appears in more isolates than S315T, spans
both resistant and susceptible samples, and co-occurs with other variants
in a pattern consistent with shared ancestry rather than drug selection.

Variants with no prediction (insufficient prior knowledge): katG A109V,
katG Q295K, katG K557N, rpoB A286V, rpoB G973D, rpoB V695L, rpoB P45S,
inhA I21T, inhA I21V, ahpC c.-142G>A, ahpC c.-48G>A, fabG1 c.-34C>G.

If any prediction above is wrong, that is recorded rather than revised.
## Outcome (checked 2026-08-23)

| Variant | Predicted | Actual (WHO v2, column 106) |
|---|---|---|
| katG S315T | 1) Assoc w R | **1) Assoc w R** |
| rpoB S450L | 1) Assoc w R | **1) Assoc w R** |
| katG R463L | 4 or 5, not associated | **5) Not assoc w R** |

All three confirmed. katG R463L received the strongest negative grading
available despite appearing in more isolates (13) than S315T (11), which is
the clearest available demonstration that presence in a resistance gene does
not imply a resistance mechanism.

## Annotator validation, revised

The custom annotator was to be cross-checked against snpEff. That is now
redundant: the WHO genomic-coordinate file is an independent derivation of
position-to-amino-acid mapping, produced by a different group with different
code, and it agreed exactly with our annotation on all three tested variants
including two on the reverse strand. The snpEff cross-check was dropped and
the reason recorded here rather than left as an unexplained gap.
