# Gene naming differs between our annotation and the WHO catalogue

Position 1673425 C>T is a well-established isoniazid resistance mutation.

- Our annotator: `fabG1 c.-15C>T` (15 bases upstream of fabG1, start 1673440)
- WHO catalogue: `inhA_c.-777C>T` (777 bases upstream of inhA, start 1674202)

Both are correct. fabG1 and inhA lie 18 bp apart in a single operon, and the
WHO numbers promoter variants relative to inhA. The string `fabG1` does not
appear anywhere in the catalogue master file.

## Why this did not cause an error

Matching is performed on CHROM, POS, REF and ALT as the catalogue authors
recommend, not on variant names. The naming difference is therefore
irrelevant to the lookup.

Had names been matched instead, this variant would have failed to match in
3 of 30 isolates, all of them isoniazid-resistant, and would have been
scored as "no resistance mutation found". Three false negatives out of 15
resistant isolates would have reduced apparent sensitivity by 20 percentage
points with no error message and no indication in the output.

## Tier 1 coverage

| Drug | WHO tier 1 genes | Annotated here |
|---|---|---|
| Isoniazid | katG, inhA, ahpC | yes (fabG1 region included, catalogued as inhA upstream) |
| Rifampicin | rpoB | yes |

Tier 2 genes were not examined. A tier 2 variant could explain a
phenotypically resistant isolate carrying no tier 1 mutation, and this
possibility must be listed when discussing any such discordance.
