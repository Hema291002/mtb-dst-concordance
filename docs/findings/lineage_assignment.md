# Lineage assignment and the variant-count clustering prediction

## Method

TB-Profiler SNP barcode (jodyphelan/tbdb), 1,111 positions, 126 labels.

Note on provenance: the widely cited Coll et al. 2014 62-SNP barcode was
superseded by Napier et al. 2020, which extended it to 90 SNPs using 35,298
isolates. The file used here is larger still, covering deep sublineages and
animal-adapted species.

Barcode coordinates are H37Rv but the sequence is named "Chromosome" rather
than NC_000962.3. Matching is therefore on position and allele only. This is
safe because both sides have a single chromosome, and is stated rather than
left implicit.

Coordinates were verified against our own reference before use: for 20 sampled
barcode rows, the reference base differed from the lineage-defining allele, as
expected since H37Rv carries the ancestral state at markers defining other
lineages.

Raw VCFs were used rather than mask-filtered ones. 13 of 1,111 barcode
positions fall inside the repeat mask, and lineage markers should not be lost
for reasons unrelated to lineage. Depth >= 10 and allele fraction >= 0.9 were
applied directly instead.

Assignment requires a complete ancestral path: an isolate assigned
lineage4.2.2.1 must also carry markers for lineage4.2.2, lineage4.2 and
lineage4. Off-path markers are reported, not discarded.

## Result

| Lineage | n |
|---|---|
| 4 | 16 |
| 2 | 6 |
| 3 | 5 |
| 1 | 2 |
| unassigned | 1 |

No isolate carried off-path markers, so no evidence of mixed infection.

## Registered prediction: outcome

The prediction (docs/findings/variant_count_clusters.md) was that the four
variant-count groups correspond to lineages ordered by distance from H37Rv,
with group A being lineage 4.

| Group | Filtered variants | n | Lineages found |
|---|---|---|---|
| A | 516-537 | 4 | lineage 4.7, 4.8 |
| B | 712-890 | 13 | lineage 4.1-4.4, plus 1 unassigned |
| C | 1253-1386 | 11 | lineage 2 (6), lineage 3 (5) |
| D | 1922-1935 | 2 | lineage 1.1.2 |

CONFIRMED: the ordering tracks phylogenetic distance. Lineage 4 nearest the
reference, lineages 2 and 3 intermediate, lineage 1 furthest, consistent with
lineage 1 having diverged earliest among these four.

CONFIRMED, and more finely than predicted: group A is lineage 4.7 and 4.8,
the sublineages nearest H37Rv (itself 4.9). Group B is lineage 4.1-4.4. The
A/B split is sublineage distance within lineage 4.

FAILED: the sub-claim that each group contains a single lineage. Group C
mixes lineages 2 and 3. This is biologically expected rather than an error:
lineages 2 and 3 are roughly equidistant from lineage 4, so variant count
resolves depth but not direction.

## Independent corroboration from our own data

Three patterns were observed in the variant tables before any lineage
information was available, and all three resolve as lineage markers:

- katG R463L (WHO grade 5, not associated with resistance) is present in
  exactly the 13 non-lineage-4 isolates and absent from all 16 lineage-4
  isolates.
- rpoB c.-61C>T and ahpC c.-88G>A co-occur in exactly 5 isolates. All 5 are
  lineage 3.
- ahpC c.-142G>A, ungraded by the WHO catalogue, occurs in exactly 2
  isolates. Both are lineage 1.1.2.

The barcode named patterns that had already been identified from the data.

## Why one isolate is unassigned

ERR4813681 carries variants at 11 barcode positions but none with a
lineage-defining allele.

The barcode contains lineage4.9 markers, labelled "Euro-American (H37Rv-like)".
At five sampled lineage4.9 positions the reference base equals the defining
allele. H37Rv is itself lineage 4.9, so a genuine 4.9 isolate matches the
reference at those positions and produces no variant call.

Lineage 4.9 is therefore undetectable when assigning lineage from a
variants-only VCF. ERR4813681 is most plausibly lineage 4.9; the absence of
markers is informative rather than missing data.

This is a limitation of the approach used here, not of the barcode.
TB-Profiler works from alignments and handles reference-matching positions
correctly. The script reports the count of such positions (20 of 1,111) so the
limitation is visible.

## Bearing on the main analysis

Resistance category is partly confounded with collecting site, and site
predicts lineage (documented at sample selection). The two robust isoniazid
discordances are ERR4814489 (lineage 1.1.2) and ERR8975559 (lineage 2.2.1),
both non-lineage-4. With two isolates this is a description, not a test, and
no association is claimed.
