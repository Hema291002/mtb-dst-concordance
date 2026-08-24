# Two isoniazid-resistant isolates with no catalogued mechanism

**Status:** observed, five explanations eliminated, three remain untestable.

## The observation

ERR4814489 and ERR8975559 are phenotypically isoniazid-resistant and carry no
WHO-graded resistance variant in any tier 1 gene. Both are discordant under
both phenotype labellings (primary and sensitivity analysis), so neither is an
artefact of the phenotype resolution rule.

Every variant either isolate carries in the five target genes:

| Isolate | Variant | INH grade | RIF grade |
|---|---|---|---|
| ERR4814489 | rpoB V228V (synonymous) | - | 5) Not assoc |
| ERR4814489 | rpoB A1075A (synonymous) | - | 5) Not assoc |
| ERR4814489 | katG R463L | 5) Not assoc | - |
| ERR4814489 | ahpC c.-142G>A | not graded | not graded |
| ERR8975559 | rpoB A1075A (synonymous) | - | 5) Not assoc |
| ERR8975559 | inhA G3G (synonymous) | 5) Not assoc | - |
| ERR8975559 | katG R463L | 5) Not assoc | - |

Nothing here can plausibly cause isoniazid resistance. katG R463L is a
phylogenetic marker present in every non-lineage-4 isolate in this set.
ahpC c.-142G>A is present only in the two lineage 1 isolates and is likewise
lineage-associated.

## Strength of the phenotype

ERR8975559: UKMYC plate assay, HIGH quality, MIC >12.8 against a critical
concentration of 0.1. More than 128-fold above the cutoff. High-level
resistance, not a borderline call.

ERR4814489: solid-media DST, WHO critical concentration, MEDIUM quality,
resistant. No MIC recorded.

## Explanations eliminated, with the evidence

| Explanation | Status | Evidence |
|---|---|---|
| Coverage gap in a resistance gene | ruled out | 100% of katG, inhA, ahpC, fabG1 covered at 30x in both |
| Gene-scale deletion | ruled out | same coverage analysis |
| Low-frequency subpopulation (heteroresistance) | ruled out | no variant below 0.9 allele fraction in any target region, in raw pre-filter calls |
| Variant lost to filtering | ruled out | raw VCFs inspected directly |
| Variant graded 3, uncertain significance | ruled out | neither isolate carries any group 3 variant |
| Phenotype error | implausible for ERR8975559 | MIC 128x above the critical concentration, HIGH quality |

## Explanations remaining

- A variant in a WHO tier 2 gene. Not examined here; this analysis was scoped
  to tier 1.
- An uncatalogued resistance mechanism.
- Non-genetic resistance, for example efflux-mediated.

None is testable with the present data.

## Lineage context

ERR4814489 is lineage 1.1.2 and ERR8975559 is lineage 2.2.1. Both are
non-lineage-4. With two isolates this is a description, not a test, and no
association between lineage and unexplained resistance is claimed.

## Why this is reported rather than dropped

Genotypic resistance prediction is being adopted clinically. Its failure modes
matter more than its success rate, and an isolate with high-level isoniazid
resistance and no detectable mechanism in any established gene is exactly the
case where a genome-based result would mislead a clinician.
