# Large deletion removing ethA and ethR in ERR8975061

**Status:** observed, not yet cross-checked against phenotype.
**Date observed:** 2026-08-21, during post-alignment coverage QC.

## Observation

Zero read coverage across NC_000962.3:4325429-4328703 (3,275 bp) in
ERR8975061, an isolate with 103.0x mean genome-wide depth.

Control: ERR181948, mean genome-wide depth 60.7x, has 70.3x mean depth
over the same interval. The absence is specific to ERR8975061, not a
property of the locus.

## Genes affected

| Locus | Gene | Coordinates | Extent of loss |
|---|---|---|---|
| Rv3852 | hns | 4325074-4325478 | partial (3' end) |
| Rv3853 | rraA | 4325495-4325968 | complete |
| Rv3854c | ethA | 4326004-4327473 | complete |
| Rv3855 | ethR | 4327549-4328199 | complete |
| Rv3856c | - | 4328401-4329408 | partial (3' end) |

## Interpretation

ethA encodes the monooxygenase that activates the prodrug ethionamide.
Loss-of-function in ethA is a documented mechanism of ethionamide
resistance. Complete deletion of the gene would be expected to abolish
activation entirely.

ethR, the transcriptional repressor of ethA, is also fully deleted.

## What this observation cannot establish

- Whether the isolate is phenotypically ethionamide-resistant. CRyPTIC
  records ETH phenotypes; this has not yet been checked. See registered
  prediction below.
- Whether loss of rraA, hns or Rv3856c has any phenotypic consequence.
- The mechanism generating the deletion. The involvement of flanking
  genes is consistent with a large recombination event, but this
  analysis provides no direct evidence.
- Whether similar deletions exist elsewhere in this or other isolates.
  Only 12 resistance-associated regions were examined for coverage.

## Registered prediction

Recorded before consulting the phenotype data: ERR8975061 is expected to
be phenotypically resistant to ethionamide (ETH). If the recorded
phenotype is susceptible, that discordance is itself a finding and must
be reported rather than explained away.

## Methodological significance

No variant caller in this pipeline emits a record for a deletion of this
size. A VCF describes differences at positions present in the data; it
has no representation for absent sequence. This deletion was detectable
only through coverage analysis.

Consequence: genotypic resistance prediction based on variant calls has a
structural blind spot for gene-scale deletions. If a comparable deletion
affected katG or rpoB in any isolate, the prediction would return
wild-type and the isolate would be scored as a false negative with no
indication of why.

## Commands used

    samtools depth -a -r NC_000962.3:4310000-4345000 \
      bam_markdup/ERR8975061.markdup.bam | awk '...'   # breakpoints
    samtools depth -a -r NC_000962.3:4325500-4328200 \
      bam_markdup/ERR181948.markdup.bam                # control
## Outcome (checked 2026-08-24)

The registered prediction is confirmed.

| Source | Result | Quality |
|---|---|---|
| DST_MEASUREMENTS, liquid media, CRyPTIC | ETH resistant | HIGH |
| UKMYC plate assay | ETH MIC >8, resistant | HIGH |

Two independent measurements agree.

### Why the deletion is the parsimonious explanation

The isolate is susceptible to 10 of the 14 drugs tested (AMI, BDQ, CFZ, DLM,
KAN, LEV, LZD, MXF, RFB, RIF), so the ethionamide resistance is specific
rather than part of a broadly resistant profile.

Isoniazid and ethionamide share a target, InhA, so inhA promoter mutations
confer cross-resistance to both. This isolate carries no inhA or fabG1
variant. Its only non-synonymous variant in the five annotated genes is
katG S315T, which blocks isoniazid activation by KatG. KatG plays no part in
activating ethionamide, which is activated by EthA. The isoniazid mechanism
therefore cannot account for the ethionamide phenotype.

### What this does not establish

- Causation. One isolate with one deletion is consistent with the mechanism,
  not a demonstration of it. ethA coverage was not examined in the other 29
  isolates, so no comparison group exists.
- The isolate is also resistant to ethambutol (EMB, MIC 8.0). embB was
  included in the coverage targets but not annotated, since this analysis was
  scoped to isoniazid and rifampicin. That resistance is unexplained here.

### Significance

The deletion was found by coverage analysis alone. No variant caller in this
pipeline emits a record for absent sequence. A genotypic prediction based on
variant calls would have reported this isolate as carrying no ethionamide
resistance mechanism, and would have been wrong.
