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
