# Why annotation was done with custom code

bcftools csq was tried first. It ran without error and annotated nothing,
reporting "Indexed 3 transcripts, 3 exons, 0 CDSs, 0 UTRs" against a GFF
containing 3,906 CDS features.

Cause: csq expects the Ensembl GFF3 chain gene -> transcript -> CDS. NCBI
RefSeq bacterial annotation links CDS directly to gene
(Parent=gene-Rv0001), because bacteria have no introns and no separate mRNA
feature is emitted. csq found no transcripts to attach CDSs to and indexed
none. This is a silent format incompatibility, not a bug in either tool.

A custom annotator was written for the five genes relevant to isoniazid and
rifampicin, rather than converting the GFF or fighting the packaging.

## Validation

bin/annotate_targets.py --selftest checks four independently published facts
before the annotator is used:

| Check | Tests |
|---|---|
| rpoB 761155 C>T -> rpoB S450L | forward-strand codon arithmetic |
| katG 2155168 C>G -> katG S315T | reverse-strand codon arithmetic |
| fabG1 1673425 C>T -> fabG1 c.-15C>T | promoter offset handling |
| katG codon 315 translates to serine | codon table and reverse complement |

During validation the reference-base check caught an error in a test case:
reverse-strand bases had been specified where VCF convention requires
forward-strand bases. The code was correct; the test input was not. Any
variant whose REF field disagrees with the reference genome is reported as
REF_MISMATCH rather than annotated.

## Scope and limitations

- Only five genes are annotated. Variants elsewhere are not interpreted.
- Indels are reported and classified as frameshift or in-frame but are not
  translated. Frameshift consequences depend on downstream sequence and
  short-read indel calls are less reliable than SNV calls.
- A snpEff cross-check on the same five genes is planned as an independent
  confirmation of the translation logic.
