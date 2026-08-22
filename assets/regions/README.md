# Region definitions

## targets.bed
12 resistance-associated genes with 200 bp upstream flank (strand-aware).
The flank is required: the commonest non-katG isoniazid resistance mutation
is a promoter variant upstream of fabG1 (position 1673425), which lies
outside the gene body.

Note: fabG1 and inhA regions overlap by 182 bp because they are adjacent in
one operon. A variant in the overlap will match both regions.

## mask.bed
297 merged regions, 364,536 bp, 8.26% of the genome. Comprises:
- PE and PPE family genes (154 gene features, includes PE_PGRS)
- transposase and insertion-sequence CDS features (37)
- mobile_genetic_element features, including 16 IS6110 copies
- repeat_region features, including MIRUs and short direct repeats

Rationale: these sequences have high internal similarity, so short reads
map ambiguously and generate variant calls that reflect mapping failure
rather than genuine sequence difference.

Limitation: this is a name-based proxy for mappability, not a computed
mappability track. It is the field-standard approach but not optimal. Some
short repeats (39-53 bp) are masked despite being spanned by 150 bp reads,
so the mask errs toward over-exclusion. Verified to have no overlap with
targets.bed (positive control applied to confirm the check functions).
