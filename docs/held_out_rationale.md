# Held-out tables

PREDICTIONS.parquet and EFFECTS.parquet contain the consortium's own
genotypic resistance predictions and variant effect annotations.

These were moved to data/meta/held_out/ before any analysis began, and
were not opened during pipeline development.

Reason: this project generates independent variant calls from raw reads
and compares them against laboratory DST phenotypes. Inspecting the
consortium's predictions beforehand would allow filtering thresholds to
be tuned, consciously or otherwise, until they reproduced a known
answer. Held-out data is used only at the final comparison stage.
