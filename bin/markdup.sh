#!/usr/bin/env bash
set -uo pipefail

BAM="$HOME/mtb-data/bam"
OUT="$HOME/mtb-data/bam_markdup"
TMP="$HOME/mtb-data/tmp"

mkdir -p "$OUT" "$TMP"

for b in "$BAM"/*.sorted.bam; do
    run=$(basename "$b" .sorted.bam)
    out="$OUT/${run}.markdup.bam"

    if [[ -s "${out}.bai" ]]; then
        echo "SKIP  $run"
        continue
    fi

    echo "MARK  $run"
    samtools collate -@ 2 -O -u "$b" "$TMP/collate_${run}" \
      | samtools fixmate -@ 2 -m -u - - \
      | samtools sort -@ 2 -m 512M -u - \
      | samtools markdup -@ 2 -f "$OUT/${run}.markdup.stats" - "$out"

    samtools index "$out"
done

echo "done: $(ls -1 "$OUT"/*.markdup.bam.bai 2>/dev/null | wc -l) / 30"
