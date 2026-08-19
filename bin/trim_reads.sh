#!/usr/bin/env bash
set -uo pipefail

RAW="$HOME/mtb-data/raw"
OUT="$HOME/mtb-data/trimmed"
REP="$HOME/mtb-data/results/qc/fastp"
THREADS=4

mkdir -p "$OUT" "$REP"

for r1 in "$RAW"/*_1.fastq.gz; do
    run=$(basename "$r1" _1.fastq.gz)
    r2="$RAW/${run}_2.fastq.gz"

    if [[ -s "$OUT/${run}_1.trim.fastq.gz" && -s "$OUT/${run}_2.trim.fastq.gz" ]]; then
        echo "SKIP  $run"
        continue
    fi

    echo "TRIM  $run"
    fastp \
        -i "$r1" -I "$r2" \
        -o "$OUT/${run}_1.trim.fastq.gz" \
        -O "$OUT/${run}_2.trim.fastq.gz" \
        --detect_adapter_for_pe \
        --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 \
        --length_required 50 \
        --n_base_limit 5 \
        --thread "$THREADS" \
        --json "$REP/${run}.fastp.json" \
        --html "$REP/${run}.fastp.html" \
        2> "$REP/${run}.fastp.log"
done

echo "done: $(ls -1 "$OUT"/*_1.trim.fastq.gz 2>/dev/null | wc -l) / 30"
