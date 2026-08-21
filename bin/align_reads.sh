#!/usr/bin/env bash
set -uo pipefail

REF="$HOME/mtb-data/ref/GCF_000195955.2_ASM19595v2_genomic.fna"
TRIM="$HOME/mtb-data/trimmed"
BAM="$HOME/mtb-data/bam"

mkdir -p "$BAM"

for r1 in "$TRIM"/*_1.trim.fastq.gz; do
    run=$(basename "$r1" _1.trim.fastq.gz)
    r2="$TRIM/${run}_2.trim.fastq.gz"
    out="$BAM/${run}.sorted.bam"

    if [[ -s "${out}.bai" ]]; then
        echo "SKIP  $run"
        continue
    fi

    echo "ALIGN $run"
    bwa-mem2 mem -t 4 \
        -R "@RG\tID:${run}\tSM:${run}\tLB:${run}\tPL:ILLUMINA" \
        "$REF" "$r1" "$r2" 2> "$BAM/${run}.bwa.log" \
      | samtools sort -@ 2 -m 512M -o "$out" - 2>> "$BAM/${run}.bwa.log"

    samtools index "$out"
done

echo "done: $(ls -1 "$BAM"/*.sorted.bam.bai 2>/dev/null | wc -l) / 30"
