#!/usr/bin/env bash
set -uo pipefail

VCF="$HOME/mtb-data/vcf"
MASK="$HOME/mtb-data/ref/mask.bed"
OUT="$HOME/mtb-data/vcf_filtered"

EXPR='QUAL>=30 && FORMAT/DP>=10 && MQ>=40 && (FORMAT/AD[0:1])/(FORMAT/AD[0:0]+FORMAT/AD[0:1])>=0.9'

mkdir -p "$OUT"
printf "sample\traw\tfiltered\tpct_kept\n" > "$OUT/filter_stats.tsv"

for v in "$VCF"/*.raw.vcf.gz; do
    run=$(basename "$v" .raw.vcf.gz)
    out="$OUT/${run}.filt.vcf.gz"

    if [[ ! -s "${out}.csi" ]]; then
        echo "FILTER $run"
        bcftools view -T "^${MASK}" -i "$EXPR" -Oz -o "$out" "$v" 2>/dev/null
        bcftools index "$out"
    fi

    raw=$(bcftools view -H "$v" 2>/dev/null | wc -l)
    filt=$(bcftools view -H "$out" 2>/dev/null | wc -l)
    awk -v s="$run" -v r="$raw" -v f="$filt" 'BEGIN{printf "%s\t%d\t%d\t%.1f\n", s, r, f, 100*f/r}' >> "$OUT/filter_stats.tsv"
done

echo "done: $(ls -1 "$OUT"/*.filt.vcf.gz.csi 2>/dev/null | wc -l) / 30"
