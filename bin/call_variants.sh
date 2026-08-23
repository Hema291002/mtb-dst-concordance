#!/usr/bin/env bash
set -uo pipefail

REF="$HOME/mtb-data/ref/GCF_000195955.2_ASM19595v2_genomic.fna"
BAM="$HOME/mtb-data/bam_markdup"
VCF="$HOME/mtb-data/vcf"

mkdir -p "$VCF"

for b in "$BAM"/*.markdup.bam; do
    run=$(basename "$b" .markdup.bam)
    out="$VCF/${run}.raw.vcf.gz"

    if [[ -s "${out}.csi" ]]; then
        echo "SKIP  $run"
        continue
    fi

    echo "CALL  $run"
    bcftools mpileup -f "$REF" -a AD,DP,SP -q 20 -Q 20 -d 500 -Ou "$b" \
      | bcftools call -mv --ploidy 1 -Oz -o "$out"
    bcftools index "$out"
done

echo "done: $(ls -1 "$VCF"/*.raw.vcf.gz.csi 2>/dev/null | wc -l) / 30"
