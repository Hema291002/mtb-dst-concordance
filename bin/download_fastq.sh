#!/usr/bin/env bash
set -uo pipefail

MANIFEST="$HOME/projects/mtb-dst-concordance/assets/download_manifest.tsv"
DEST="$HOME/mtb-data/raw"
JOBS=4

mkdir -p "$DEST"
cd "$DEST"

fetch_one() {
    local filename="$1" url="$2" md5="$3"

    if [[ -f "$filename" ]]; then
        if echo "$md5  $filename" | md5sum -c --status -; then
            echo "OK (cached)   $filename"
            return 0
        fi
        echo "RESUMING      $filename"
    fi

    curl -sL -C - --retry 5 --retry-delay 10 --retry-all-errors -o "$filename" "$url"

    if echo "$md5  $filename" | md5sum -c --status -; then
        echo "OK            $filename"
    else
        echo "CHECKSUM FAIL $filename"
        return 1
    fi
}
export -f fetch_one

tail -n +2 "$MANIFEST" \
  | awk -F'\t' '{print $2, $3, $4}' \
  | xargs -P "$JOBS" -n 3 bash -c 'fetch_one "$0" "$1" "$2"'

echo
echo "--- summary ---"
echo "files present: $(ls -1 *.fastq.gz 2>/dev/null | wc -l) / 60"
