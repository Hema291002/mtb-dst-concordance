#!/usr/bin/env python3
"""
Assign M. tuberculosis lineage from the TB-Profiler SNP barcode.

The barcode (jodyphelan/tbdb, 1,111 positions) lists, for each position, the
allele that defines a given (sub)lineage. Labels are hierarchical:
lineage4.2.2.1 lies within lineage4.2.2, within lineage4.2, within lineage4.

Assignment therefore is not "find a marker, read its label". An isolate in a
deep sublineage should carry the defining marker at every level above it,
because each arose in an ancestor. We collect all markers present, then report
the deepest label whose full ancestral path is also present.

Markers found OFF that path are reported rather than discarded. They indicate
mixed infection, convergent mutation, or a calling error, and hiding them would
turn a signal into silence.

Raw VCFs are used, not mask-filtered ones: 13 of 1,111 barcode positions fall
inside the repeat mask, and lineage markers should not be lost for reasons
unrelated to lineage. Depth and allele-fraction thresholds are applied here
instead.

Usage:
    python3 bin/assign_lineage.py \
        --barcode ~/mtb-data/barcode/barcode.bed \
        --vcfdir  ~/mtb-data/vcf \
        --out     ~/mtb-data/results/lineage
"""

import argparse
import glob
import os
import subprocess
import sys

MIN_DP = 10
MIN_AF = 0.9


def load_barcode(path):
    """position -> {allele: lineage}. One position can define different
    lineages via different alleles, so the inner dict is keyed by allele."""
    bc = {}
    labels = set()
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 5:
                continue
            pos = int(f[2])          # BED end == 1-based position
            lineage, allele = f[3], f[4].upper()
            bc.setdefault(pos, {})[allele] = lineage
            labels.add(lineage)
    return bc, labels


def check_reference_collisions(bc, ref_path):
    """Warn if any barcode allele equals the reference base.

    At such a position every isolate would appear to carry the marker, since a
    VCF records only differences from the reference. Those positions cannot
    discriminate and must be known about rather than silently trusted.
    """
    if not os.path.exists(ref_path):
        return None
    seq = []
    with open(ref_path) as fh:
        for line in fh:
            if not line.startswith(">"):
                seq.append(line.strip())
    genome = "".join(seq).upper()
    collisions = 0
    for pos, alleles in bc.items():
        if pos - 1 < len(genome) and genome[pos - 1] in alleles:
            collisions += 1
    return collisions


def genotype_sample(vcf, bc):
    """Return the set of lineage labels whose defining allele this isolate
    carries, subject to depth and allele-fraction thresholds."""
    cmd = ["bcftools", "query", "-f",
           "%POS\t%REF\t%ALT\t[%DP]\t[%AD]\n", vcf]
    try:
        out = subprocess.run(cmd, capture_output=True, text=True,
                             check=True).stdout
    except subprocess.CalledProcessError as e:
        sys.exit(f"bcftools failed on {vcf}:\n{e.stderr}")

    found = {}
    for line in out.splitlines():
        f = line.split("\t")
        if len(f) < 5:
            continue
        pos = int(f[0])
        if pos not in bc:
            continue
        alt = f[2].upper()
        if alt not in bc[pos]:
            continue
        try:
            dp = int(f[3])
            ad = f[4].split(",")
            ad_ref, ad_alt = int(ad[0]), int(ad[1])
        except (ValueError, IndexError):
            continue
        if dp < MIN_DP:
            continue
        total = ad_ref + ad_alt
        if total == 0 or ad_alt / total < MIN_AF:
            continue
        found[bc[pos][alt]] = (pos, dp, ad_alt / total)
    return found


def prefixes(label):
    """Ancestral path of a hierarchical label, deepest last.
    lineage4.2.2.1 -> [lineage4, lineage4.2, lineage4.2.2, lineage4.2.2.1]"""
    parts = label.split(".")
    return [".".join(parts[:i + 1]) for i in range(len(parts))]


def resolve(found):
    """Deepest label with a complete ancestral path, plus off-path markers."""
    labels = set(found)
    if not labels:
        return None, [], []

    # Candidates ordered by depth, deepest first.
    ranked = sorted(labels, key=lambda l: (len(l.split(".")), l), reverse=True)

    best = None
    for cand in ranked:
        path = prefixes(cand)
        # A parent level may legitimately have no marker in the barcode at all,
        # so require only that every parent present in the barcode is matched.
        if all(p in labels or p not in ALL_LABELS for p in path):
            best = cand
            break

    if best is None:
        best = ranked[0]

    on_path = set(prefixes(best))
    off = sorted(l for l in labels if l not in on_path)
    return best, sorted(on_path & labels), off


def main():
    global ALL_LABELS
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--barcode", required=True)
    ap.add_argument("--vcfdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--ref", default=os.path.expanduser(
        "~/mtb-data/ref/GCF_000195955.2_ASM19595v2_genomic.fna"))
    args = ap.parse_args()

    bc, ALL_LABELS = load_barcode(os.path.expanduser(args.barcode))
    sys.stderr.write(f"barcode: {len(bc)} positions, "
                     f"{len(ALL_LABELS)} distinct lineage labels\n")

    col = check_reference_collisions(bc, os.path.expanduser(args.ref))
    if col is not None:
        sys.stderr.write(
            f"positions where a barcode allele equals the reference base: "
            f"{col}\n")
        if col:
            sys.stderr.write(
                "  these cannot discriminate, because a VCF records only "
                "differences from the reference\n")

    outdir = os.path.expanduser(args.out)
    os.makedirs(outdir, exist_ok=True)

    vcfs = sorted(glob.glob(os.path.join(
        os.path.expanduser(args.vcfdir), "*.raw.vcf.gz")))
    sys.stderr.write(f"{len(vcfs)} VCFs\n\n")

    rows = []
    for v in vcfs:
        sample = os.path.basename(v).replace(".raw.vcf.gz", "")
        found = genotype_sample(v, bc)
        best, path, off = resolve(found)
        rows.append((sample, best or "unassigned", len(found),
                     ";".join(path), ";".join(off)))
        print(f"{sample:<14} {best or 'unassigned':<20} "
              f"markers={len(found):<3} off_path={len(off)}")

    out = os.path.join(outdir, "lineage_assignments.tsv")
    with open(out, "w") as fh:
        fh.write("sample\tlineage\tn_markers\tancestral_path\t"
                 "off_path_markers\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")

    print(f"\nwrote {out}")

    # Summary by top-level lineage
    from collections import Counter
    tops = Counter(r[1].split(".")[0] for r in rows)
    print("\ntop-level lineage counts:")
    for k, n in sorted(tops.items()):
        print(f"  {k:<14} {n}")

    noff = sum(1 for r in rows if r[4])
    if noff:
        print(f"\n{noff} isolate(s) carry off-path markers. These are not "
              "errors to ignore:\nthey indicate mixed infection, convergent "
              "mutation, or miscalled positions.")


if __name__ == "__main__":
    main()
