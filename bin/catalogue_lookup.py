#!/usr/bin/env python3
"""
Match annotated target-gene variants against the WHO mutation catalogue
(2nd edition, 2023) and produce a genotypic resistance prediction per isolate.

Matching follows the catalogue authors' own instructions: exact match on
CHROM, POS, REF and ALT against their genomic-coordinates VCF, then lookup of
the resulting graded-variant name in the master file. Variant names are never
matched directly, because naming conventions differ (our fabG1 c.-15C>T is
the catalogue's inhA_c.-777C>T for the same nucleotide).

This script does NOT read phenotype data. Prediction and truth are joined in a
separate step so the prediction cannot be influenced by the answer.

Outputs two files:
  catalogue_matches.tsv      one row per annotated variant, with gradings
  genotypic_predictions.tsv  one row per isolate per drug

Usage:
    python3 bin/catalogue_lookup.py \
        --variants docs/target_variants.tsv \
        --vcf   ~/mtb-data/catalogue/Genomic_coordinates_7May2024.vcf.gz \
        --master ~/mtb-data/catalogue/WHO-UCN-TB-2023.6-eng_catalogue_master_file.txt \
        --outdir ~/mtb-data/results/catalogue
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict

# Drugs we are predicting. Names must match the catalogue's 'drug' column.
DRUGS = ["Isoniazid", "Rifampicin"]

# WHO confidence gradings that count as predicting resistance.
#
# The catalogue grades variants 1 to 5. Groups 1 and 2 are the two levels of
# "associated with resistance"; group 3 is uncertain; groups 4 and 5 are the
# two levels of "not associated". Treating 1 and 2 as resistant and everything
# else as not is the standard interpretation, but it IS a choice, and moving
# group 3 into the resistant set would change sensitivity and specificity.
# The choice is made here, in one place, and stated in the output header.
RESISTANT_GRADES = ("1)", "2)")


def load_catalogue_vcf(path):
    """Map (pos, ref, alt) -> list of graded-variant names.

    Each genomic variant may map to several graded-variant names, concatenated
    with '&' in the INFO field, per the catalogue's documented format.
    """
    index = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 8:
                continue
            pos, ref, alt, info = int(f[1]), f[3], f[4], f[7]
            names = []
            for part in info.split(";"):
                if part.startswith("graded_variant="):
                    names = part.split("=", 1)[1].split("&")
                    break
            if names:
                index[(pos, ref, alt)] = names
    return index


def load_master(path, drugs):
    """Map (drug, graded-variant name) -> final confidence grading.

    The grading column is located by name rather than by index, so a change in
    column order upstream produces a clear error instead of silently wrong
    gradings.
    """
    grades = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            i_drug = header.index("drug")
            i_var = header.index("variant")
            i_tier = header.index("tier")
            i_grade = header.index("FINAL CONFIDENCE GRADING")
        except ValueError as e:
            sys.exit(f"ERROR: expected column not found in master file: {e}")

        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= i_grade:
                continue
            if f[i_drug] not in drugs:
                continue
            grades[(f[i_drug], f[i_var])] = (f[i_grade], f[i_tier])
    return grades, i_grade + 1


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--variants", required=True,
                    help="annotated target variants TSV")
    ap.add_argument("--vcf", required=True,
                    help="catalogue genomic coordinates VCF")
    ap.add_argument("--master", required=True,
                    help="catalogue master file TSV")
    ap.add_argument("--samples", default=None,
                    help="file listing all expected sample IDs, one per line")
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    outdir = os.path.expanduser(args.outdir)
    os.makedirs(outdir, exist_ok=True)

    sys.stderr.write("loading catalogue coordinates ...\n")
    coords = load_catalogue_vcf(os.path.expanduser(args.vcf))
    sys.stderr.write(f"  {len(coords)} genomic variants indexed\n")

    sys.stderr.write("loading catalogue gradings ...\n")
    grades, grade_col = load_master(os.path.expanduser(args.master), DRUGS)
    sys.stderr.write(f"  {len(grades)} drug/variant gradings loaded "
                     f"(grading read from column {grade_col})\n")

    # ---------------------------------------------------------------
    # Per-variant matching
    # ---------------------------------------------------------------
    matched = 0
    unmatched = 0
    rows = []

    with open(os.path.expanduser(args.variants)) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < len(header):
                continue
            sample = f[idx["sample"]]
            pos = int(f[idx["pos"]])
            ref = f[idx["ref"]]
            alt = f[idx["alt"]]
            gene = f[idx["gene"]]
            notation = f[idx["notation"]]
            effect = f[idx["effect"]]

            names = coords.get((pos, ref, alt))
            if names is None:
                unmatched += 1
                rows.append([sample, str(pos), ref, alt, gene, notation,
                             effect, "NOT_IN_CATALOGUE", "", "", ""])
                continue

            matched += 1
            for name in names:
                inh = grades.get(("Isoniazid", name), ("", ""))
                rif = grades.get(("Rifampicin", name), ("", ""))
                rows.append([sample, str(pos), ref, alt, gene, notation,
                             effect, name, inh[0], rif[0],
                             inh[1] or rif[1]])

    match_path = os.path.join(outdir, "catalogue_matches.tsv")
    with open(match_path, "w") as out:
        out.write("\t".join(["sample", "pos", "ref", "alt", "gene",
                             "our_notation", "effect", "who_variant",
                             "INH_grade", "RIF_grade", "tier"]) + "\n")
        for r in rows:
            out.write("\t".join(r) + "\n")

    sys.stderr.write(f"  {matched} variants matched, "
                     f"{unmatched} not in catalogue\n")

    # ---------------------------------------------------------------
    # Per-isolate prediction
    # ---------------------------------------------------------------
    # A sample is predicted resistant to a drug if it carries at least one
    # variant graded 1 or 2 for that drug.

    samples = set()
    support = defaultdict(list)          # (sample, drug) -> [variant names]

    for r in rows:
        sample, _, _, _, _, notation, _, who_var, inh, rif, _ = r
        samples.add(sample)
        if inh.startswith(RESISTANT_GRADES):
            support[(sample, "INH")].append(f"{who_var} [{inh}]")
        if rif.startswith(RESISTANT_GRADES):
            support[(sample, "RIF")].append(f"{who_var} [{rif}]")

    # Every sample in the annotation file is represented, including samples
    # with no graded variant at all. A sample missing from the output would be
    # indistinguishable from a sample predicted susceptible, and the two mean
    # very different things.
    if args.samples:
        with open(os.path.expanduser(args.samples)) as fh:
            expected = {l.strip() for l in fh if l.strip()}
        missing = expected - samples
        if missing:
            sys.stderr.write(
                f"note: {len(missing)} sample(s) carry no variant in any "
                f"target gene and are recorded as predicted susceptible: "
                f"{', '.join(sorted(missing))}\n")
        samples |= expected
    all_samples = sorted(samples)

    pred_path = os.path.join(outdir, "genotypic_predictions.tsv")
    with open(pred_path, "w") as out:
        out.write("# genotypic resistance prediction from WHO catalogue v2 "
                  "(2023)\n")
        out.write("# predicted resistant if any variant graded "
                  "1) Assoc w R or 2) Assoc w R - Interim\n")
        out.write("# phenotype data was not read by this script\n")
        out.write("\t".join(["sample", "drug", "predicted",
                             "n_variants", "variants"]) + "\n")
        for s in all_samples:
            for drug in ("INH", "RIF"):
                vs = support.get((s, drug), [])
                out.write("\t".join([
                    s, drug,
                    "R" if vs else "S",
                    str(len(vs)),
                    ";".join(sorted(set(vs))) if vs else "-"]) + "\n")

    sys.stderr.write(f"\nwrote {match_path}\n")
    sys.stderr.write(f"wrote {pred_path}\n")
    sys.stderr.write(f"{len(all_samples)} samples represented\n")

    if len(all_samples) < 30:
        sys.stderr.write(
            f"\nWARNING: only {len(all_samples)} samples appear in the "
            "annotation file. Samples with no variant in any target gene "
            "produce no rows and will be missing here. Check against the "
            "full sample list before interpreting.\n")


if __name__ == "__main__":
    main()
