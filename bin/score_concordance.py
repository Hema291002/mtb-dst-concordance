#!/usr/bin/env python3
"""
Score genotypic resistance predictions against laboratory phenotypes.

Primary analysis uses assets/phenotypes_resolved.csv, in which conflicting DST
records were resolved by a quality-ranked rule stated in
docs/findings/phenotype_conflicts.md.

Sensitivity analysis repeats the scoring using the original samplesheet labels
(drop_duplicates keep="first"), so the dependence of the conclusions on that
choice is visible rather than assumed away.

Confidence intervals are Wilson score intervals. At n of 13 to 17 the normal
approximation is unreliable near 0 and 1, and an interval that cannot exceed
1.0 matters when several cells are perfect.

Usage:
    python3 bin/score_concordance.py
"""

import csv
import math
import os

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESOLVED = os.path.join(REPO, "assets", "phenotypes_resolved.csv")
SHEET = os.path.join(REPO, "assets", "samplesheet.csv")
PRED = os.path.expanduser(
    "~/mtb-data/results/catalogue/genotypic_predictions.tsv")
OUTDIR = os.path.expanduser("~/mtb-data/results/catalogue")

Z = 1.959963985  # 95%


def wilson(x, n, z=Z):
    """Wilson score interval for a binomial proportion.

    Unlike the normal approximation this never returns bounds outside [0, 1]
    and behaves sensibly when x equals 0 or n, which happens repeatedly here.
    """
    if n == 0:
        return (float("nan"), float("nan"))
    p = x / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = (z / d) * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (max(0.0, centre - half), min(1.0, centre + half))


def load_predictions(path):
    pred = {}
    detail = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("sample"):
                continue
            f = line.rstrip("\n").split("\t")
            pred[(f[0], f[1])] = f[2]
            detail[(f[0], f[1])] = f[4]
    return pred, detail


def load_resolved(path):
    ph = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            ph[(r["run_accession"], r["drug"])] = r["phenotype"]
    return ph


def load_original(path):
    ph = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            ph[(r["run_accession"], "INH")] = r["INH"]
            ph[(r["run_accession"], "RIF")] = r["RIF"]
    return ph


def score(pred, pheno, label):
    print(f"\n{'=' * 60}")
    print(label)
    print("=" * 60)
    results = {}
    for drug in ("INH", "RIF"):
        tp = fp = fn = tn = 0
        for (s, d), p in pred.items():
            if d != drug:
                continue
            o = pheno.get((s, d))
            if o is None:
                continue
            if p == "R" and o == "R":
                tp += 1
            elif p == "R" and o == "S":
                fp += 1
            elif p == "S" and o == "R":
                fn += 1
            else:
                tn += 1

        sens_lo, sens_hi = wilson(tp, tp + fn)
        spec_lo, spec_hi = wilson(tn, tn + fp)

        print(f"\n{drug}")
        print(f"                  phenotype R    phenotype S")
        print(f"  predicted R  {tp:>12}   {fp:>12}")
        print(f"  predicted S  {fn:>12}   {tn:>12}")
        if tp + fn:
            print(f"  sensitivity  {tp/(tp+fn):.3f}  ({tp}/{tp+fn})"
                  f"   95% CI {sens_lo:.3f} to {sens_hi:.3f}")
        if tn + fp:
            print(f"  specificity  {tn/(tn+fp):.3f}  ({tn}/{tn+fp})"
                  f"   95% CI {spec_lo:.3f} to {spec_hi:.3f}")
        results[drug] = (tp, fp, fn, tn)
    return results


def main():
    pred, detail = load_predictions(PRED)
    resolved = load_resolved(RESOLVED)
    original = load_original(SHEET)

    primary = score(pred, resolved, "PRIMARY: quality-ranked phenotypes")
    sens = score(pred, original, "SENSITIVITY: original keep=first labels")

    # ---- per-sample table, primary analysis ----
    out = os.path.join(OUTDIR, "concordance_primary.tsv")
    os.makedirs(OUTDIR, exist_ok=True)
    with open(out, "w") as fh:
        fh.write("sample\tdrug\tpredicted\tphenotype\tagreement\t"
                 "phenotype_changed\tvariants\n")
        for (s, d), p in sorted(pred.items()):
            o = resolved.get((s, d), "")
            orig = original.get((s, d), "")
            if p == o:
                agree = "concordant"
            elif o == "R":
                agree = "false_negative"
            else:
                agree = "false_positive"
            fh.write(f"{s}\t{d}\t{p}\t{o}\t{agree}\t"
                     f"{'yes' if o != orig else 'no'}\t{detail[(s,d)]}\n")

    print(f"\n\nwrote {out}")

    # ---- what changed between the two analyses ----
    print("\nDiscordances under each labelling:")
    for name, ph in (("primary  ", resolved), ("sensitivity", original)):
        disc = [(s, d) for (s, d), p in pred.items()
                if ph.get((s, d)) and p != ph[(s, d)]]
        print(f"  {name}: {len(disc)}  "
              f"{', '.join(f'{s}/{d}' for s, d in sorted(disc))}")

    print("\nNote: sensitivity and specificity are estimable from this "
          "balanced-by-design sample.\nPositive and negative predictive value "
          "are not, because they depend on\nresistance prevalence in the "
          "population tested, which was set by design.")


if __name__ == "__main__":
    main()
