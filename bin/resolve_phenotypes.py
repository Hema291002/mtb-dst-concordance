#!/usr/bin/env python3
"""Resolve conflicting phenotype records using a quality-ranked rule.

Some isolates have multiple DST records from different sources. Where these
disagree, selection must be principled rather than arbitrary. The original
build_pool.py used drop_duplicates(keep="first"), which selects by file order.

Rule, stated on quality grounds:
  1. prefer PHENOTYPE_QUALITY == HIGH over MEDIUM
  2. among equals, prefer records carrying an MIC value

This rule is corroborated by the independent UKMYC plate assay, which agrees
with every record the rule selects (see docs/findings/phenotype_conflicts.md).

Sample selection is NOT repeated. The same 30 isolates are retained; only the
phenotype label is corrected.
"""
import csv
import pandas as pd

SHEET = "assets/samplesheet.csv"
DST = "data/meta/DST_MEASUREMENTS.parquet"
OUT = "assets/phenotypes_resolved.csv"

rows = list(csv.DictReader(open(SHEET)))
byuid = {r["UNIQUEID"]: r for r in rows}

d = pd.read_parquet(DST).reset_index()
d = d[d["UNIQUEID"].isin(byuid) & d["DRUG"].isin(["INH", "RIF"])]
d = d[d["PHENOTYPE"].isin(["R", "S"])]
d["qr"] = d["QUALITY"].map({"HIGH": 0, "MEDIUM": 1, "LOW": 2})
d["hm"] = d["METHOD_MIC"].notna().map({True: 0, False: 1})
best = (d.sort_values(["UNIQUEID", "DRUG", "qr", "hm"])
          .drop_duplicates(subset=["UNIQUEID", "DRUG"], keep="first"))

resolved = {(r["UNIQUEID"], r["DRUG"]): (r["PHENOTYPE"], r["QUALITY"],
                                         r["SOURCE"])
            for _, r in best.iterrows()}

changed = 0
with open(OUT, "w", newline="") as fh:
    w = csv.writer(fh, lineterminator="\n")
    w.writerow(["run_accession", "UNIQUEID", "drug", "phenotype",
                "quality", "source", "changed_from_original"])
    for uid, r in sorted(byuid.items()):
        for drug in ("INH", "RIF"):
            ph, q, src = resolved[(uid, drug)]
            orig = r[drug]
            flag = "yes" if ph != orig else "no"
            if flag == "yes":
                changed += 1
            w.writerow([r["run_accession"], uid, drug, ph, q, src, flag])

print(f"wrote {OUT}")
print(f"{changed} label(s) differ from the original samplesheet")
