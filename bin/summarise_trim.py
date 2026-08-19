import glob, json, os

rows = []
for f in sorted(glob.glob(os.path.expanduser("~/mtb-data/results/qc/fastp/*.fastp.json"))):
    d = json.load(open(f))
    run = os.path.basename(f).replace(".fastp.json", "")
    before = d["summary"]["before_filtering"]
    after = d["summary"]["after_filtering"]
    filt = d["filtering_result"]
    rows.append({
        "run": run,
        "pct_passed": 100 * after["total_reads"] / before["total_reads"],
        "pct_adapter": 100 * d["adapter_cutting"]["adapter_trimmed_reads"] / before["total_reads"]
                       if "adapter_cutting" in d else 0.0,
        "insert": d.get("insert_size", {}).get("peak", None),
        "q30_after": 100 * after["q30_rate"],
        "meanlen_after": after["read1_mean_length"],
        "cov_after": after["total_bases"] / 4411532,
    })

rows.sort(key=lambda r: r["cov_after"])
print(f"{'run':<14}{'pass%':>7}{'adapt%':>8}{'insert':>8}{'Q30%':>7}{'len':>6}{'cov':>7}")
for r in rows:
    print(f"{r['run']:<14}{r['pct_passed']:>7.1f}{r['pct_adapter']:>8.1f}"
          f"{str(r['insert']):>8}{r['q30_after']:>7.1f}{r['meanlen_after']:>6}{r['cov_after']:>7.0f}")

worst = min(rows, key=lambda r: r["pct_passed"])
print(f"\nlowest pass rate: {worst['run']} at {worst['pct_passed']:.1f}%")
print(f"lowest coverage:  {rows[0]['run']} at {rows[0]['cov_after']:.0f}x")
