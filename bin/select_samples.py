import pandas as pd

SEED = 42
DEPTH_MIN, DEPTH_MAX = 50, 120
TARGETS = {("R", "R"): 8, ("R", "S"): 7, ("S", "R"): 7, ("S", "S"): 8}

pool = pd.read_parquet("data/meta/candidate_pool.parquet")
pool = pool[pool["TB_DEPTH"].between(DEPTH_MIN, DEPTH_MAX)]
print("pool within depth band:", len(pool))
print(pd.crosstab(pool["INH"], pool["RIF"]))
print()

picks = []
for (inh, rif), n in TARGETS.items():
    cell = pool[(pool["INH"] == inh) & (pool["RIF"] == rif)]
    if len(cell) < n:
        raise SystemExit(f"cell {inh}/{rif} has only {len(cell)} isolates")
    picks.append(cell.sample(n=n, random_state=SEED))

sel = pd.concat(picks).sort_values("UNIQUEID").reset_index(drop=True)

# split the semicolon-separated FASTQ fields into columns
ftp = sel["fastq_ftp"].str.split(";", expand=True)
md5 = sel["fastq_md5"].str.split(";", expand=True)
sel["fastq_1"] = "https://" + ftp[0]
sel["fastq_2"] = "https://" + ftp[1]
sel["md5_1"] = md5[0]
sel["md5_2"] = md5[1]

cols = ["UNIQUEID", "run_accession", "sample_accession", "INH", "RIF",
        "TB_DEPTH", "TB_COVERAGE", "dataset", "center_name",
        "fastq_1", "fastq_2", "md5_1", "md5_2", "MB"]
sel[cols].to_csv("assets/samplesheet.csv", index=False)

print("selected:", len(sel))
print(pd.crosstab(sel["INH"], sel["RIF"]))
print()
print("total download GB:", round(sel["MB"].sum() / 1000, 2))
print()
print("datasets represented:")
print(sel["dataset"].value_counts())
print()
print("wrote assets/samplesheet.csv")
