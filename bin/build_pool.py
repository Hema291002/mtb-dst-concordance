import pandas as pd
import re

# ---- phenotypes ----
dst = pd.read_parquet("data/meta/DST_MEASUREMENTS.parquet").reset_index()
dst = dst[dst["DRUG"].isin(["INH", "RIF"])]
dst = dst[dst["PHENOTYPE"].isin(["R", "S"])]
dst = dst[dst["QUALITY"].isin(["HIGH", "MEDIUM"])]

# one row per isolate, columns INH and RIF
dst = dst.drop_duplicates(subset=["UNIQUEID", "DRUG"], keep="first")
pheno = dst.pivot(index="UNIQUEID", columns="DRUG", values="PHENOTYPE")
pheno = pheno.dropna(subset=["INH", "RIF"])
print("isolates with both INH and RIF:", len(pheno))

# ---- sequencing ----
wgs = pd.read_parquet("data/meta/WGS_SAMPLES.parquet").reset_index()
wgs = wgs[wgs["status"] == "complete"]
wgs = wgs[~wgs["has_multiple_ena_run_accessions"]]
wgs = wgs[wgs["fastq_ftp"].notna()]
wgs = wgs[wgs["fastq_ftp"].str.count(";") == 1]     # exactly 2 files = paired
print("usable sequencing runs:", len(wgs))

# ---- genome QC (QC columns only) ----
gen = pd.read_parquet("data/meta/GENOMES.parquet").reset_index()
gen = gen[["UNIQUEID", "SPECIES", "N_LINEAGES", "TB_COVERAGE", "TB_DEPTH"]]
gen = gen[gen["SPECIES"] == "M. tuberculosis"]
gen = gen[gen["N_LINEAGES"] == 1]
gen = gen[gen["TB_COVERAGE"] >= 95]
gen = gen[gen["TB_DEPTH"] >= 40]
print("isolates passing genome QC:", len(gen))

# ---- join ----
pool = pheno.reset_index().merge(wgs, on="UNIQUEID").merge(gen, on="UNIQUEID")
print("pool after join:", len(pool))

# ---- one isolate per subject ----
pool["SUBJECT"] = pool["UNIQUEID"].str.extract(r"subj\.([^.]+)")
pool = pool.sort_values("TB_DEPTH", ascending=False)
pool = pool.drop_duplicates(subset=["SUBJECT"], keep="first")
print("pool after one-per-subject:", len(pool))

# ---- download size ----
pool["MB"] = pool["fastq_bytes"].apply(
    lambda s: sum(int(x) for x in s.split(";")) / 1e6
)

print()
print("INH x RIF contingency:")
print(pd.crosstab(pool["INH"], pool["RIF"]))
print()
print("median MB per isolate:", round(pool["MB"].median(), 1))
print()
print("depth summary:")
print(pool["TB_DEPTH"].describe())

pool.to_parquet("data/meta/candidate_pool.parquet")
print()
print("saved candidate_pool.parquet")
