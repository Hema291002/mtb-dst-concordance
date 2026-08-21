import pandas as pd

pd.set_option("display.width", 200)
pd.set_option("display.max_columns", 60)

for t in ["WGS_SAMPLES", "DST_SAMPLES", "DST_MEASUREMENTS", "GENOMES"]:
    df = pd.read_parquet(f"data/meta/{t}.parquet")
    print("=" * 70)
    print(t, df.shape)
    print(list(df.columns))
    print(df.head(3).to_string())
    print()

dst = pd.read_parquet("data/meta/DST_MEASUREMENTS.parquet").reset_index()

print("=" * 70)
print("Columns after reset_index:", list(dst.columns))
print()

sub = dst[dst["DRUG"].isin(["INH", "RIF"])]
print("INH/RIF rows:", len(sub))
print()
print(pd.crosstab(sub["DRUG"], sub["PHENOTYPE"], dropna=False))
print()
print(pd.crosstab(sub["DRUG"], sub["QUALITY"], dropna=False))
print()
print("SOURCE values for INH/RIF:")
print(sub["SOURCE"].value_counts().head(10))
