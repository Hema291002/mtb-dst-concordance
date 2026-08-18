import pandas as pd

s = pd.read_csv("assets/samplesheet.csv")

rows = []
for _, r in s.iterrows():
    for n in (1, 2):
        url = r[f"fastq_{n}"]
        rows.append({
            "run": r["run_accession"],
            "filename": url.rsplit("/", 1)[-1],
            "url": url,
            "md5": r[f"md5_{n}"],
        })

m = pd.DataFrame(rows)
m.to_csv("assets/download_manifest.tsv", sep="\t", index=False)

print("files:", len(m))
print("unique runs:", m["run"].nunique())
print("duplicate filenames:", m["filename"].duplicated().sum())
print()
print(m.head(4).to_string(index=False))
