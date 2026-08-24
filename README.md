# Genotypic prediction of isoniazid and rifampicin resistance in *Mycobacterium tuberculosis*: concordance with laboratory phenotype in 30 clinical isolates

Variant calls generated independently from raw sequencing reads are compared against WHO-catalogued resistance markers and matched laboratory drug-susceptibility results, to quantify how well genotype predicts phenotype and to characterise every case where it does not.

---

## Biological question

**Primary.** In 30 clinical *M. tuberculosis* isolates with matched phenotypic drug-susceptibility testing, what are the sensitivity and specificity of resistance prediction based on the WHO mutation catalogue, using variant calls produced by this pipeline rather than by the consortium that generated the data?

**Secondary.** For every discordant isolate, which explanations can be ruled out with the available evidence, and which cannot?

**Tertiary, descriptive only.** Does lineage, assigned independently from the same sequence data, co-occur with resistance or with discordance? Underpowered at this sample size and reported as description, not as a test.

## Background

Tuberculosis treatment depends on knowing which drugs still work. The reference method is to culture the organism and expose it to each drug, which is reliable but slow because *M. tuberculosis* grows slowly. Sequencing the genome and looking for known resistance mutations is faster, but only useful if the prediction can be trusted.

The two drugs studied here behave differently. **Rifampicin** resistance arises almost entirely from a short, well-characterised region of *rpoB*, and genotypic prediction performs near-perfectly in published work. **Isoniazid** is a prodrug requiring activation by the *katG* catalase-peroxidase; resistance arises through *katG* loss of function, through upregulation of the drug's target via the *fabG1*-*inhA* operon promoter, and through mechanisms that remain incompletely catalogued. Concordance is consistently lower.

This project reproduces that contrast from raw reads, and treats the disagreements as the object of study rather than as residual error.

## Dataset

**Source:** CRyPTIC consortium, release v3.4.0 (21 May 2025), Zenodo record [15680920](https://zenodo.org/records/15680920), CC-BY 4.0.

CRyPTIC is unusual in providing, for the same isolates, both whole-genome sequence and matched laboratory drug-susceptibility results produced by independent methods. That independence is what makes this an experiment rather than a demonstration.

| | |
|---|---|
| Organism | *Mycobacterium tuberculosis* |
| Isolates | 30, selected from 37,962 passing all criteria |
| Sequencing | Illumina paired-end, 80-151 bp reads |
| Raw data | 60 FASTQ files, 9.24 GB, all MD5-verified |
| Depth after alignment | 52.2x to 118.6x (median ~90x) |
| Reference | H37Rv, RefSeq GCF_000195955.2 (ASM19595v2), NC_000962.3, 4,411,532 bp |
| Geography | 10 collecting sites plus 6 direct ENA submissions, spanning Peru, South Africa, India, Vietnam, Germany, UK, Canada, Italy, China |
| Accessions | listed in [`assets/samplesheet.csv`](assets/samplesheet.csv) |

Isolates were selected by committed script ([`bin/build_pool.py`](bin/build_pool.py), [`bin/select_samples.py`](bin/select_samples.py)) with a fixed random seed, balanced 8/7/7/8 across the four isoniazid-by-rifampicin resistance categories, giving 15 resistant and 15 susceptible per drug.

**Design consequence.** Because the resistant proportion was set by design, sensitivity and specificity are estimable but positive and negative predictive value are not: those depend on resistance prevalence in the population tested. PPV and NPV are therefore not reported.

## Workflow

```
CRyPTIC metadata (Zenodo)
        |
        v
Sample selection  ---- explicit criteria, fixed seed, committed script
        |
        v
FASTQ download  ---- MD5-verified (2 of 60 files arrived corrupted)
        |
        v
Raw read QC  ---- FastQC + MultiQC
        |
        v
Trimming  ---- fastp, adapter auto-detection
        |
        v
Alignment to H37Rv  ---- bwa-mem2, per-sample read groups
        |
        v
Duplicate marking  ---- samtools markdup (marked, not removed)
        |
        v
Coverage assessment  ---- mosdepth, per resistance gene
        |
        v
Variant calling  ---- bcftools, --ploidy 1, no BQSR
        |
        v
Filtering  ---- QUAL, depth, mapping quality, allele fraction, repeat mask
        |
        v
Annotation  ---- custom, validated against published positions
        |
        v
WHO catalogue lookup  ---- coordinate matching, not name matching
        |
        v
Phenotype comparison  ---- primary + sensitivity analysis
        |
        v
Lineage assignment  ---- SNP barcode, same VCFs
```

## Tools

| Tool | Version | Purpose |
|---|---|---|
| FastQC | 0.12.1 | Per-file read quality assessment |
| MultiQC | 1.35 | Aggregation across the cohort |
| fastp | 1.3.6 | Adapter removal and quality trimming |
| bwa-mem2 | 2.2.1 | Read alignment to H37Rv |
| samtools | 1.24 | BAM handling, duplicate marking, indexing |
| mosdepth | 0.3.14 | Per-region coverage |
| bcftools | 1.24 | Variant calling and filtering |
| pandas / pyarrow | | Metadata handling (Parquet) |

Environments are declared in [`envs/`](envs/) and each stage runs in its own.

## Installation

```bash
git clone https://github.com/Hema291002/mtb-dst-concordance.git
cd mtb-dst-concordance

conda env create -f envs/mtb-meta.yml     # metadata, selection
conda env create -f envs/mtb-qc.yml       # QC and trimming
conda env create -f envs/mtb-align.yml    # alignment, coverage
conda env create -f envs/mtb-call.yml     # variant calling
```

Reference genome:

```bash
mkdir -p ~/mtb-data/ref && cd ~/mtb-data/ref
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.fna.gz
curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/195/955/GCF_000195955.2_ASM19595v2/GCF_000195955.2_ASM19595v2_genomic.gff.gz
gunzip -k GCF_000195955.2_ASM19595v2_genomic.fna.gz
bwa-mem2 index GCF_000195955.2_ASM19595v2_genomic.fna
samtools faidx GCF_000195955.2_ASM19595v2_genomic.fna
```

## Usage

```bash
conda activate mtb-meta
python3 bin/build_pool.py                 # apply inclusion criteria
python3 bin/select_samples.py             # balanced selection, seed 42
python3 bin/make_manifest.py              # download manifest with checksums
bash    bin/download_fastq.sh             # resumable, MD5-verified

conda activate mtb-qc
bash    bin/trim_reads.sh

conda activate mtb-align
bash    bin/align_reads.sh
bash    bin/markdup.sh

conda activate mtb-call
bash    bin/call_variants.sh
bash    bin/filter_variants.sh
python3 bin/annotate_targets.py --selftest    # validate before use
python3 bin/catalogue_lookup.py --variants docs/target_variants.tsv ...
python3 bin/assign_lineage.py --barcode ... --vcfdir ... --out ...

conda activate mtb-meta
python3 bin/resolve_phenotypes.py
python3 bin/score_concordance.py
```

Every script is resumable: completed samples are detected and skipped.

---

## Results

### Genotype-phenotype concordance

Primary analysis. 95% Wilson score intervals.

| Drug | Sensitivity | 95% CI | Specificity | 95% CI |
|---|---|---|---|---|
| Isoniazid | 0.824 (14/17) | 0.590 – 0.938 | 1.000 (13/13) | 0.772 – 1.000 |
| Rifampicin | 1.000 (13/13) | 0.772 – 1.000 | 0.882 (15/17) | 0.657 – 0.967 |

**The intervals are the result, not the point estimates.** Every lower bound falls between 0.59 and 0.77. This study cannot distinguish good performance from perfect performance, and cannot independently establish that rifampicin prediction outperforms isoniazid prediction, although the literature supports that relationship.

A sensitivity analysis using the alternative phenotype labelling gives intervals that overlap the primary analysis throughout, so the conclusions do not depend on that choice. See [`docs/findings/phenotype_conflicts.md`](docs/findings/phenotype_conflicts.md).

### Variants found

90 variants across *katG*, *rpoB*, *inhA*, *fabG1* and *ahpC*: 48 missense, 25 synonymous, 17 promoter. 88 of 90 matched the WHO catalogue.

Signature resistance mutations recovered: `katG S315T` (11 isolates), `rpoB S450L` (13), `rpoB H445Y` (2), `fabG1 c.-15C>T` (3).

### Lineage

| Lineage | n |
|---|---|
| 4 | 16 |
| 2 | 6 |
| 3 | 5 |
| 1 | 2 |
| unassigned | 1 |

No isolate carried off-path barcode markers, so no evidence of mixed infection. See [`docs/findings/lineage_assignment.md`](docs/findings/lineage_assignment.md) for why one isolate is unassigned and why that absence is informative.

---

## Biological interpretation

**Presence in a resistance gene does not imply a resistance mechanism.** `katG R463L` appeared in 13 isolates, more than carried `katG S315T`. The WHO catalogue grades it **5, not associated with resistance**. It is present in exactly the 13 non-lineage-4 isolates and absent from all 16 lineage-4 isolates: a phylogenetic marker, not a resistance mutation. Scoring every non-synonymous *katG* change as resistance would have produced 13 false positives.

**Promoter regions must be included.** The second commonest isoniazid resistance mutation lies 15 bases upstream of *fabG1*, outside the gene body. Target regions here extend 200 bp upstream, strand-aware. Restricting to coding sequence would have missed it in 3 isolates.

**Two discordances are robust and unexplained.** ERR4814489 and ERR8975559 are phenotypically isoniazid-resistant with no graded tier 1 variant, under both phenotype labellings. Both have complete 30x coverage across all isoniazid genes, ruling out a coverage gap or a gene deletion. Neither carries any sub-threshold variant, ruling out heteroresistance or a filtering artefact. ERR8975559 has an MIC greater than 12.8 against a critical concentration of 0.1, making phenotype error implausible. **A tier 2 gene or an uncatalogued mechanism remains.**

**One finding is invisible to variant calling entirely.** ERR8975061 carries a 3,275 bp homozygous deletion removing *rraA*, *ethA* and *ethR*. Loss of *ethA*, which activates ethionamide, is a documented resistance mechanism. No VCF can express absent sequence; this was found only by coverage analysis. See [`docs/findings/ethA_deletion_ERR8975061.md`](docs/findings/ethA_deletion_ERR8975061.md).

**Ground truth is not singular.** For 4 of 60 isolate-drug pairs, two accredited laboratories using different methods reached opposite conclusions about the same isolate. This bounds what any genotypic prediction can achieve: it cannot exceed the reproducibility of the standard it is scored against.

---

## Reproducibility

- Every analytical decision exists as a committed script, not as a manual step.
- Sample selection uses a fixed seed, so the same 30 isolates are recoverable.
- All 60 FASTQ files are MD5-verified against published checksums ([`assets/fastq_checksums.md5`](assets/fastq_checksums.md5)). Two arrived corrupted on first download; one failed silently and was caught only by a systematic verification pass.
- Software environments are declared per stage in [`envs/`](envs/).
- The custom annotator refuses to run until it reproduces three independently published variant positions ([`bin/annotate_targets.py --selftest`](bin/annotate_targets.py)).
- [`docs/ANALYSIS_STATUS.md`](docs/ANALYSIS_STATUS.md) records what has been run versus what is planned, and every limitation was recorded when identified rather than assembled afterwards.

### Pre-registration

Three predictions were committed to git **before** the relevant data was consulted, and all three were confirmed:

| Prediction | Committed in | Outcome |
|---|---|---|
| `katG R463L` grades as not associated with resistance despite appearing in 13 isolates | [`catalogue_predictions.md`](docs/findings/catalogue_predictions.md) | Confirmed, grade 5 |
| Variant-count clusters correspond to lineages ordered by distance from H37Rv | [`variant_count_clusters.md`](docs/findings/variant_count_clusters.md) | Confirmed; one sub-claim failed and is recorded |
| ERR8975061 is ethionamide-resistant given its *ethA* deletion | [`ethA_deletion_ERR8975061.md`](docs/findings/ethA_deletion_ERR8975061.md) | Pending phenotype check |

The consortium's own resistance predictions and lineage assignments were quarantined before analysis began ([`docs/held_out_rationale.md`](docs/held_out_rationale.md)) so that filtering thresholds could not be tuned toward a known answer.

---

## Limitations

- **Sample size.** 95% confidence intervals span roughly 0.59 to 1.00. The study cannot distinguish good from perfect performance for any measure.
- **Selection bias.** CRyPTIC deliberately enriched for resistant isolates. These estimates do not generalise to a routine diagnostic case mix.
- **Confounding.** Resistance category is partly confounded with collecting site, and site predicts lineage. Some sites contributed exclusively resistant or exclusively susceptible isolates.
- **Scope.** Only WHO tier 1 genes were examined. A tier 2 variant could explain a phenotypically resistant isolate with no graded tier 1 mutation.
- **Structural variation.** Gene-scale deletions are invisible to variant calling. One was found by coverage analysis; genes not examined for coverage could carry undetected deletions.
- **Reference bias.** H37Rv is lineage 4. Three isolates mapped 1-3% lower than the rest, consistent with divergence from the reference. Accessory sequence absent from H37Rv cannot be assessed at all.
- **Masking.** The repeat mask removes both mapping artefacts and genuine variation. PE/PPE gene families are hypervariable, so some masked variants are real. No claims are made about PE/PPE loci.
- **Technical heterogeneity.** Library preparation differs by era (80/100/150 bp reads with correspondingly different insert sizes), so mapping resolution is not uniform across the cohort.
- **Indels.** Filters were applied identically to SNPs and indels. Short-read indel calling is less reliable, and some *katG* resistance mutations are frameshifts. No indels were found in the target genes here.
- **Lineage method.** Assignment from a variants-only VCF cannot detect lineage 4.9, because H37Rv belongs to it and carries its defining alleles.

---

## Status

**Completed:** sample selection, data acquisition, QC, trimming, alignment, duplicate marking, coverage assessment, variant calling, filtering, annotation, WHO catalogue interpretation, phenotype comparison, lineage assignment.

**Not yet done:** R analysis and publication-quality figures; Nextflow implementation of the pipeline.

The bash and Python scripts here are the working pipeline. Nextflow would add resumability across machines, containerisation and portability, and is planned rather than claimed. See [`docs/ANALYSIS_STATUS.md`](docs/ANALYSIS_STATUS.md) for the per-stage record.

## Future improvements

- Nextflow DSL2 implementation with containerised processes and a test profile
- R analysis: variant frequency by lineage, coverage distributions, transition/transversion ratios, confidence interval visualisation
- Extension to WHO tier 2 genes, to test whether either robust discordance is explained
- Independent variant calling with freebayes as a concordance check
- Computed mappability track in place of the gene-family-based mask

---

## Citations

**Dataset.** The CRyPTIC Consortium. CRyPTIC Consortium Dataset v3.4.0. Zenodo, 2025. https://doi.org/10.5281/zenodo.15680920

**Reference genome.** Cole ST et al. Deciphering the biology of *Mycobacterium tuberculosis* from the complete genome sequence. *Nature* 393:537-544 (1998). Assembly GCF_000195955.2.

**Resistance catalogue.** World Health Organization. *Catalogue of mutations in Mycobacterium tuberculosis complex and their association with drug resistance*, 2nd edition. Geneva, 2023. ISBN 978-92-4-008241-0.

**Lineage barcode.** Coll F et al. A robust SNP barcode for typing *Mycobacterium tuberculosis* complex strains. *Nat Commun* 5:4812 (2014). Updated in Napier G et al. Robust barcoding and identification of *Mycobacterium tuberculosis* lineages for epidemiological and clinical studies. *Genome Med* 12:114 (2020). Barcode file from jodyphelan/tbdb.

**Tools.** Li H. Aligning sequence reads with BWA-MEM (2013), and Vasimuddin M et al. Efficient architecture-aware acceleration of BWA-MEM (IPDPS 2019). Danecek P et al. Twelve years of SAMtools and BCFtools. *GigaScience* 10:giab008 (2021). Chen S et al. fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics* 34:i884 (2018). Pedersen BS, Quinlan AR. Mosdepth. *Bioinformatics* 34:867 (2018). Ewels P et al. MultiQC. *Bioinformatics* 32:3047 (2016).

## Licence

Code in this repository is released under the MIT Licence. CRyPTIC data is CC-BY 4.0; the WHO catalogue is subject to WHO terms of use.
