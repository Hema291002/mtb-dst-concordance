#!/usr/bin/env python3
"""
Annotate variants falling in resistance-associated genes with codon number,
reference and alternate amino acids, or promoter position.

Written because bcftools csq silently indexes zero CDS features from an NCBI
RefSeq bacterial GFF: RefSeq links CDS directly to gene (Parent=gene-RvXXXX)
with no intervening transcript, and csq requires gene -> transcript -> CDS.

Reads VCF records on stdin in the format produced by:
    bcftools query -f '%POS\t%REF\t%ALT\t%QUAL\t[%DP]\t[%AD]\n'

Writes a tab-separated table on stdout.

Usage:
    bcftools query -f '%POS\t%REF\t%ALT\t%QUAL\t[%DP]\t[%AD]\n' sample.vcf.gz \
      | python3 bin/annotate_targets.py --sample SAMPLENAME

    python3 bin/annotate_targets.py --selftest
"""

import argparse
import gzip
import os
import re
import sys

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

REF_PATH = os.path.expanduser(
    "~/mtb-data/ref/GCF_000195955.2_ASM19595v2_genomic.fna")
GFF_PATH = os.path.expanduser(
    "~/mtb-data/ref/GCF_000195955.2_ASM19595v2_genomic.gff.gz")

# Genes relevant to isoniazid and rifampicin. Others can be added, but every
# gene added is a gene whose annotation must be trusted.
GENES = ["katG", "rpoB", "inhA", "fabG1", "ahpC"]

# How far upstream of a gene to treat a variant as promoter-region.
PROMOTER = 200

# ---------------------------------------------------------------------------
# Genetic code
# ---------------------------------------------------------------------------
#
# NCBI translation table 11 (bacterial). For amino acid assignment this is
# identical to the standard table; table 11 differs only in which codons are
# permitted as translation starts, which does not affect substitutions.
#
# The table is built rather than typed out, which removes the possibility of a
# typo in 64 hand-entered entries. Codons are generated in NCBI's canonical
# order (T, C, A, G at each position, first base varying slowest) so they line
# up with the standard amino acid string.

_BASES = "TCAG"
_AAS = ("FFLLSSSSYY**CC*W"
        "LLLLPPPPHHQQRRRR"
        "IIIMTTTTNNKKSSRR"
        "VVVVAAAADDEEGGGG")

_CODONS = [x + y + z for x in _BASES for y in _BASES for z in _BASES]
CODON_TABLE = dict(zip(_CODONS, _AAS))

# Three-letter to one-letter is not needed; WHO catalogue notation uses
# one-letter codes for substitutions (e.g. katG S315T).

_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def revcomp(seq):
    """Reverse complement of a DNA string."""
    return seq.translate(_COMPLEMENT)[::-1]


def translate(codon):
    """Translate a three-base codon. Returns '?' if the codon is not clean."""
    codon = codon.upper()
    if len(codon) != 3 or any(b not in "ACGT" for b in codon):
        return "?"
    return CODON_TABLE[codon]


# ---------------------------------------------------------------------------
# Loading reference data
# ---------------------------------------------------------------------------

def load_reference(path):
    """Load a single-contig FASTA into one uppercase string.

    The genome is 4.4 Mb, so holding it in memory costs a few megabytes and
    makes every lookup a simple string slice.
    """
    parts = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            parts.append(line.strip())
    return "".join(parts).upper()


def load_genes(path, wanted):
    """Extract (start, end, strand) for the named genes from a GFF3 file.

    Coordinates stay 1-based inclusive, as GFF defines them. Conversion to
    0-based happens only at the point of slicing the sequence, in one place,
    which keeps the off-by-one risk contained.
    """
    genes = {}
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[2] != "gene":
                continue
            m = re.search(r"gene=([^;]+)", f[8])
            if not m:
                continue
            name = m.group(1)
            if name in wanted:
                genes[name] = (int(f[3]), int(f[4]), f[6])
    missing = set(wanted) - set(genes)
    if missing:
        sys.exit(f"ERROR: genes not found in GFF: {sorted(missing)}")
    return genes


# ---------------------------------------------------------------------------
# Annotation
# ---------------------------------------------------------------------------

def locate(pos, genes):
    """Which gene, if any, does this position belong to?

    Returns (gene_name, region) where region is 'cds' or 'promoter',
    or (None, None) if the position is outside every target.
    """
    for name, (start, end, strand) in genes.items():
        if start <= pos <= end:
            return name, "cds"
        if strand == "+" and start - PROMOTER <= pos < start:
            return name, "promoter"
        if strand == "-" and end < pos <= end + PROMOTER:
            return name, "promoter"
    return None, None


def annotate_snv(pos, ref, alt, gene, region, genes, genome):
    """Annotate a single-nucleotide variant.

    Returns a dict of annotation fields.
    """
    start, end, strand = genes[gene]

    # ---- sanity: does the reference base in the VCF match the genome? ----
    genome_base = genome[pos - 1]          # 1-based -> 0-based
    if genome_base != ref.upper():
        return {"effect": "REF_MISMATCH",
                "detail": f"VCF says {ref}, genome has {genome_base}"}

    # ---- promoter variants ----
    if region == "promoter":
        if strand == "+":
            offset = start - pos           # 15 means 15 bases before the gene
            r, a = ref, alt
        else:
            offset = pos - end
            r, a = revcomp(ref), revcomp(alt)
        return {"effect": "promoter",
                "codon": "",
                "aa_ref": "",
                "aa_alt": "",
                "notation": f"{gene} c.-{offset}{r}>{a}"}

    # ---- coding variants ----
    if strand == "+":
        offset = pos - start
        codon_index = offset // 3          # 0-based codon number
        in_codon = offset % 3              # 0, 1 or 2
        codon_start = start + codon_index * 3
        ref_codon = genome[codon_start - 1: codon_start + 2]
        alt_base = alt.upper()
    else:
        offset = end - pos
        codon_index = offset // 3
        in_codon = offset % 3
        # Codon occupies forward coordinates [codon_hi-2 .. codon_hi]
        codon_hi = end - codon_index * 3
        fwd_codon = genome[codon_hi - 3: codon_hi]
        ref_codon = revcomp(fwd_codon)
        alt_base = revcomp(alt.upper())

    if len(ref_codon) != 3:
        return {"effect": "TRUNCATED_CODON", "detail": ""}

    alt_codon = ref_codon[:in_codon] + alt_base + ref_codon[in_codon + 1:]

    aa_ref = translate(ref_codon)
    aa_alt = translate(alt_codon)
    codon_number = codon_index + 1

    if aa_ref == "?" or aa_alt == "?":
        effect = "unknown"
    elif aa_ref == aa_alt:
        effect = "synonymous"
    elif aa_alt == "*":
        effect = "stop_gained"
    elif aa_ref == "*":
        effect = "stop_lost"
    else:
        effect = "missense"

    return {"effect": effect,
            "codon": f"{ref_codon}>{alt_codon}",
            "aa_ref": aa_ref,
            "aa_alt": aa_alt,
            "notation": f"{gene} {aa_ref}{codon_number}{aa_alt}",
            "codon_number": codon_number}


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

def selftest(genes, genome):
    """Check the annotator against positions whose answers are published.

    These are run before the annotator is trusted on real data. A failure here
    means the coordinate arithmetic, strand handling or codon table is wrong,
    and any downstream result would be silently incorrect.
    """
    cases = [
        # (position, ref, alt, expected notation, description)
        (761155, "C", "T", "rpoB S450L",
         "forward strand, middle base of codon"),
        (2155168, "C", "G", "katG S315T",
         "reverse strand, middle base of codon"),
        (1673425, "C", "T", "fabG1 c.-15C>T",
         "forward strand promoter, 15 bases upstream"),
    ]

    failures = 0
    print("Self-test")
    print("-" * 72)
    for pos, ref, alt, expected, desc in cases:
        gene, region = locate(pos, genes)
        if gene is None:
            print(f"FAIL  {pos}  not inside any target gene")
            failures += 1
            continue
        res = annotate_snv(pos, ref, alt, gene, region, genes, genome)
        got = res.get("notation", res.get("effect", "?"))
        status = "PASS" if got == expected else "FAIL"
        if status == "FAIL":
            failures += 1
        print(f"{status}  {pos:>9}  {ref}>{alt}")
        print(f"        expected: {expected}")
        print(f"        got:      {got}")
        print(f"        ({desc})")
        print()

    # Extra check: the reference codon at katG 315 must be serine.
    start, end, strand = genes["katG"]
    offset = end - 2155168
    codon_hi = end - (offset // 3) * 3
    codon = revcomp(genome[codon_hi - 3: codon_hi])
    aa = translate(codon)
    status = "PASS" if aa == "S" else "FAIL"
    if status == "FAIL":
        failures += 1
    print(f"{status}  katG codon 315 reference sequence")
    print(f"        forward strand: {genome[codon_hi - 3: codon_hi]}")
    print(f"        reverse comp:   {codon}  ->  {aa} (expect S)")
    print()

    print("-" * 72)
    if failures:
        print(f"{failures} CHECK(S) FAILED. Do not use this annotator.")
        return 1
    print("All checks passed.")
    return 0


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", default="NA",
                    help="sample name to include in the output table")
    ap.add_argument("--selftest", action="store_true",
                    help="run validation checks and exit")
    ap.add_argument("--ref", default=REF_PATH)
    ap.add_argument("--gff", default=GFF_PATH)
    ap.add_argument("--header", action="store_true",
                    help="print a header line")
    args = ap.parse_args()

    genome = load_reference(args.ref)
    genes = load_genes(args.gff, GENES)

    if args.selftest:
        sys.exit(selftest(genes, genome))

    if args.header:
        print("\t".join(["sample", "pos", "ref", "alt", "gene", "region",
                         "effect", "codon", "aa_ref", "aa_alt", "notation",
                         "qual", "dp", "ad"]))

    for line in sys.stdin:
        line = line.rstrip("\n")
        if not line:
            continue
        f = line.split("\t")
        if len(f) < 3:
            continue
        pos = int(f[0])
        ref = f[1]
        alt = f[2]
        qual = f[3] if len(f) > 3 else ""
        dp = f[4] if len(f) > 4 else ""
        ad = f[5] if len(f) > 5 else ""

        gene, region = locate(pos, genes)
        if gene is None:
            continue                      # outside every target gene

        # Indels are reported but not translated. Short-read indel calling is
        # less reliable than SNV calling, and frameshift consequences depend on
        # downstream sequence in a way that deserves separate treatment.
        if len(ref) != 1 or len(alt) != 1:
            kind = "insertion" if len(alt) > len(ref) else "deletion"
            shift = "frameshift" if abs(len(alt) - len(ref)) % 3 else "in_frame"
            print("\t".join([args.sample, str(pos), ref, alt, gene, region,
                             f"{kind}_{shift}", "", "", "",
                             f"{gene} {ref}>{alt} @ {pos}", qual, dp, ad]))
            continue

        res = annotate_snv(pos, ref, alt, gene, region, genes, genome)
        print("\t".join([args.sample, str(pos), ref, alt, gene, region,
                         res.get("effect", ""), res.get("codon", ""),
                         res.get("aa_ref", ""), res.get("aa_alt", ""),
                         res.get("notation", res.get("detail", "")),
                         qual, dp, ad]))


if __name__ == "__main__":
    main()
