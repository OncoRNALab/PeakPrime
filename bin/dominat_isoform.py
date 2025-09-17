#!/usr/bin/env python3

import pandas as pd
import argparse
import glob
import os

def parse_gtf_attributes(attr_str):
    """Parse GTF attribute string into dictionary"""
    attrs = {}
    for entry in attr_str.strip().split(";"):
        if entry.strip() == "":
            continue
        key, value = entry.strip().split(" ", 1)
        attrs[key] = value.strip('"')
    return attrs

def load_quant_file(path):
    """Load Salmon quant.sf or Kallisto abundance.tsv"""
    df = pd.read_csv(path, sep="\t")
    if "Name" in df.columns and "TPM" in df.columns:
        df = df[["Name", "TPM"]].rename(columns={"Name": "transcript_id", "TPM": os.path.basename(path)})
    elif "target_id" in df.columns and "tpm" in df.columns:
        df = df[["target_id", "tpm"]].rename(columns={"target_id": "transcript_id", "tpm": os.path.basename(path)})
    else:
        raise ValueError(f"File {path} not recognized as Salmon or Kallisto format")
    return df.set_index("transcript_id")

def main():
    parser = argparse.ArgumentParser(
        description="Generate a dominant isoform GTF using multiple Salmon or Kallisto quantifications."
    )
    parser.add_argument("-g", "--gtf", required=True, help="Input GTF file")
    parser.add_argument("-q", "--quant", nargs="+", required=True,
                        help="One or more quant.sf (Salmon) or abundance.tsv (Kallisto) files")
    parser.add_argument("-o", "--out", required=True, help="Output GTF file")

    args = parser.parse_args()

    # Load all quantifications and merge
    dfs = [load_quant_file(q) for q in args.quant]
    quant_all = pd.concat(dfs, axis=1).fillna(0)

    # Compute mean TPM across samples
    quant_all["mean_tpm"] = quant_all.mean(axis=1)

    # Map transcript -> gene from GTF
    transcript_to_gene = {}
    with open(args.gtf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[2] != "transcript":
                continue
            attrs = parse_gtf_attributes(fields[8])
            if "transcript_id" in attrs and "gene_id" in attrs:
                transcript_to_gene[attrs["transcript_id"]] = attrs["gene_id"]

    # Pick dominant isoform per gene (highest mean TPM)
    gene_to_best_tx = {}
    for tx, gene in transcript_to_gene.items():
        if tx not in quant_all.index:
            continue
        tpm = quant_all.loc[tx, "mean_tpm"]
        if gene not in gene_to_best_tx or tpm > gene_to_best_tx[gene][1]:
            gene_to_best_tx[gene] = (tx, tpm)

    dominant_tx = set(tx for tx, _ in gene_to_best_tx.values())

    # Filter GTF to only dominant isoforms
    with open(args.out, "w") as out_f:
        with open(args.gtf) as f:
            for line in f:
                if line.startswith("#"):
                    out_f.write(line)
                    continue
                fields = line.strip().split("\t")
                attrs = parse_gtf_attributes(fields[8])
                if "transcript_id" in attrs:
                    if attrs["transcript_id"] in dominant_tx:
                        out_f.write(line)

    print(f"✅ Wrote dominant isoform GTF with {len(dominant_tx)} transcripts to {args.out}")

if __name__ == "__main__":
    main()
