#!/usr/bin/env python3
"""
Select best primers from alignment reports.

Specificity policy (3′-end protocols):
  1. Infer intended gene from 0-mismatch transcriptome hits.
  2. Consider *relevant* alignments: mismatches ≤ max_mismatches AND
     distance_to_end ≤ distance_threshold.
  3. Keep primer only if the intended gene is present in relevant hits and
     no OTHER gene appears in relevant hits.

Off-targets far from the transcript 3′ end are ignored (not harmful for
3′-end sequencing). Off-targets near the 3′ end with ≤ max_mismatches reject
the primer.

Writes 'best_primers.tsv' (or --out).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

# Allow importing sibling helper when run as bin/select_best_primer.py
sys.path.insert(0, str(Path(__file__).resolve().parent))
from primer_specificity import (  # noqa: E402
    classify_primer,
    group_alignments_by_primer,
)


def _ensure_primer_id(summary: pd.DataFrame) -> pd.DataFrame:
    if "primer_id" in summary.columns:
        return summary
    summary = summary.copy()
    if "peak_id" in summary.columns and summary["peak_id"].notna().any():
        need = ["gene_id", "peak_id", "primer_index", "primer_type", "gene_strand"]
        if all(col in summary.columns for col in need):
            summary["primer_id"] = (
                summary["gene_id"].astype(str)
                + "|"
                + summary["peak_id"].astype(str)
                + "|idx"
                + summary["primer_index"].astype(str)
                + "|"
                + summary["primer_type"].astype(str)
                + "|strand"
                + summary["gene_strand"].astype(str)
            )
            return summary
    need = ["gene_id", "primer_index", "primer_type", "gene_strand"]
    if all(col in summary.columns for col in need):
        summary["primer_id"] = (
            summary["gene_id"].astype(str)
            + "|idx"
            + summary["primer_index"].astype(str)
            + "|"
            + summary["primer_type"].astype(str)
            + "|strand"
            + summary["gene_strand"].astype(str)
        )
        return summary
    raise ValueError("Cannot construct primer_id from summary columns")


def _empty_output(report: pd.DataFrame, out_path: str, extra_cols: list) -> None:
    cols = list(report.columns) + extra_cols
    for drop in ("alignment_quality", "best_mapq"):
        if drop in cols:
            cols.remove(drop)
    # de-duplicate while preserving order
    seen = set()
    ordered = []
    for c in cols:
        if c not in seen:
            seen.add(c)
            ordered.append(c)
    pd.DataFrame(columns=ordered).to_csv(out_path, sep="\t", index=False)


def main() -> None:
    p = argparse.ArgumentParser(description="Select best primers with 3′-aware specificity filtering")
    p.add_argument("--report", required=True)
    p.add_argument("--summary", required=True)
    p.add_argument("--out", default="best_primers.tsv")
    p.add_argument(
        "--distance_threshold",
        type=float,
        default=1000.0,
        help="Max distance from transcript 3′ end (bp) for an alignment to count "
        "as relevant for on-target usefulness or off-target rejection.",
    )
    p.add_argument(
        "--max_mismatches",
        type=int,
        default=0,
        help="Max mismatches for a relevant hit (0-3). Other genes with a relevant "
        "hit within this mismatch budget reject the primer.",
    )
    args = p.parse_args()

    report = pd.read_csv(args.report, sep="\t", dtype=str)
    extra_cols = [
        "aligned_gene_name",
        "alignments_within_threshold",
        "distance_to_end_min",
        "offtarget_genes",
        "specificity_reason",
    ]

    try:
        summary = pd.read_csv(args.summary, sep="\t", dtype=str)
    except Exception:
        out = report[report.get("num_alignments", "0").astype(int) == 1].copy()
        for drop in ("alignment_quality", "best_mapq"):
            if drop in out.columns:
                out = out.drop(drop, axis=1)
        out.to_csv(args.out, sep="\t", index=False)
        print(f"Wrote {args.out} (fallback quick filter)")
        return

    try:
        summary = _ensure_primer_id(summary)
    except ValueError as exc:
        print(f"{exc}; falling back to report heuristic", file=sys.stderr)
        out = report[report.get("num_alignments", "0").astype(int) == 1].copy()
        for drop in ("alignment_quality", "best_mapq"):
            if drop in out.columns:
                out = out.drop(drop, axis=1)
        out.to_csv(args.out, sep="\t", index=False)
        print(f"Wrote {args.out} (fallback quick filter)")
        return

    print("\n=== 3′-AWARE SPECIFICITY FILTERING ===")
    print(f"Total alignments in summary: {len(summary)}")
    print(f"max_mismatches={args.max_mismatches}  distance_threshold={args.distance_threshold}bp")

    rows = summary.to_dict(orient="records")
    by_primer = group_alignments_by_primer(rows)
    print(f"Unique primers in summary: {len(by_primer)}")

    selected = []
    reject_reasons = {}
    for pid, alns in by_primer.items():
        verdict = classify_primer(alns, args.max_mismatches, args.distance_threshold)
        if not verdict["pass"]:
            reject_reasons[verdict["reason"].split(":")[0]] = (
                reject_reasons.get(verdict["reason"].split(":")[0], 0) + 1
            )
            continue

        rpt = report[report["primer_id"] == pid]
        if rpt.empty:
            reject_reasons["missing_report_row"] = reject_reasons.get("missing_report_row", 0) + 1
            continue

        row = rpt.iloc[0].to_dict()
        for drop in ("alignment_quality", "best_mapq"):
            row.pop(drop, None)
        row["aligned_gene_name"] = verdict["intended_gene"]
        row["alignments_within_threshold"] = verdict["n_relevant"]
        row["distance_to_end_min"] = verdict["distance_to_end_min"]
        row["offtarget_genes"] = ""
        row["specificity_reason"] = "pass"
        selected.append(row)

    print(f"PASS: {len(selected)} primers")
    print("Reject reason counts:")
    for reason, count in sorted(reject_reasons.items(), key=lambda kv: (-kv[1], kv[0])):
        print(f"  {reason}: {count}")

    if not selected:
        _empty_output(report, args.out, extra_cols)
        print(f"Wrote {args.out} (no primers passed specificity filtering)")
        return

    outdf = pd.DataFrame(selected)
    outdf.to_csv(args.out, sep="\t", index=False)
    print(f"Wrote {args.out} with {len(outdf)} primers")
    print("  aligned_gene_name = intended gene (from 0-mismatch hits)")
    print(
        f"  alignments_within_threshold = relevant hits "
        f"(≤{args.max_mismatches} mm and ≤{args.distance_threshold} bp from 3′)"
    )


if __name__ == "__main__":
    main()
