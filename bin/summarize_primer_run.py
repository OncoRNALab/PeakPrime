#!/usr/bin/env python3
"""
Generate a PeakPrime primer-run summary (always published by the pipeline).

Outputs (in --outdir):
  gene_funnel.tsv              — per-gene stage / counts
  primer_specificity.tsv       — per-primer specificity audit (if alignments present)
  final_specific_primers.tsv   — primers that pass shared specificity policy
  run_summary.txt              — human-readable funnel + notes
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from primer_specificity import classify_primer, group_alignments_by_primer  # noqa: E402


def _exists(path: str | None) -> bool:
    return bool(path) and path != "NO_FILE" and Path(path).exists() and Path(path).stat().st_size > 0


def _read_tsv(path: str) -> list:
    with open(path, newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def _read_genes(path: str) -> list:
    return [g.strip() for g in Path(path).read_text().splitlines() if g.strip()]


def main() -> int:
    p = argparse.ArgumentParser(description="Summarize a PeakPrime primer design run")
    p.add_argument("--genes", required=True, help="Target gene list")
    p.add_argument("--qc", default="NO_FILE", help="peaks_qc_summary.tsv")
    p.add_argument("--peaks", default="NO_FILE", help="selected_peaks.tsv")
    p.add_argument("--cdna_primers", default="NO_FILE", help="cdna_primers.tsv")
    p.add_argument("--best_primers", default="NO_FILE", help="best_primers.tsv")
    p.add_argument("--optimized_primers", default="NO_FILE", help="optimized_primers*.tsv")
    p.add_argument("--alignment_summary", default="NO_FILE", help="primer_alignment_summary.tsv")
    p.add_argument("--max_mismatches", type=int, default=3)
    p.add_argument("--distance_threshold", type=float, default=1000.0)
    p.add_argument("--outdir", required=True, help="Directory for summary outputs")
    args = p.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    genes = _read_genes(args.genes)
    qc_by = {}
    if _exists(args.qc):
        for r in _read_tsv(args.qc):
            qc_by[r["gene_id"]] = r

    peaks_n = Counter()
    if _exists(args.peaks):
        for r in _read_tsv(args.peaks):
            peaks_n[r.get("gene") or r.get("gene_id")] += 1

    cdna_n = Counter()
    if _exists(args.cdna_primers):
        for r in _read_tsv(args.cdna_primers):
            cdna_n[r["gene_id"]] += 1

    best_n = Counter()
    best_rows = []
    if _exists(args.best_primers):
        best_rows = _read_tsv(args.best_primers)
        for r in best_rows:
            best_n[r["gene_id"]] += 1

    opt_n = Counter()
    opt_rows = []
    if _exists(args.optimized_primers):
        opt_rows = _read_tsv(args.optimized_primers)
        for r in opt_rows:
            opt_n[r["gene_id"]] += 1

    # Specificity audit on optimized (fallback: best)
    primer_source = opt_rows if opt_rows else best_rows
    source_label = "optimized" if opt_rows else "best"
    by_primer_aln = {}
    if _exists(args.alignment_summary):
        aln_rows = _read_tsv(args.alignment_summary)
        by_primer_aln = group_alignments_by_primer(aln_rows)

    specificity_rows = []
    final_specific = []
    for r in primer_source:
        pid = r.get("primer_id")
        if not pid:
            continue
        alns = by_primer_aln.get(pid, [])
        if alns:
            verdict = classify_primer(alns, args.max_mismatches, args.distance_threshold)
        else:
            # No alignment data: cannot audit; mark unknown
            verdict = {
                "pass": False,
                "reason": "no_alignment_data",
                "intended_gene": r.get("aligned_gene_name"),
                "offtarget_genes": [],
                "distance_to_end_min": r.get("distance_to_end_min"),
                "n_relevant": None,
            }
        rec = {
            "gene_id": r.get("gene_id"),
            "primer_id": pid,
            "peak_id": r.get("peak_id", ""),
            "primer_sequence": r.get("primer_sequence", ""),
            "primer_type": r.get("primer_type", ""),
            "source": source_label,
            "intended_gene": verdict.get("intended_gene"),
            "pass": verdict.get("pass"),
            "reason": verdict.get("reason"),
            "offtarget_genes": ",".join(verdict.get("offtarget_genes") or []),
            "distance_to_end_min": verdict.get("distance_to_end_min"),
            "n_relevant": verdict.get("n_relevant"),
        }
        specificity_rows.append(rec)
        if verdict.get("pass"):
            final_specific.append(rec)

    # Gene funnel
    funnel = []
    stage_counts = Counter()
    for g in genes:
        q = qc_by.get(g, {})
        if opt_n[g]:
            stage = "optimized_primers"
        elif best_n[g]:
            stage = "best_primers_only"
        elif cdna_n[g]:
            stage = "primer3_no_best"
        elif peaks_n[g] or str(q.get("final_selection", "")).upper() == "TRUE":
            stage = "peaks_no_primers"
        elif str(q.get("has_macs2_peaks", "")).upper() == "FALSE":
            stage = "no_macs2_peaks"
        else:
            stage = q.get("failure_reason") or "unknown"
        stage_counts[stage] += 1
        funnel.append(
            {
                "gene_id": g,
                "fail_or_final_stage": stage,
                "qc_failure_reason": q.get("failure_reason", ""),
                "has_macs2_peaks": q.get("has_macs2_peaks", ""),
                "final_selection": q.get("final_selection", ""),
                "n_selected_peaks": peaks_n[g],
                "n_cdna_primers": cdna_n[g],
                "n_best_primers": best_n[g],
                "n_optimized_primers": opt_n[g],
                "n_final_specific_primers": sum(1 for r in final_specific if r["gene_id"] == g),
            }
        )

    def write_tsv(path: Path, rows: list, fieldnames: list) -> None:
        with open(path, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
            w.writeheader()
            for row in rows:
                w.writerow(row)

    write_tsv(
        outdir / "gene_funnel.tsv",
        funnel,
        [
            "gene_id",
            "fail_or_final_stage",
            "qc_failure_reason",
            "has_macs2_peaks",
            "final_selection",
            "n_selected_peaks",
            "n_cdna_primers",
            "n_best_primers",
            "n_optimized_primers",
            "n_final_specific_primers",
        ],
    )
    write_tsv(
        outdir / "primer_specificity.tsv",
        specificity_rows,
        [
            "gene_id",
            "primer_id",
            "peak_id",
            "primer_type",
            "primer_sequence",
            "source",
            "intended_gene",
            "pass",
            "reason",
            "offtarget_genes",
            "distance_to_end_min",
            "n_relevant",
        ],
    )
    write_tsv(
        outdir / "final_specific_primers.tsv",
        final_specific,
        [
            "gene_id",
            "intended_gene",
            "primer_id",
            "peak_id",
            "primer_type",
            "primer_sequence",
            "distance_to_end_min",
            "source",
        ],
    )

    lines = []
    lines.append("PeakPrime primer run summary")
    lines.append("=" * 60)
    lines.append(f"Target genes: {len(genes)}")
    lines.append(f"Specificity: relevant = mismatches≤{args.max_mismatches} AND distance≤{args.distance_threshold} bp")
    lines.append(f"  Reject if any OTHER gene has a relevant hit; require intended gene near 3′.")
    lines.append("")
    lines.append("Gene-level stages:")
    for stage, n in sorted(stage_counts.items(), key=lambda kv: (-kv[1], kv[0])):
        lines.append(f"  {stage}: {n}")
    lines.append("")
    lines.append(f"Primers audited ({source_label}): {len(specificity_rows)}")
    lines.append(f"Final specific primers: {len(final_specific)} across {len({r['gene_id'] for r in final_specific})} genes")
    if specificity_rows:
        reasons = Counter(r["reason"] for r in specificity_rows if not r["pass"])
        if reasons:
            lines.append("Rejected primer reasons:")
            for reason, n in reasons.most_common():
                lines.append(f"  {reason}: {n}")
    lines.append("")
    lines.append("Outputs:")
    lines.append("  gene_funnel.tsv")
    lines.append("  primer_specificity.tsv")
    lines.append("  final_specific_primers.tsv")
    lines.append("  run_summary.txt")
    text = "\n".join(lines) + "\n"
    (outdir / "run_summary.txt").write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
