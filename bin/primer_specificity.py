#!/usr/bin/env python3
"""
Shared primer specificity helpers for PeakPrime.

Policy (3′-end protocols):
  An alignment is *relevant* if mismatches ≤ max_mismatches AND
  distance_to_end ≤ distance_threshold.

  A primer PASSES when:
    1. It has ≥1 zero-mismatch alignment (defines the intended gene symbol)
    2. Zero-mismatch hits are to a single gene
    3. Among relevant alignments, that intended gene is present
    4. Among relevant alignments, no OTHER gene is present

  Off-targets far from the transcript 3′ end are ignored (not sequenced in
  3′-end protocols). Off-targets near the 3′ end with ≤ max_mismatches reject
  the primer — even if the on-target hit itself is far from the 3′ end.
"""

from __future__ import annotations

from typing import Any, Dict, Iterable, List, Optional, Sequence, Set, Tuple


def _as_int(value: Any, default: int = 9999) -> int:
    try:
        if value is None or value == "" or str(value).upper() == "NA":
            return default
        return int(float(value))
    except (TypeError, ValueError):
        return default


def _as_float(value: Any, default: float = float("inf")) -> float:
    try:
        if value is None or value == "" or str(value).upper() == "NA":
            return default
        return float(value)
    except (TypeError, ValueError):
        return default


def intended_gene_from_alignments(alignments: Sequence[Dict[str, Any]]) -> Tuple[Optional[str], str]:
    """
    Infer intended gene symbol from 0-mismatch alignments.

    Returns (gene_symbol_or_None, reason_if_none).
    """
    zero = []
    for row in alignments:
        if _as_int(row.get("mismatches")) != 0:
            continue
        gene = row.get("aligned_gene_name")
        if gene is None or str(gene).strip() == "" or str(gene).upper() == "NA":
            continue
        zero.append(str(gene))
    if not zero:
        return None, "no_zero_mismatch_alignment"
    genes = sorted(set(zero))
    if len(genes) > 1:
        return None, f"zero_mismatch_multi_gene:{','.join(genes)}"
    return genes[0], "ok"


def relevant_alignments(
    alignments: Sequence[Dict[str, Any]],
    max_mismatches: int,
    distance_threshold: float,
) -> List[Dict[str, Any]]:
    """Alignments with mismatches ≤ max_mismatches and distance ≤ threshold."""
    out = []
    for row in alignments:
        mm = _as_int(row.get("mismatches"))
        dist = _as_float(row.get("distance_to_end"))
        if mm <= max_mismatches and dist <= distance_threshold:
            out.append(row)
    return out


def genes_in_alignments(alignments: Iterable[Dict[str, Any]]) -> List[str]:
    genes: Set[str] = set()
    for row in alignments:
        gene = row.get("aligned_gene_name")
        if gene is None or str(gene).strip() == "" or str(gene).upper() == "NA":
            continue
        genes.add(str(gene))
    return sorted(genes)


def classify_primer(
    alignments: Sequence[Dict[str, Any]],
    max_mismatches: int,
    distance_threshold: float,
    intended_gene: Optional[str] = None,
) -> Dict[str, Any]:
    """
    Classify one primer's transcriptome alignments.

    If intended_gene is None, it is inferred from 0-mismatch hits.
    """
    result: Dict[str, Any] = {
        "pass": False,
        "reason": "",
        "intended_gene": None,
        "relevant_genes": [],
        "offtarget_genes": [],
        "n_relevant": 0,
        "distance_to_end_min": None,
        "n_zero_mismatch": 0,
    }

    if intended_gene is None:
        intended_gene, reason = intended_gene_from_alignments(alignments)
        if intended_gene is None:
            result["reason"] = reason
            return result
    result["intended_gene"] = intended_gene
    result["n_zero_mismatch"] = sum(
        1
        for row in alignments
        if _as_int(row.get("mismatches")) == 0
        and str(row.get("aligned_gene_name", "")) == intended_gene
    )

    relevant = relevant_alignments(alignments, max_mismatches, distance_threshold)
    result["n_relevant"] = len(relevant)
    relevant_genes = genes_in_alignments(relevant)
    result["relevant_genes"] = relevant_genes

    if not relevant:
        result["reason"] = "no_relevant_alignment_near_3prime"
        return result

    offtargets = [g for g in relevant_genes if g != intended_gene]
    result["offtarget_genes"] = offtargets

    if offtargets:
        result["reason"] = f"offtarget_near_3prime:{','.join(offtargets)}"
        return result

    if intended_gene not in relevant_genes:
        result["reason"] = "intended_gene_not_near_3prime"
        return result

    # Min distance among relevant on-target hits
    ontarget_dists = [
        _as_float(row.get("distance_to_end"))
        for row in relevant
        if str(row.get("aligned_gene_name", "")) == intended_gene
    ]
    result["distance_to_end_min"] = int(min(ontarget_dists)) if ontarget_dists else None
    result["pass"] = True
    result["reason"] = "pass"
    return result


def group_alignments_by_primer(
    rows: Iterable[Dict[str, Any]],
    primer_id_col: str = "primer_id",
) -> Dict[str, List[Dict[str, Any]]]:
    grouped: Dict[str, List[Dict[str, Any]]] = {}
    for row in rows:
        pid = row.get(primer_id_col)
        if pid is None or str(pid).strip() == "":
            continue
        grouped.setdefault(str(pid), []).append(row)
    return grouped
