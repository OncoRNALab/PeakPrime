#!/usr/bin/env python3
"""
Reformat GTF chromosome names to match a BAM file's contig naming.

PeakPrime joins MACS2 peaks (BAM contig names) to genes from the GTF using
exact chromosome string equality. Ensembl-style GTFs often use "19" while
BAM/MACS2 use "chr19" (or vice versa), which yields "No MACS2 peaks detected"
for every gene even when peaks exist.

This helper is intentionally opt-in: inspect the BAM header, rewrite the GTF
seqname column to the matching BAM contig names, and use the rewritten GTF
as pipeline input.

Requires: samtools on PATH

Examples:
  # Check whether GTF and BAM chrom styles already match
  bin/harmonize_gtf_to_bam.py --bam sample.bam --gtf genes.gtf --check

  # Write a BAM-matched GTF
  bin/harmonize_gtf_to_bam.py --bam sample.bam --gtf genes.gtf --out genes.bam_chroms.gtf

  # Fail if any GTF chromosome cannot be mapped uniquely
  bin/harmonize_gtf_to_bam.py --bam sample.bam --gtf genes.gtf --out genes.bam_chroms.gtf --strict
"""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Set


def _load_chrom_lib():
    """Load shared helper from the same bin/ directory (works without PYTHONPATH)."""
    lib_path = Path(__file__).resolve().parent / "chrom_harmonize_lib.py"
    if not lib_path.is_file():
        sys.exit(
            f"Error: missing shared helper: {lib_path}\n"
            "Copy or sync bin/chrom_harmonize_lib.py next to this script "
            "(same PeakPrime checkout / git pull)."
        )
    spec = importlib.util.spec_from_file_location("chrom_harmonize_lib", lib_path)
    if spec is None or spec.loader is None:
        sys.exit(f"Error: could not load {lib_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


_chrom = _load_chrom_lib()
bam_contigs = _chrom.bam_contigs
build_mapping = _chrom.build_mapping
open_text = _chrom.open_text
summarize_style = _chrom.summarize_style


def collect_gtf_chroms(gtf_path: str) -> List[str]:
    chroms: Set[str] = set()
    with open_text(gtf_path, "r") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            chrom = line.split("\t", 1)[0]
            if chrom:
                chroms.add(chrom)
    return sorted(chroms)


def rewrite_gtf(
    gtf_path: str,
    out_path: str,
    mapping: Dict[str, str],
) -> Counter:
    changed = Counter()
    with open_text(gtf_path, "r") as inn, open_text(out_path, "w") as out:
        for line in inn:
            if not line or line.startswith("#"):
                out.write(line)
                continue
            parts = line.split("\t")
            if len(parts) < 9:
                out.write(line)
                continue
            old = parts[0]
            new = mapping.get(old, old)
            if new != old:
                changed[f"{old}->{new}"] += 1
                parts[0] = new
                out.write("\t".join(parts))
            else:
                out.write(line)
    return changed


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Rewrite GTF chromosome names to match BAM contig names.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--bam", required=True, help="BAM/CRAM/SAM whose contig names are authoritative")
    parser.add_argument("--gtf", required=True, help="Input GTF (.gtf or .gtf.gz)")
    parser.add_argument(
        "--out",
        help="Output GTF path (.gtf or .gtf.gz). Required unless --check is set.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Only report chromosome compatibility; do not write an output GTF",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit with error if any GTF chromosome cannot be mapped uniquely to the BAM",
    )
    args = parser.parse_args(argv)

    if not args.check and not args.out:
        parser.error("--out is required unless --check is used")
    if not os.path.exists(args.bam):
        sys.exit(f"Error: BAM not found: {args.bam}")
    if not os.path.exists(args.gtf):
        sys.exit(f"Error: GTF not found: {args.gtf}")

    print("============================================================")
    print("HARMONIZE GTF CHROMOSOMES TO BAM")
    print("============================================================")

    bam_chroms = bam_contigs(args.bam)
    gtf_chroms = collect_gtf_chroms(args.gtf)
    mapping, unmapped, ambiguous = build_mapping(gtf_chroms, bam_chroms)

    needs_change = sum(1 for src, dst in mapping.items() if src != dst)
    already_ok = sum(1 for src, dst in mapping.items() if src == dst)

    print(f"BAM contigs : {len(bam_chroms)}  [{summarize_style(bam_chroms)}]")
    print(f"GTF chroms  : {len(gtf_chroms)}  [{summarize_style(gtf_chroms)}]")
    print(f"Exact match : {already_ok}")
    print(f"Will rename : {needs_change}")
    print(f"Unmapped    : {len(unmapped)}")
    print(f"Ambiguous   : {len(ambiguous)}")

    if needs_change:
        print("\nSample renames:")
        shown = 0
        for src, dst in sorted(mapping.items(), key=lambda kv: (kv[0] == kv[1], kv[0])):
            if src == dst:
                continue
            print(f"  {src} -> {dst}")
            shown += 1
            if shown >= 15:
                remaining = needs_change - shown
                if remaining > 0:
                    print(f"  ... and {remaining} more")
                break

    if unmapped:
        preview = ", ".join(unmapped[:10])
        more = "" if len(unmapped) <= 10 else f" (+{len(unmapped) - 10} more)"
        print(f"\nWarning: no BAM match for GTF chrom(s): {preview}{more}")
        print("  These names will be left unchanged.")

    if ambiguous:
        preview = ", ".join(ambiguous[:10])
        print(f"\nWarning: ambiguous BAM matches for GTF chrom(s): {preview}")
        print("  These names will be left unchanged.")

    if args.strict and (unmapped or ambiguous):
        sys.exit(
            "Error: --strict set and some GTF chromosomes could not be mapped uniquely."
        )

    if args.check:
        if needs_change == 0 and not unmapped and not ambiguous:
            print("\n✓ GTF chromosome names already match the BAM.")
            return 0
        if needs_change:
            print("\n✗ GTF needs reformatting to match the BAM.")
            print("  Re-run without --check and with --out to write a fixed GTF.")
            return 2
        print("\n⚠ Some GTF chromosomes have no/ambiguous BAM match; review warnings above.")
        return 3

    changed = rewrite_gtf(args.gtf, args.out, mapping)
    total_rows = sum(changed.values())
    print("\n============================================================")
    print(f"✓ Wrote: {args.out}")
    print(f"✓ Feature rows renamed: {total_rows}")
    if changed:
        print("Rename counts:")
        for rule, count in changed.most_common(20):
            print(f"  {rule}: {count}")
    print("============================================================")
    print(
        "Note: genome FASTA used with PeakPrime should use the same contig names "
        "as the BAM. See bin/harmonize_fasta_to_bam.py"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
