#!/usr/bin/env python3
"""
Reformat genome FASTA sequence IDs to match a BAM file's contig naming.

PeakPrime extracts primer target sequences with Rsamtools::getSeq using
chromosome names from BAM/MACS2/GTF. If the FASTA still uses Ensembl-style
IDs ("19") while peaks use "chr19", sequence extraction fails, e.g.:

  Error: record 1 (chr19:48966583-48966879) failed

This helper only rewrites FASTA header lines (">19 ..." -> ">chr19 ...").
Sequence bases are copied unchanged.

Requires: samtools on PATH (to read the BAM header)

Examples:
  # Check whether FASTA and BAM already match
  bin/harmonize_fasta_to_bam.py --bam sample.bam --fasta genome.fa --check

  # Write a BAM-matched FASTA
  bin/harmonize_fasta_to_bam.py --bam sample.bam --fasta genome.fa \\
      --out genome.bam_chroms.fa

  # Then rebuild the index (required)
  samtools faidx genome.bam_chroms.fa
"""

from __future__ import annotations

import argparse
import os
import sys
from collections import Counter
from pathlib import Path
from typing import Dict, List, Optional, Set

sys.path.insert(0, str(Path(__file__).resolve().parent))
from chrom_harmonize_lib import (  # noqa: E402
    bam_contigs,
    build_mapping,
    open_text,
    parse_fasta_header,
    summarize_style,
)


def collect_fasta_ids(fasta_path: str) -> List[str]:
    ids: List[str] = []
    seen: Set[str] = set()
    with open_text(fasta_path, "r") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            seq_id, _ = parse_fasta_header(line)
            if not seq_id:
                continue
            if seq_id in seen:
                sys.exit(f"Error: duplicate FASTA sequence ID: {seq_id}")
            seen.add(seq_id)
            ids.append(seq_id)
    if not ids:
        sys.exit(f"Error: no FASTA headers found in {fasta_path}")
    return ids


def rewrite_fasta(
    fasta_path: str,
    out_path: str,
    mapping: Dict[str, str],
) -> Counter:
    """Rewrite sequence IDs in FASTA headers; copy sequence lines unchanged."""
    changed = Counter()
    n_headers = 0
    with open_text(fasta_path, "r") as inn, open_text(out_path, "w") as out:
        for line in inn:
            if not line.startswith(">"):
                out.write(line)
                continue
            n_headers += 1
            old_id, rest = parse_fasta_header(line)
            new_id = mapping.get(old_id, old_id)
            if new_id != old_id:
                changed[f"{old_id}->{new_id}"] += 1
            out.write(f">{new_id}{rest}")
    if n_headers == 0:
        sys.exit(f"Error: no FASTA headers written from {fasta_path}")
    return changed


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Rewrite FASTA sequence IDs to match BAM contig names.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("--bam", required=True, help="BAM/CRAM/SAM whose contig names are authoritative")
    parser.add_argument("--fasta", required=True, help="Input genome FASTA (.fa/.fasta, optional .gz)")
    parser.add_argument(
        "--out",
        help="Output FASTA path. Required unless --check is set. Prefer uncompressed .fa for PeakPrime.",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="Only report chromosome compatibility; do not write an output FASTA",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Exit with error if any FASTA sequence ID cannot be mapped uniquely to the BAM",
    )
    args = parser.parse_args(argv)

    if not args.check and not args.out:
        parser.error("--out is required unless --check is used")
    if not os.path.exists(args.bam):
        sys.exit(f"Error: BAM not found: {args.bam}")
    if not os.path.exists(args.fasta):
        sys.exit(f"Error: FASTA not found: {args.fasta}")

    print("============================================================")
    print("HARMONIZE FASTA SEQUENCE IDS TO BAM")
    print("============================================================")

    bam_chroms = bam_contigs(args.bam)
    fasta_ids = collect_fasta_ids(args.fasta)
    mapping, unmapped, ambiguous = build_mapping(fasta_ids, bam_chroms)

    needs_change = sum(1 for src, dst in mapping.items() if src != dst)
    already_ok = sum(1 for src, dst in mapping.items() if src == dst)

    print(f"BAM contigs : {len(bam_chroms)}  [{summarize_style(bam_chroms)}]")
    print(f"FASTA seqs  : {len(fasta_ids)}  [{summarize_style(fasta_ids)}]")
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
        print(f"\nWarning: no BAM match for FASTA id(s): {preview}{more}")
        print("  These IDs will be left unchanged.")

    if ambiguous:
        preview = ", ".join(ambiguous[:10])
        print(f"\nWarning: ambiguous BAM matches for FASTA id(s): {preview}")
        print("  These IDs will be left unchanged.")

    if args.strict and (unmapped or ambiguous):
        sys.exit(
            "Error: --strict set and some FASTA sequence IDs could not be mapped uniquely."
        )

    if args.check:
        if needs_change == 0 and not unmapped and not ambiguous:
            print("\n✓ FASTA sequence IDs already match the BAM.")
            return 0
        if needs_change:
            print("\n✗ FASTA needs reformatting to match the BAM.")
            print("  Re-run without --check and with --out to write a fixed FASTA.")
            return 2
        print("\n⚠ Some FASTA IDs have no/ambiguous BAM match; review warnings above.")
        return 3

    if os.path.abspath(args.fasta) == os.path.abspath(args.out):
        sys.exit("Error: --out must differ from --fasta (refusing to overwrite input in-place)")

    print("\nRewriting FASTA headers (sequence bases copied unchanged)...")
    changed = rewrite_fasta(args.fasta, args.out, mapping)

    print("\n============================================================")
    print(f"✓ Wrote: {args.out}")
    print(f"✓ Headers renamed: {sum(changed.values())} / {len(fasta_ids)}")
    if changed:
        print("Rename counts:")
        for rule, count in changed.most_common(30):
            print(f"  {rule}: {count}")
    print("============================================================")
    print("Next steps:")
    print(f"  samtools faidx {args.out}")
    print("  Point PeakPrime --fasta (or params.fasta) at this new file.")
    print("  Keep using your BAM-matched GTF (genes_bam_chroms.gtf).")
    return 0


if __name__ == "__main__":
    sys.exit(main())
