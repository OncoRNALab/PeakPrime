#!/usr/bin/env python3
"""Shared chromosome-name helpers for PeakPrime BAM/GTF/FASTA harmonization."""

from __future__ import annotations

import gzip
import subprocess
import sys
from typing import Dict, Iterable, List, TextIO, Tuple

# Canonical aliases for mitochondrial contig names across NCBI / Ensembl / UCSC
_MT_ALIASES = {"m", "mt", "chrm", "chrmt"}


def open_text(path: str, mode: str = "r") -> TextIO:
    if path == "-":
        return sys.stdin if "r" in mode else sys.stdout
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t")  # type: ignore[return-value]
    return open(path, mode, encoding="utf-8")


def canonical_chrom(name: str) -> str:
    """Collapse UCSC/Ensembl/NCBI primary contig aliases to one key."""
    key = name.strip()
    if key.lower().startswith("chr"):
        key = key[3:]
    key_l = key.lower()
    if key_l in _MT_ALIASES:
        return "MT"
    return key_l


def bam_contigs(bam_path: str) -> List[str]:
    """Read contig names from BAM/CRAM/SAM header via samtools."""
    try:
        result = subprocess.run(
            ["samtools", "view", "-H", bam_path],
            check=True,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError:
        sys.exit("Error: samtools not found on PATH. Install samtools to use this helper.")
    except subprocess.CalledProcessError as exc:
        sys.exit(
            f"Error: failed to read BAM header from {bam_path}\n"
            f"{exc.stderr.strip() or exc}"
        )

    contigs: List[str] = []
    for line in result.stdout.splitlines():
        if not line.startswith("@SQ"):
            continue
        for field in line.split("\t")[1:]:
            if field.startswith("SN:"):
                contigs.append(field[3:])
                break
    if not contigs:
        sys.exit(f"Error: no @SQ contigs found in BAM header: {bam_path}")
    return contigs


def index_by_canonical(names: Iterable[str]) -> Dict[str, List[str]]:
    index: Dict[str, List[str]] = {}
    for name in names:
        index.setdefault(canonical_chrom(name), []).append(name)
    return index


def build_mapping(
    source_chroms: List[str],
    bam_chroms: List[str],
) -> Tuple[Dict[str, str], List[str], List[str]]:
    """
    Map each source chrom/seq id to a BAM contig.

    Returns: (mapping, unmapped, ambiguous)
    """
    bam_exact = set(bam_chroms)
    bam_by_key = index_by_canonical(bam_chroms)

    mapping: Dict[str, str] = {}
    unmapped: List[str] = []
    ambiguous: List[str] = []

    for chrom in source_chroms:
        if chrom in bam_exact:
            mapping[chrom] = chrom
            continue

        candidates = bam_by_key.get(canonical_chrom(chrom), [])
        unique = list(dict.fromkeys(candidates))
        if len(unique) == 1:
            mapping[chrom] = unique[0]
        elif len(unique) > 1:
            ambiguous.append(chrom)
        else:
            unmapped.append(chrom)

    return mapping, unmapped, ambiguous


def summarize_style(names: Iterable[str]) -> str:
    names = list(names)
    if not names:
        return "empty"
    with_chr = sum(1 for n in names if n.startswith("chr"))
    if with_chr == len(names):
        return "UCSC-like (all 'chr' prefix)"
    if with_chr == 0:
        return "Ensembl/NCBI-like (no 'chr' prefix)"
    return f"mixed ({with_chr}/{len(names)} have 'chr' prefix)"


def parse_fasta_header(line: str) -> Tuple[str, str]:
    """
    Split a FASTA header into (seq_id, remainder_including_leading_whitespace).

    '>19 dna:chromosome ...' -> ('19', ' dna:chromosome ...')
    '>chr19'                 -> ('chr19', '')
    """
    if not line.startswith(">"):
        raise ValueError(f"Not a FASTA header: {line[:40]!r}")
    body = line[1:].rstrip("\n")
    if not body:
        return "", "\n" if line.endswith("\n") else ""
    parts = body.split(None, 1)
    seq_id = parts[0]
    rest = (" " + parts[1]) if len(parts) > 1 else ""
    newline = "\n" if line.endswith("\n") else ""
    return seq_id, rest + newline
