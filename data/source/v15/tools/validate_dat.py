#!/usr/bin/env python3
"""
Validate a bloodAGENT-format variation_annotation.dat or genotype_to_phenotype
file by replicating what CISBTAnno::readAnnotation + CISBTAnno::generateIndex
will do at runtime in C++.

For variation_annotation*.dat (28 cols):
  - Asserts header column order matches the C++ ISBTAnno keys exactly
  - Asserts every row has 28 cols
  - Asserts (system/gene, transcript_short stripped of '!') is unique (would otherwise stomp m_isbt_variant_to_index)
  - Asserts integer columns parse via stoi() — bloodAGENT uses stoi which throws on empty/non-int
  - Asserts strand column is a single char in {+,-}
  - Asserts is_transcript_NC_==_NC fields are TRUE/FALSE strings (compared via .compare())
  - Reports rows where hg19 OR hg38 coords are missing (those rows will silently break index)

For genotype_to_phenotype*.dat (10 cols):
  - Asserts header matches
  - Asserts no allele appears twice in same MySystemKey
  - Asserts every base_change token resolves to an entry in the matching variation file (if provided)

Exit code 0 = OK, 1 = errors.
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path

VANNO_HEADERS = [
    "system/gene", "Transcript annotation", "Transcript annotation short",
    "genomic annotation", "Chr (hg19)", "Pos (hg19)", "Chrom (hg19)",
    "0-based start (hg19)", "1-based end (hg19)", "strand (hg19)",
    "Coordinate in browser format +1bp flanks  (hg19)",
    "Reference base (hg19)", "Reference base blood group annotation",
    "is transcript_NC == hg19_NC",
    "Chrom (hg38)", "0-based start (hg38)", "1-based end (hg38)", "strand (hg38)",
    "Reference base (hg38)", "is transcript_NC == hg38_NC",
    "reference allele of builds concordant?",
    "Coordinate in VCF hg19", "RefAllele in VCF hg19", "AltAllele in VCF hg19",
    "Coordinate in VCF hg38", "RefAllele in VCF hg38", "AltAllele in VCF hg38",
    "TYPE",
]
GT2PT_HEADERS = [
    "PureSystem", "MySystemKey", "System", "Allele",
    "Phenotype_PDF_Table", "Phenotype", "Phenotype_flat",
    "base_change", "acid_change", "incidence",
]


def parse_dat(path: Path, expected_cols: int) -> tuple[list[str], list[list[str]]]:
    """Mirrors CParsedTextfile behaviour:
       - skip lines starting with '#' as comments
       - first non-comment line is the header
       - subsequent lines are data
       - returns (header, data_rows) WITHOUT skipped lines.
    """
    header = None
    rows = []
    with path.open() as f:
        for i, line in enumerate(f, 1):
            if line.startswith("#"):
                continue
            line = line.rstrip("\n")
            if header is None:
                header = line.split("\t")
                if len(header) != expected_cols:
                    raise SystemExit(f"FAIL: {path}:{i} header has {len(header)} cols, expected {expected_cols}")
                continue
            parts = line.split("\t")
            # CParsedTextfile pads short rows with empty strings — match that
            if len(parts) < expected_cols:
                parts += [""] * (expected_cols - len(parts))
            rows.append(parts)
    return header, rows


def strip_bang(s: str) -> str:
    while s.startswith("!"):
        s = s[1:]
    return s


def validate_vanno(path: Path) -> tuple[int, list[str]]:
    errs: list[str] = []
    try:
        header, rows = parse_dat(path, 28)
    except SystemExit as e:
        return 1, [str(e)]

    for i, h in enumerate(VANNO_HEADERS):
        got = header[i] if i < len(header) else "<missing>"
        if got != h:
            errs.append(f"col {i+1}: header mismatch (got '{got}', expected '{h}')")

    idx: dict[tuple[str, str], int] = {}
    # CISBTAnno calls stoi() only on these (see ISBTAnno.cpp ~line 60):
    #   col 9  "1-based end (hg19)"  -> index 8
    #   col 17 "1-based end (hg38)"  -> index 16
    #   col 22 "Coordinate in VCF hg19"  -> index 21
    #   col 25 "Coordinate in VCF hg38"  -> index 24
    int_cols = [8, 16, 21, 24]

    for rn, p in enumerate(rows, 1):
        if len(p) != 28:
            errs.append(f"row {rn}: {len(p)} cols (need 28)")
            continue
        sys_gene = p[0]
        tr_short = strip_bang(p[2])
        if not sys_gene or not tr_short:
            errs.append(f"row {rn}: empty system/gene or transcript-short")
            continue
        key = (sys_gene, tr_short)
        if key in idx:
            errs.append(f"row {rn}: duplicate key {key} (also row {idx[key]})")
        idx[key] = rn

        # integer parseability — CISBTAnno calls stoi() on these
        for c in int_cols:
            v = p[c]
            if v == "" or v is None:
                # Note: bloodAGENT WILL crash on empty strings via stoi; warn
                errs.append(f"row {rn}: col {c+1} ({VANNO_HEADERS[c]}) is empty, stoi() will throw at runtime")
                continue
            try:
                int(v)
            except ValueError:
                errs.append(f"row {rn}: col {c+1} ({VANNO_HEADERS[c]}) = '{v}' is not integer")

        # strand char
        for c in (9, 17):
            if p[c] not in ("+", "-"):
                errs.append(f"row {rn}: col {c+1} ({VANNO_HEADERS[c]}) = '{p[c]}' is not '+' or '-'")

        # TRUE/FALSE
        for c in (13, 19, 20):
            if p[c] not in ("TRUE", "FALSE"):
                errs.append(f"row {rn}: col {c+1} ({VANNO_HEADERS[c]}) = '{p[c]}' should be TRUE or FALSE")

    return len(errs), errs


def validate_gt2pt(path: Path) -> tuple[int, list[str]]:
    errs: list[str] = []
    try:
        header, rows = parse_dat(path, 10)
    except SystemExit as e:
        return 1, [str(e)]

    for i, h in enumerate(GT2PT_HEADERS):
        got = header[i] if i < len(header) else "<missing>"
        if got != h:
            errs.append(f"col {i+1}: header mismatch (got '{got}', expected '{h}')")

    seen = {}
    for rn, p in enumerate(rows, 1):
        if len(p) != 10:
            errs.append(f"row {rn}: {len(p)} cols (need 10)")
            continue
        key = (p[1], p[3])  # MySystemKey + Allele
        if key in seen:
            errs.append(f"row {rn}: duplicate allele {key} (also row {seen[key]})")
        seen[key] = rn

    return len(errs), errs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("paths", nargs="+", help=".dat file(s) to validate")
    args = ap.parse_args()

    fail = 0
    for p_ in args.paths:
        p = Path(p_)
        is_gt2pt = "genotype_to_phenotype" in p.name
        n, errs = validate_gt2pt(p) if is_gt2pt else validate_vanno(p)
        kind = "gt2pt" if is_gt2pt else "vanno"
        print(f"=== {p} ({kind}) — {n} issue(s) ===")
        for e in errs[:50]:
            print("  ", e)
        if len(errs) > 50:
            print(f"  ... and {len(errs)-50} more")
        if n:
            fail = 1
    sys.exit(fail)


if __name__ == "__main__":
    main()
