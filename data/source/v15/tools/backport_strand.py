#!/usr/bin/env python3
"""
Backport the `strand (hg19)` and `strand (hg38)` columns into the V15-derived
variation_annotation.v15.dat, using two sources:
  1. The existing bloodAGENT master (data/config/variation_annotation.dat)
     -- canonical: same project, same hg coords convention.
  2. A small hand-curated dictionary of RefSeq-known strands for genes that
     don't appear in #1 (i.e., the 12 new V15 systems).

After this, the V15 vanno has every row's strand populated (was '+' default
from build_variation_annotation.py — wrong for many genes).

Run after build_variation_annotation.py. Idempotent.
"""
from __future__ import annotations
import argparse, sys
from pathlib import Path

# Hand-curated strand for V15 genes NOT present in the old master.
# Source: NCBI Gene / RefSeq (verified against the corresponding hg38 entries).
EXTRA_STRAND = {
    # New ISBT V15 systems (12 missing)
    "C4A":     "+",   # chr6, C4A
    "C4B":     "+",   # chr6, C4B (CH/RG)
    "PRNP":    "+",   # chr20, KANNO
    "B4GALNT2":"-",   # chr17, SID
    "SLC44A2": "+",   # chr19, CTL2
    "ABCC4":   "-",   # chr13, PEL
    "EMP3":    "+",   # chr19, MAM
    "PIGG":    "+",   # chr4,  EMM
    "ABCC1":   "+",   # chr16, ABCC1
    "PIEZO1":  "-",   # chr16, ER
    "CD36":    "-",   # chr7,  CD36
    "ATP11C":  "-",   # chrX,  ATP11C
    "MAL":     "-",   # chr2,  MAL
    # Other gene-name keys that appear in V15 but not in old master
    "MNS_HYBRID": "-",  # GYPA/GYPB hybrids on chr4 minus strand
    "RSRP1":   "+",   # Arginine and serine rich protein 1 — appears as gene_name for one variant
    "CD99":    "-",   # chrX/Y (PAR1) — used by XG system
}


def load_old_strand(old_master_path: Path) -> dict[str, str]:
    """Return {system_or_gene_key -> strand} from existing variation_annotation.dat."""
    result: dict[str, str] = {}
    with old_master_path.open() as f:
        header = next(f).rstrip("\n").split("\t")
        try:
            sys_col = header.index("system/gene")
            strand_col = header.index("strand (hg19)")
        except ValueError as e:
            raise SystemExit(f"FAIL: header missing expected column: {e}")
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) <= max(sys_col, strand_col):
                continue
            key, strand = p[sys_col], p[strand_col]
            if strand in ("+", "-") and key and key not in result:
                result[key] = strand
    return result


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--vanno", required=True, help="V15-derived variation_annotation.v15.dat to update IN-PLACE")
    ap.add_argument("--old-master", default="data/config/variation_annotation.dat")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    vanno = Path(args.vanno)
    old = Path(args.old_master)
    old_strand = load_old_strand(old)
    # ISBT V15 uses gene names; some old keys are system symbols. Build alias map.
    # P1PK -> A4GALT, LU -> BCAM, etc. Hardcoded because there are <20 such cases.
    SYM_TO_GENE = {
        "P1PK":"A4GALT", "LU":"BCAM", "LE":"FUT3", "FY":"ACKR1", "JK":"SLC14A1",
        "DI":"SLC4A1", "YT":"ACHE", "SC":"ERMAP", "DO":"ART4", "CO":"AQP1",
        "LW":"ICAM4", "H":"FUT1", "GE":"GYPC", "CROM":"CD55", "KN":"CR1",
        "IN":"CD44", "OK":"BSG", "RAPH":"CD151", "JMH":"SEMA7A", "I":"GCNT2",
        "GLOB":"B3GALNT1", "GIL":"AQP3", "FORS":"GBGT1", "JR":"ABCG2",
        "LAN":"ABCB6", "VEL":"SMIM1", "AUG":"SLC29A1",
    }
    # propagate from old's system-key to corresponding gene name
    for sym, gene in SYM_TO_GENE.items():
        if sym in old_strand and gene not in old_strand:
            old_strand[gene] = old_strand[sym]

    combined = {**EXTRA_STRAND, **old_strand}  # old_strand wins
    print(f"strand map: {len(combined)} genes", file=sys.stderr)

    rows = vanno.read_text().splitlines()
    header = rows[0].split("\t")
    sys_col = header.index("system/gene")
    strand_hg19_col = header.index("strand (hg19)")
    strand_hg38_col = header.index("strand (hg38)")

    out_lines = [rows[0]]
    changed = 0
    missing: set[str] = set()
    for line in rows[1:]:
        p = line.split("\t")
        if len(p) < len(header):
            p += [""] * (len(header) - len(p))
        gene_name = p[sys_col]
        strand = combined.get(gene_name)
        if strand:
            if p[strand_hg19_col] != strand:
                changed += 1
            p[strand_hg19_col] = strand
            p[strand_hg38_col] = strand
        else:
            missing.add(gene_name)
        out_lines.append("\t".join(p))

    print(f"updated: {changed} rows", file=sys.stderr)
    if missing:
        print(f"WARN: no strand for {len(missing)} gene(s): {sorted(missing)}", file=sys.stderr)
    if args.dry_run:
        print("(dry-run, no write)", file=sys.stderr)
        return
    vanno.write_text("\n".join(out_lines) + "\n")
    print(f"OK: wrote {vanno}", file=sys.stderr)


if __name__ == "__main__":
    main()
