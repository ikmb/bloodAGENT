#!/usr/bin/env python3
"""
Generate exonic_annotation supplement rows for V15 systems NOT in the existing
bloodAGENT exonic_annotation.hg{19,38}.BGStarget.txt.

Output columns (12, TAB):
  knownToRefSeq.value · refGene.name · refGene.chrom · refGene.strand ·
  refGene.txStart · refGene.txEnd · refGene.cdsStart · refGene.cdsEnd ·
  refGene.exonStarts · refGene.exonEnds · refGene.name2 · System

Source: V15 /api/gene gives chromosome, hg19/hg38 start/end, and exon_count
but does NOT expose per-exon coordinates. We populate everything that's known
and leave `cdsStart`, `cdsEnd`, `exonStarts`, `exonEnds` empty with a comment
in the supplementary `exonic_annotation_TODO.tsv` listing what to fetch from
UCSC refGene before the file is production-ready.

For pure SNV typing the empty exon arrays do not break bloodAGENT (SNVs are
located via variation_annotation.dat coordinates, not the exon arrays). The
coverage-based RHD detection and similar features do need exon arrays — but
none of the 12 new systems are coverage-detected by upstream yet, so this is
acceptable as the initial drop-in.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path

# System -> primary gene name (for new V15 systems)
NEW_SYSTEM_GENES = {
    "CH_RG":   ["C4A", "C4B"],
    "KANNO":   ["PRNP"],
    "SID":     ["B4GALNT2"],
    "CTL2":    ["SLC44A2"],
    "PEL":     ["ABCC4"],
    "MAM":     ["EMP3"],
    "EMM":     ["PIGG"],
    "ABCC1":   ["ABCC1"],
    "ER":      ["PIEZO1"],
    "CD36":    ["CD36"],
    "ATP11C":  ["ATP11C"],
    "MAL":     ["MAL"],
}

# Hand-curated strand & RefSeq accession for these 12 systems (matches build_variation_annotation.py)
GENE_META = {
    "C4A":      {"strand":"+", "refseq":"NM_007293", "chrom":"chr6"},
    "C4B":      {"strand":"+", "refseq":"NM_001002029", "chrom":"chr6"},
    "PRNP":     {"strand":"+", "refseq":"NM_000311", "chrom":"chr20"},
    "B4GALNT2": {"strand":"-", "refseq":"NM_153446", "chrom":"chr17"},
    "SLC44A2":  {"strand":"+", "refseq":"NM_020428", "chrom":"chr19"},
    "ABCC4":    {"strand":"-", "refseq":"NM_005845", "chrom":"chr13"},
    "EMP3":     {"strand":"+", "refseq":"NM_001425", "chrom":"chr19"},
    "PIGG":     {"strand":"+", "refseq":"NM_001127178", "chrom":"chr4"},
    "ABCC1":    {"strand":"+", "refseq":"NM_004996", "chrom":"chr16"},
    "PIEZO1":   {"strand":"-", "refseq":"NM_001142864", "chrom":"chr16"},
    "CD36":     {"strand":"-", "refseq":"NM_001001548", "chrom":"chr7"},
    "ATP11C":   {"strand":"-", "refseq":"NM_173694", "chrom":"chrX"},
    "MAL":      {"strand":"-", "refseq":"NM_002371", "chrom":"chr2"},
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", default="data/source/v15/raw")
    ap.add_argument("--out-dir", default="data/source/v15/derived")
    args = ap.parse_args()

    genes = {g["name"]: g for g in json.loads(Path(args.raw_dir, "gene.json").read_text())}

    todo_rows = []  # for the gap report

    for build in ("hg19", "hg38"):
        out = Path(args.out_dir, f"exonic_annotation.{build}.BGStarget.supplement.txt")
        with out.open("w") as f:
            # NO header (this is a SUPPLEMENT meant to be appended to the existing file)
            for sys_sym, gene_list in NEW_SYSTEM_GENES.items():
                for gn in gene_list:
                    meta = GENE_META.get(gn, {})
                    g = genes.get(gn, {})
                    chrom = meta.get("chrom") or (f"chr{g.get('chromosome', '?')}" if g.get('chromosome') else "?")
                    strand = meta.get("strand") or "+"
                    refseq = meta.get("refseq") or "NM_TODO"
                    start_key, end_key = f"{build}_start", f"{build}_end"
                    tx_start = g.get(start_key) or 0
                    tx_end = g.get(end_key) or 0
                    # Per-exon coords unknown from V15 API. Leave empty; populate from UCSC refGene.
                    cds_start = ""
                    cds_end = ""
                    exon_starts = ""
                    exon_ends = ""
                    row = [
                        refseq, refseq, chrom, strand,
                        str(tx_start), str(tx_end),
                        cds_start, cds_end, exon_starts, exon_ends,
                        gn, sys_sym,
                    ]
                    f.write("\t".join(row) + "\n")
                    if build == "hg38":
                        todo_rows.append((sys_sym, gn, refseq, chrom, strand,
                                          tx_start, tx_end,
                                          "needs cdsStart, cdsEnd, exonStarts, exonEnds from UCSC refGene"))
        print(f"wrote {out}")

    todo_out = Path(args.out_dir, "exonic_annotation_TODO.tsv")
    with todo_out.open("w") as f:
        f.write("system\tgene\trefseq\tchrom\tstrand\ttxStart_hg38\ttxEnd_hg38\tnote\n")
        for r in todo_rows:
            f.write("\t".join(map(str, r)) + "\n")
    print(f"wrote {todo_out} ({len(todo_rows)} gene rows pending UCSC backfill)")


if __name__ == "__main__":
    main()
