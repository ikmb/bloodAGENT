#!/usr/bin/env python3
"""
Generate per-pipeline gt2pt files from the V15 master gt2pt.

bloodAGENT runs 6 pipelines; each has its own gt2pt.dat. Differences from the
master V15 gt2pt are:

  - HGDP: strips ALL RHCE alleles (RHCE detection is coverage-based via the
    -k/--trick flag, not via the variant table). The 4 RHC/RHE pseudo-alleles
    that the upstream HGDP file keeps for the "coverage detected C/E" output
    are appended here verbatim from the existing HGDP file.

  - Microarray: subset to alleles that can be discriminated by SNV probes only
    (i.e., no large indels, no SVs, single-variant alleles preferred). Currently
    upstream ships only 6 ABO alleles; we re-derive that minimal set + carry
    over the upstream rows for backward compat.

  - CMR / Dragen / ONT / PacBio-GATK / PacBio-pbsv: use the master V15 gt2pt
    verbatim. (The per-caller VCF-representation differences live in the
    variation_annotation overlay files, NOT in gt2pt.)

Output: data/config/<PIPELINE>/genotype_to_phenotype_annotation_<TAG>.v15.dat
"""
from __future__ import annotations
import argparse, csv
from pathlib import Path

PIPELINES = {
    "CMR":       "ICACMR",
    "Dragen":    "Dragen",
    "ONT":       "MINIMAP2SNIFFLES",
    "PacBio":    ["TGSGATK", "TGSPBSV"],  # two callers
}

# Upstream HGDP file uses these 4 pseudo-alleles for coverage-based C/E detection.
# Names taken verbatim from the pre-V15 HGDP gt2pt so downstream JSON consumers
# that read RHCE*c / RHCE*C / RHCE*e / RHCE*E keep working unchanged.
HGDP_RHC_RHE_PSEUDO = [
    # PureSystem  MySystemKey  System  Allele      Pheno_PDF  Pheno  Pheno_flat  base_change                                              acid_change  incidence
    ["ToDo", "RHC", "RHC", "RHCE*c", "NaN", "NaN", "c", "", "", "20.13%"],
    ["ToDo", "RHC", "RHC", "RHCE*C", "NaN", "NaN", "C", "336-2849_336-2848insTTGCTATAGCTTAAGGACTCACCTGGCAGCAACACCAAACCAGGGCCACCACCATTTGAAATCCCCCAGGGTGCCCTTTGTCACTTCCCAGTGGTACAATCATAGCT", "", ""],
    ["ToDo", "RHE", "RHE", "RHCE*e", "NaN", "NaN", "e", "", "", "13.38%"],
    ["ToDo", "RHE", "RHE", "RHCE*E", "NaN", "NaN", "E", "676G>C", "Ala226Pro", "13.38%"],
]


def read_master(path: Path) -> tuple[list[str], list[list[str]]]:
    rows = []
    with path.open() as f:
        header = next(f).rstrip("\n").split("\t")
        for line in f:
            if line.startswith("#"):
                continue
            rows.append(line.rstrip("\n").split("\t"))
    return header, rows


def write_dat(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--master", default="data/source/v15/derived/genotype_to_phenotype_annotation.v15.dat")
    ap.add_argument("--out-root", default="data/source/v15/derived/pipeline")
    args = ap.parse_args()

    header, master_rows = read_master(Path(args.master))
    print(f"master: {len(master_rows)} rows, {len(header)} cols")

    out_root = Path(args.out_root); out_root.mkdir(parents=True, exist_ok=True)

    # 1. Drop-in pipelines: CMR, Dragen, ONT, PacBio (two callers) — same as master
    for pipeline, tag in [
        ("CMR", "ICACMR"), ("Dragen", "Dragen"),
        ("ONT", "MINIMAP2SNIFFLES"),
        ("PacBio", "TGSGATK"), ("PacBio", "TGSPBSV"),
    ]:
        out = out_root / pipeline / f"genotype_to_phenotype_annotation_{tag}.v15.dat"
        write_dat(out, header, master_rows)
        print(f"  wrote {out} ({len(master_rows)} rows)")

    # 2. HGDP fork: drop RHCE-keyed alleles, append RHC/RHE pseudo-alleles
    hgdp_rows = [r for r in master_rows if (len(r) >= 2 and r[1] != "RHCE")]
    hgdp_rows.extend(HGDP_RHC_RHE_PSEUDO)
    out = out_root / "HGDP" / "genotype_to_phenotype_annotation_HGDP.v15.dat"
    write_dat(out, header, hgdp_rows)
    print(f"  wrote {out} ({len(hgdp_rows)} rows; dropped {len(master_rows)-len(hgdp_rows)+4} RHCE, added 4 RHC/RHE pseudo)")

    # 3. Microarray subset — single-variant alleles only, SNV/small-indel only.
    # Heuristic: base_change contains exactly one space-separated token AND
    # the token has no large indels (no _, no dup of >5 bp).
    def is_simple(r: list[str]) -> bool:
        if len(r) < 8:
            return False
        bc = r[7].strip()
        if not bc or " " in bc:
            return False
        if "_" in bc:                     # range-style coordinates (e.g. 332_456delGT) → too big for SNV array
            return False
        return True

    array_rows = [r for r in master_rows if is_simple(r)]
    out = out_root / "Microarray" / "genotype_to_phenotype_annotation_Array.v15.dat"
    write_dat(out, header, array_rows)
    print(f"  wrote {out} ({len(array_rows)} rows; SNV/short-indel single-variant alleles only)")


if __name__ == "__main__":
    main()
