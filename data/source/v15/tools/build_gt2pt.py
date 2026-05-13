#!/usr/bin/env python3
"""
Build bloodAGENT-format genotype_to_phenotype_annotation*.dat from ISBT V15 API dumps.

Schema (10 cols, TAB):
  1 PureSystem         -- "ToDo" placeholder (preserved from upstream convention)
  2 MySystemKey        -- ISBT system symbol, e.g. ABO, P1PK, RH (was previously gene name; we now align with ISBT)
  3 System             -- gene name, e.g. ABO, A4GALT, RHCE, RHD
  4 Allele             -- isbt_allele, e.g. ABO*A1.01
  5 Phenotype_PDF_Table -- isbt_phenotype (verbatim from V15)
  6 Phenotype          -- isbt_phenotype (same as col 5, kept for backward-compat downstream)
  7 Phenotype_flat     -- isbt_phenotype with HTML stripped + spaces normalised
  8 base_change        -- space-joined list of variant 'Transcript annotation short' (c.* without prefix)
  9 acid_change        -- space-joined list of hgvs_predicted_protein, simplified
 10 incidence          -- gnomad_all (max across variants) OR blank

Inputs (in --raw-dir):
  - alleles/<id>.json (per-allele detail with variants[] embedded)
  - allele.json (bulk list for cross-check)
  - system.json (id -> symbol)

Output (in --out-dir):
  - genotype_to_phenotype_annotation.v15.dat
"""
from __future__ import annotations
import argparse, json, re
from pathlib import Path

COLUMNS = [
    "PureSystem", "MySystemKey", "System", "Allele",
    "Phenotype_PDF_Table", "Phenotype", "Phenotype_flat",
    "base_change", "acid_change", "incidence",
]


def short_hgvs(hgvs_transcript: str) -> str:
    if not hgvs_transcript:
        return ""
    suffix = hgvs_transcript.split(":", 1)[-1]  # works whether or not ':' is present
    return re.sub(r"^[a-z]\.", "", suffix)


def short_protein(hgvs_protein: str) -> str:
    """'NP_002091.4:p.(Ile66_Ile78dup)' -> 'Ile66_Ile78dup'."""
    if not hgvs_protein or ":" not in hgvs_protein:
        return ""
    suffix = hgvs_protein.split(":", 1)[1]
    suffix = re.sub(r"^p\.", "", suffix)
    return suffix.strip("()")


def strip_html(s: str) -> str:
    if not s:
        return ""
    s = re.sub(r"<[^>]+>", "", s)
    # Collapse all whitespace runs (incl. newlines/tabs) to single space — TSV-safe
    return re.sub(r"\s+", " ", s).strip()


def sanitize(s: str) -> str:
    """TSV-safe: collapse all whitespace incl newlines/tabs into single space."""
    if not s:
        return ""
    return re.sub(r"\s+", " ", str(s)).strip()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", default="raw")
    ap.add_argument("--allele-detail-dir", default="raw/alleles")
    ap.add_argument("--out-dir", default="derived")
    args = ap.parse_args()

    raw = Path(args.raw_dir); out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    allele_detail_dir = Path(args.allele_detail_dir)
    files = sorted(allele_detail_dir.glob("*.json"))

    rows = []
    for fp in files:
        try:
            a = json.loads(fp.read_text())
        except Exception:
            continue
        if a.get("deleted") or a.get("obsolete"):
            continue
        isbt_allele = a.get("isbt_allele")
        if not isbt_allele:
            continue
        sys_sym = (a.get("system") or {}).get("symbol", "")
        gene_name = (a.get("gene") or {}).get("name", "")
        phen = strip_html(a.get("isbt_phenotype") or "")
        variants = a.get("variants") or []
        base_changes = []
        acid_changes = []
        gnomad_max = 0.0
        for v in variants:
            sh = short_hgvs(v.get("hgvs_transcript") or "")
            if not sh and v.get("freeform_dna_change"):
                sh = v["freeform_dna_change"]
            if sh:
                base_changes.append(sh)
            sp = short_protein(v.get("hgvs_predicted_protein") or "")
            if not sp and v.get("freeform_predicted_protein"):
                sp = v["freeform_predicted_protein"]
            if sp:
                acid_changes.append(sp)
            try:
                af = v.get("gnomad_all")
                if af is not None:
                    af = float(af)
                    if af > gnomad_max:
                        gnomad_max = af
            except Exception:
                pass

        incidence = ""
        if gnomad_max > 0:
            incidence = f"{gnomad_max*100:.2f}%"

        # MySystemKey alignment — MUST match build_variation_annotation.py:system_key_for()
        # exactly, because the C++ joins vanno[system/gene] <-> gt2pt[MySystemKey].
        if sys_sym in ("RH", "MNS", "CH_RG"):
            my_system_key = gene_name           # RHD/RHCE, GYPA/GYPB/..., C4A/C4B
        elif sys_sym == "H" and gene_name == "FUT2":
            my_system_key = "FUT2"              # Secretor split off as its own pseudo-system
        else:
            my_system_key = sys_sym or gene_name
        rows.append([
            "ToDo",
            sanitize(my_system_key),
            sanitize(gene_name),
            sanitize(isbt_allele),
            sanitize(phen),
            sanitize(phen),
            sanitize(phen),
            sanitize(" ".join(base_changes)),
            sanitize(" ".join(acid_changes)),
            sanitize(incidence),
        ])

    # Sort: system, allele
    rows.sort(key=lambda r: (r[1], r[3]))

    out_path = out / "genotype_to_phenotype_annotation.v15.dat"
    with out_path.open("w") as f:
        f.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    print(f"OK: wrote {len(rows)} rows to {out_path}")


if __name__ == "__main__":
    main()
