#!/usr/bin/env python3
"""
Diff the V15-derived gt2pt against the current bloodAGENT gt2pt to surface
added / renamed / retired alleles and per-system delta counts.

Output (in --out-dir):
  - allele_diff.tsv: allele, status (added/retired/unchanged/renamed?), v15_system, current_system, v15_phen, current_phen
  - system_summary.tsv: system, current_count, v15_count, added, retired
"""
from __future__ import annotations
import argparse, csv
from pathlib import Path


def read_gt2pt(path: Path) -> dict[str, dict]:
    """Returns {allele: {system, gene, phenotype}}."""
    rows = {}
    with path.open() as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < len(header) or p[0].startswith("#"):
                continue
            try:
                allele = p[idx.get("Allele", 3)]
                sys_ = p[idx.get("MySystemKey", 1)]
                gene = p[idx.get("System", 2)]
                phen = p[idx.get("Phenotype", 5)] if "Phenotype" in idx else ""
            except Exception:
                continue
            rows[allele] = {"system": sys_, "gene": gene, "phenotype": phen}
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--current", required=True, help="Path to current gt2pt .dat (e.g. data/config/HGDP/genotype_to_phenotype_annotation_HGDP.dat)")
    ap.add_argument("--v15", required=True, help="Path to V15-derived .dat (e.g. derived/genotype_to_phenotype_annotation.v15.dat)")
    ap.add_argument("--out-dir", default="derived")
    args = ap.parse_args()

    cur = read_gt2pt(Path(args.current))
    new = read_gt2pt(Path(args.v15))

    cur_alleles = set(cur); new_alleles = set(new)
    added = new_alleles - cur_alleles
    retired = cur_alleles - new_alleles
    common = cur_alleles & new_alleles
    phen_changed = [a for a in common if cur[a]["phenotype"] != new[a]["phenotype"]]

    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
    with (out / "allele_diff.tsv").open("w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["allele", "status", "v15_system", "current_system", "v15_phenotype", "current_phenotype"])
        for a in sorted(added):
            r = new[a]
            w.writerow([a, "added", r["system"], "", r["phenotype"], ""])
        for a in sorted(retired):
            r = cur[a]
            w.writerow([a, "retired", "", r["system"], "", r["phenotype"]])
        for a in sorted(phen_changed):
            w.writerow([a, "phen_changed", new[a]["system"], cur[a]["system"], new[a]["phenotype"], cur[a]["phenotype"]])

    # Per-system summary
    sys_cur = {}; sys_new = {}
    for r in cur.values():
        sys_cur[r["system"]] = sys_cur.get(r["system"], 0) + 1
    for r in new.values():
        sys_new[r["system"]] = sys_new.get(r["system"], 0) + 1
    all_sys = sorted(set(sys_cur) | set(sys_new))

    with (out / "system_summary.tsv").open("w") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["system", "current_count", "v15_count", "delta"])
        for s in all_sys:
            c, n = sys_cur.get(s, 0), sys_new.get(s, 0)
            w.writerow([s, c, n, n - c])

    print(f"added: {len(added)}")
    print(f"retired: {len(retired)}")
    print(f"phenotype_changed (same allele name): {len(phen_changed)}")
    print(f"reports written under {out}/")


if __name__ == "__main__":
    main()
