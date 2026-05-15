#!/usr/bin/env python3
"""
Compare V15-run bloodAGENT phenotype JSONs against the pre-V15 baseline JSONs
bundled in data/testdata/. Emits a per-sample diff highlighting allele renames
and phenotype changes.

Usage:
  python3 data/source/v15/tools/compare_regression.py \
    --baseline-dir data/testdata \
    --v15-dir data/source/v15/derived/regression

For each test sample (HGDP00001/00003/00005, NA24143):
  - Match the V15 *.v15.json against the pre-V15 *.phased.json (or *.json)
  - For each system: compare called allele name(s) + phenotype(s)
  - Report:
      [unchanged]  same allele + same phenotype       OK
      [renamed]    different allele name, same phenotype string  -> likely V15 rename
      [phen_diff]  same allele but phenotype string changed
      [new_sys]    system present in V15 but not in baseline (the 12 new systems)
      [lost_sys]   system was in baseline but missing in V15 (should be empty)
"""
from __future__ import annotations
import argparse, json, sys
from pathlib import Path


def parse_loci(json_path: Path) -> dict[str, dict]:
    """Returns {system_key: {alleles: [(haplotype_idx, names_list)], phenotypes: list, score: float}}."""
    d = json.loads(json_path.read_text())
    out = {}
    for sys_key, info in d.get("loci", {}).items():
        calls = info.get("calls", [])
        if not calls:
            continue
        call = calls[0]  # primary call (first-ranked by score)
        alleles = []
        for hap_i, a in enumerate(call.get("alleles", []) or []):
            names = a.get("names") or ([a["name"]] if a.get("name") else [])
            alleles.append((hap_i, tuple(sorted(names))))
        phenotypes = call.get("flat_phenotypes") or call.get("phenotypes") or []
        out[sys_key] = {
            "alleles": alleles,
            "phenotypes": phenotypes,
            "score": call.get("score"),
        }
    return out


def fmt_alleles(loc: dict) -> str:
    return " | ".join("+".join(names) for _, names in loc["alleles"])


def fmt_phen(loc: dict) -> str:
    p = loc["phenotypes"]
    if not p:
        return "-"
    flat = []
    for hap in p:
        if isinstance(hap, list):
            flat.append("/".join(map(str, hap)) or "-")
        else:
            flat.append(str(hap))
    return " | ".join(flat)


def diff_sample(baseline: Path, v15: Path):
    print(f"\n=== {baseline.name} vs {v15.name} ===")
    b = parse_loci(baseline) if baseline.exists() else {}
    n = parse_loci(v15) if v15.exists() else {}
    if not b:
        print(f"  (no baseline at {baseline})")
    if not n:
        print(f"  (no V15 at {v15} — skip)")
        return
    all_sys = sorted(set(b) | set(n))
    tallies = {"unchanged": 0, "renamed": 0, "phen_diff": 0, "new_sys": 0, "lost_sys": 0}
    for s in all_sys:
        in_b = s in b
        in_n = s in n
        if in_n and not in_b:
            print(f"  [new_sys] {s}: alleles={fmt_alleles(n[s])} phen={fmt_phen(n[s])}")
            tallies["new_sys"] += 1
        elif in_b and not in_n:
            print(f"  [lost_sys] {s}: alleles={fmt_alleles(b[s])} phen={fmt_phen(b[s])}")
            tallies["lost_sys"] += 1
        else:
            same_alleles = b[s]["alleles"] == n[s]["alleles"]
            same_phen = b[s]["phenotypes"] == n[s]["phenotypes"]
            if same_alleles and same_phen:
                tallies["unchanged"] += 1
            elif same_phen:
                print(f"  [renamed] {s}: {fmt_alleles(b[s])}  ->  {fmt_alleles(n[s])}    phen={fmt_phen(n[s])}")
                tallies["renamed"] += 1
            else:
                print(f"  [phen_diff] {s}: alleles={fmt_alleles(b[s])} -> {fmt_alleles(n[s])}    phen {fmt_phen(b[s])} -> {fmt_phen(n[s])}")
                tallies["phen_diff"] += 1
    print(f"  totals: {tallies}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline-dir", default="data/testdata")
    ap.add_argument("--v15-dir", default="data/source/v15/derived/regression")
    args = ap.parse_args()

    samples = [
        ("HGDP00001", "HGDP00001.phased.json"),
        ("HGDP00003", "HGDP00003.phased.json"),
        ("HGDP00005", "HGDP00005.phased.json"),
        ("NA24143",   "NA24143.json"),  # only flat JSON in repo for this one
    ]
    for s, baseline_name in samples:
        baseline = Path(args.baseline_dir) / s / baseline_name
        v15 = Path(args.v15_dir) / f"{s}.v15.json"
        diff_sample(baseline, v15)


if __name__ == "__main__":
    main()
