#!/usr/bin/env python3
"""
Regression analysis V15 vs pre-V15:
  - For each baseline-called allele, check that it still appears in V15's call set.
  - For each baseline phenotype, check that it still appears in V15's phenotype set.
  - Quantify call-set expansion (more alleles tying at top score is a known
    consequence of V15's larger allele table when input VCF can't discriminate).
  - List the new V15 systems that produced calls.

Usage:
  python3 data/source/v15/tools/analyze_regression.py
"""
from __future__ import annotations
import argparse, json
from collections import defaultdict
from pathlib import Path


def parse_loci(p: Path) -> dict:
    if not p.exists():
        return {}
    d = json.loads(p.read_text())
    out = {}
    for sys_key, info in d.get("loci", {}).items():
        calls = info.get("calls", [])
        if not calls:
            continue
        call = calls[0]
        # Collect all allele names across the 2 haplotypes (set of names per haplotype)
        haps_alleles = []
        haps_phens = []
        for a in call.get("alleles", []) or []:
            names = set(a.get("names") or ([] if a.get("name") is None else [a["name"]]))
            haps_alleles.append(names)
        for ph in call.get("flat_phenotypes") or call.get("phenotypes") or []:
            haps_phens.append(set(ph) if isinstance(ph, list) else {ph})
        out[sys_key] = {
            "alleles_per_hap": haps_alleles,
            "phens_per_hap": haps_phens,
            "score": call.get("score"),
        }
    return out


def hap_match(baseline_set: set, v15_set: set) -> bool:
    """True if every baseline call is preserved in V15 (V15 may have MORE alleles)."""
    return bool(baseline_set) and baseline_set.issubset(v15_set)


def analyze(sample: str, baseline_path: Path, v15_path: Path) -> dict:
    base = parse_loci(baseline_path)
    v15 = parse_loci(v15_path)
    if not base or not v15:
        return {"sample": sample, "skip": True, "base_path": str(baseline_path), "v15_path": str(v15_path)}

    common = set(base) & set(v15)
    only_v15 = set(v15) - set(base)
    only_base = set(base) - set(v15)
    allele_preserved = []   # systems where every baseline allele is still in V15's top set
    allele_dropped = []     # systems where V15 LOST a baseline allele
    phen_preserved = []
    phen_dropped = []
    expansion = defaultdict(int)  # how many extra alleles per haplotype V15 introduces

    for s in common:
        # diploid: 2 haplotypes; require BOTH haps to preserve baseline allele
        hap_count = max(len(base[s]["alleles_per_hap"]), len(v15[s]["alleles_per_hap"]))
        if hap_count == 0:
            continue
        a_ok = all(
            hap_match(base[s]["alleles_per_hap"][i], v15[s]["alleles_per_hap"][i])
            for i in range(min(len(base[s]["alleles_per_hap"]), len(v15[s]["alleles_per_hap"])))
        )
        if a_ok:
            allele_preserved.append(s)
        else:
            allele_dropped.append(s)
        p_ok = all(
            hap_match(base[s]["phens_per_hap"][i], v15[s]["phens_per_hap"][i])
            for i in range(min(len(base[s]["phens_per_hap"]), len(v15[s]["phens_per_hap"])))
        )
        if p_ok:
            phen_preserved.append(s)
        else:
            phen_dropped.append(s)
        for i in range(min(len(base[s]["alleles_per_hap"]), len(v15[s]["alleles_per_hap"]))):
            expansion[s] = max(expansion[s], len(v15[s]["alleles_per_hap"][i]) - len(base[s]["alleles_per_hap"][i]))

    return {
        "sample": sample,
        "common_systems": len(common),
        "only_v15_systems": sorted(only_v15),
        "only_base_systems": sorted(only_base),
        "allele_preserved": len(allele_preserved),
        "allele_dropped_systems": sorted(allele_dropped),
        "phen_preserved": len(phen_preserved),
        "phen_dropped_systems": sorted(phen_dropped),
        "expansion_top": sorted(expansion.items(), key=lambda x: -x[1])[:10],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline-dir", default="data/testdata")
    ap.add_argument("--v15-dir", default="data/source/v15/derived/regression")
    args = ap.parse_args()

    samples = [
        ("HGDP00001", "HGDP00001.phased.json"),
        ("HGDP00003", "HGDP00003.phased.json"),
        ("HGDP00005", "HGDP00005.phased.json"),
        ("NA24143",   "NA24143.json"),
    ]
    for s, baseline_name in samples:
        b = Path(args.baseline_dir) / s / baseline_name
        v = Path(args.v15_dir) / f"{s}.v15.json"
        r = analyze(s, b, v)
        if r.get("skip"):
            print(f"=== {s} === SKIP (missing baseline or V15)")
            continue
        print(f"=== {s} ===")
        print(f"  common systems: {r['common_systems']}")
        print(f"  systems gained in V15: {r['only_v15_systems'] or '-'}")
        print(f"  systems lost in V15:   {r['only_base_systems'] or '-'}")
        print(f"  baseline allele PRESERVED in V15 top set: {r['allele_preserved']}/{r['common_systems']}")
        print(f"  baseline allele DROPPED:                  {r['allele_dropped_systems'] or '-'}")
        print(f"  baseline phenotype PRESERVED:             {r['phen_preserved']}/{r['common_systems']}")
        print(f"  baseline phenotype DROPPED:               {r['phen_dropped_systems'] or '-'}")
        if r["expansion_top"]:
            print(f"  call-set expansion (more alleles tied at top score) — top 5:")
            for sys_k, extra in r["expansion_top"][:5]:
                print(f"    {sys_k}: +{extra} extra alleles per haplotype")


if __name__ == "__main__":
    main()
