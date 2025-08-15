#!/usr/bin/env python3
"""
Append PS to FORMAT and add a PS value to each genotype field in a VCF.

- Converts FORMAT "GT" -> "GT:PS" (if PS not already present).
- Appends ":<PS>" to every sample's genotype (e.g., "0|1" -> "0|1:1234567").
- Adds a proper FORMAT header line for PS if missing.
- Works line-by-line; supports .vcf and .vcf.gz I/O.

Usage:
  vcf_add_ps.py -i in.vcf[.gz] -o out.vcf[.gz] --ps 1234567
  vcf_add_ps.py -i in.vcf.gz -o out.vcf.gz --ps-from pos           # use POS as PS
  vcf_add_ps.py -i - -o - --ps 1234567                              # stdin -> stdout

Options:
  --overwrite-ps    If PS already exists in FORMAT, overwrite sample PS values.
"""

import argparse
import gzip
import io
import sys
from typing import TextIO

PS_HEADER = '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set identifier">'

def open_read(path: str) -> TextIO:
    if path == "-":
        return io.TextIOWrapper(sys.stdin.buffer, encoding="utf-8")
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8")
    return open(path, "rt", encoding="utf-8", newline="")

def open_write(path: str) -> TextIO:
    if path == "-":
        return io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "wb"), encoding="utf-8")
    return open(path, "wt", encoding="utf-8", newline="")

def process(in_f: TextIO, out_f: TextIO, ps_mode: str, ps_const: str, overwrite_ps: bool):
    header_lines = []
    header_done = False
    have_ps_header = False

    for raw in in_f:
        if not header_done:
            if raw.startswith("##"):
                if raw.startswith("##FORMAT=<ID=PS"):
                    have_ps_header = True
                header_lines.append(raw)
                continue
            if raw.startswith("#CHROM"):
                # inject PS header if missing
                if not have_ps_header:
                    header_lines.append(PS_HEADER + "\n")
                # flush header
                for h in header_lines:
                    out_f.write(h)
                out_f.write(raw)
                header_done = True
                continue

        # Body line
        if not raw or raw[0] == "#":
            out_f.write(raw)
            continue

        cols = raw.rstrip("\n").split("\t")
        if len(cols) < 10:
            # No sample columns; write as-is
            out_f.write(raw)
            continue

        chrom, pos = cols[0], cols[1]
        fmt = cols[8]
        samples = cols[9:]

        tags = fmt.split(":")
        have_ps = "PS" in tags
        if not have_ps:
            tags.append("PS")
            cols[8] = ":".join(tags)

        # PS value
        ps_val = pos if ps_mode == "pos" else ps_const

        # update each sample
        for i, s in enumerate(samples):
            if s == ".":
                # fill GT missing with . and PS value
                samples[i] = ".:" + ps_val
                continue

            vals = s.split(":")
            if have_ps:
                # ensure we have slot for PS
                idx = tags.index("PS")
                if len(vals) <= idx:
                    # extend with missing fields
                    vals.extend(["."] * (idx - len(vals) + 1))
                if overwrite_ps or vals[idx] in (".", ""):
                    vals[idx] = ps_val
            else:
                # PS newly added at the end
                vals.append(ps_val)
            samples[i] = ":".join(vals)

        cols[9:] = samples
        out_f.write("\t".join(cols) + "\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input", required=True, help="Input VCF (.vcf or .vcf.gz or - for stdin)")
    ap.add_argument("-o", "--output", required=True, help="Output VCF (.vcf or .vcf.gz or - for stdout)")
    g = ap.add_mutually_exclusive_group()
    g.add_argument("--ps", default="1234567", help="Constant PS value to append (default: 1234567)")
    g.add_argument("--ps-from", choices=["pos"], help="Use POS as PS (ps=position)")
    ap.add_argument("--overwrite-ps", action="store_true", help="Overwrite existing PS values if present")
    args = ap.parse_args()

    ps_mode = "pos" if args.ps_from == "pos" else "const"
    ps_const = args.ps if ps_mode == "const" else None
    if ps_mode == "const" and (ps_const is None or ps_const == ""):
        ap.error("--ps must be a non-empty value")

    with open_read(args.input) as inf, open_write(args.output) as outf:
        process(inf, outf, ps_mode, ps_const or "", args.overwrite_ps)

if __name__ == "__main__":
    main()

