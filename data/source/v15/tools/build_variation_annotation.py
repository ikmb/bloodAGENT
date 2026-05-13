#!/usr/bin/env python3
"""
Build bloodAGENT-format variation_annotation.dat from ISBT V15 API dumps.

Inputs (in --raw-dir, default raw/):
  - allele.json: bulk allele list (id, system, gene, isbt_allele, isbt_snp, ...)
  - alleles/<id>.json: per-allele detail with embedded variants[] (call fetch_alleles.sh first)
  - variant.json: bulk variant list (hgvs_*, grch37_*, grch38_*)
  - gene.json: gene metadata (chromosome, hg19/hg38 start/end, transcript)
  - system.json: system metadata (id -> symbol)

Output (in --out-dir, default derived/):
  - variation_annotation.v15.dat: tab-separated, schema = bloodAGENT master schema (28 cols)
  - allele_to_variants.tsv: allele_id, isbt_allele, system, gene, variant_id_csv (for gt2pt build)
  - missing_coords_report.tsv: variants without complete hg19/hg38 coords (need liftOver)
  - unparseable_alleles.tsv: alleles with no variants[] or only freeform changes

bloodAGENT schema columns (TAB):
  1  system/gene                          -- gene.name (e.g. ABO, RHD, SLC44A2)
  2  Transcript annotation                -- hgvs_transcript (NM_*:c.*)
  3  Transcript annotation short          -- c.*-only suffix (after ':')
  4  genomic annotation                   -- hgvs_genomic_grch37 (NC_*:g.*) -- legacy use
  5  Chr (hg19)                           -- grch37_chr without 'chr'
  6  Pos (hg19)                           -- grch37_pos (1-based, VCF style)
  7  Chrom (hg19)                         -- 'chr' + grch37_chr
  8  0-based start (hg19)                 -- grch37_pos - 1
  9  1-based end (hg19)                   -- grch37_pos + len(ref) - 1
 10  strand (hg19)                        -- gene.strand (looked up from gene.json)
 11  Coordinate in browser format +1bp flanks (hg19)
                                          -- "chrN:start-end"
 12  Reference base (hg19)                -- transcript-strand base; for - strand revcomp(ref)
 13  Reference base blood group annotation -- same as col 12 (PDF heritage)
 14  is transcript_NC == hg19_NC          -- "TRUE" if hgvs_genomic_grch37 NC matches expected
 15  Chrom (hg38)                         -- 'chr' + grch38_chr
 16  0-based start (hg38)                 -- grch38_pos - 1
 17  1-based end (hg38)                   -- grch38_pos + len(ref) - 1
 18  strand (hg38)                        -- gene.strand
 19  Reference base (hg38)                -- transcript-strand base; for - strand revcomp
 20  is transcript_NC == hg38_NC          -- "TRUE" if hgvs_genomic_grch38 NC matches
 21  reference allele of builds concordant? -- "TRUE" if hg19_ref == hg38_ref
 22  Coordinate in VCF hg19               -- grch37_pos
 23  RefAllele in VCF hg19                -- grch37_ref
 24  AltAllele in VCF hg19                -- grch37_alt
 25  Coordinate in VCF hg38               -- grch38_pos
 26  RefAllele in VCF hg38                -- grch38_ref
 27  AltAllele in VCF hg38                -- grch38_alt
 28  TYPE                                 -- "SNV" | "ins" | "del" | "delins"

NOTE: This builds the *master* file. Pipeline-specific overlay files (e.g. for the
RHCE 109bp insertion variant) are NOT generated here — they remain hand-curated.
"""
from __future__ import annotations
import argparse, json, os, re, sys
from pathlib import Path

COLUMNS = [
    "system/gene", "Transcript annotation", "Transcript annotation short",
    "genomic annotation", "Chr (hg19)", "Pos (hg19)", "Chrom (hg19)",
    "0-based start (hg19)", "1-based end (hg19)", "strand (hg19)",
    "Coordinate in browser format +1bp flanks  (hg19)",  # note: two spaces, matches existing file
    "Reference base (hg19)", "Reference base blood group annotation",
    "is transcript_NC == hg19_NC",
    "Chrom (hg38)", "0-based start (hg38)", "1-based end (hg38)", "strand (hg38)",
    "Reference base (hg38)", "is transcript_NC == hg38_NC",
    "reference allele of builds concordant?",
    "Coordinate in VCF hg19", "RefAllele in VCF hg19", "AltAllele in VCF hg19",
    "Coordinate in VCF hg38", "RefAllele in VCF hg38", "AltAllele in VCF hg38",
    "TYPE",
]
assert len(COLUMNS) == 28

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

def revcomp(s: str) -> str:
    return s.translate(COMPLEMENT)[::-1]


def sanitize(s) -> str:
    """TSV-safe: collapse all whitespace runs into single space."""
    if s is None:
        return ""
    return re.sub(r"\s+", " ", str(s)).strip()


def classify_type(ref: str, alt: str) -> str:
    if not ref or not alt:
        return "?"
    if len(ref) == 1 and len(alt) == 1:
        return "SNV"
    if len(ref) < len(alt) and alt.startswith(ref):
        return "ins"
    if len(ref) > len(alt) and ref.startswith(alt):
        return "del"
    return "delins"


def short_hgvs(hgvs_transcript: str) -> str:
    """'NM_020469.2:c.1061delC' -> '1061delC' (strip transcript prefix + 'c.' prefix).
    Also handles already-stripped inputs like 'c.182-184_323-389del'."""
    if not hgvs_transcript:
        return ""
    suffix = hgvs_transcript.split(":", 1)[-1]  # works even if no ':' present
    # Strip leading 'c.' / 'n.' / 'g.' etc.
    return re.sub(r"^[a-z]\.", "", suffix)


# Curator-marked "high impact" variants from the upstream bloodAGENT data.
# These are the canonical antigen-defining SNVs (ABO A/B/O, JK Jk(a)/Jk(b),
# P1PK P1/P2, RHCE C/c+E/e). bloodAGENT's CISBTAnno strips the '!' at index
# time but downstream curators rely on it as a marker. Re-add when regenerating.
HIGH_IMPACT_MARKERS = {
    # keys are (system_key, short_hgvs) — system_key matches what system_key_for() returns
    # ('JK' not 'SLC14A1', 'P1PK' not 'A4GALT' — ISBT system symbol; RHCE/RHD split out from RH).
    ("ABO", "803G>C"), ("ABO", "802G>A"), ("ABO", "467C>T"),
    ("ABO", "261del"),  # V15 short HGVS drops the deleted base (was 261delG)
    ("JK", "499A>G"),
    ("P1PK", "109A>G"),
    ("RHCE", "676G>C"), ("RHCE", "122A>G"),
}


def extract_NC(hgvs_genomic: str | None) -> str | None:
    if not hgvs_genomic:
        return None
    m = re.match(r"(NC_\d+\.\d+):", hgvs_genomic)
    return m.group(1) if m else None


def system_key_for(gene_name: str, system_symbol: str) -> str:
    """Compute the bloodAGENT `system/gene` column value, matching the
    upstream convention: most systems use the ISBT symbol, but RH / MNS
    are split by gene name and FUT2 is its own pseudo-system. Must match
    the MySystemKey logic in build_gt2pt.py exactly (joined at runtime
    via C++ m_allele_vector[MySystemKey] = m_isbt_variant_to_index[system/gene])."""
    if system_symbol == "RH":           # RHD vs RHCE
        return gene_name
    if system_symbol == "MNS":          # GYPA vs GYPB vs GYPE vs MNS_HYBRID
        return gene_name
    if system_symbol == "CH_RG":        # C4A vs C4B (paralog collision)
        return gene_name
    if system_symbol == "H" and gene_name == "FUT2":
        return "FUT2"
    return system_symbol or gene_name


def build_row(v: dict, gene_record: dict | None, system_symbol: str, strand_map: dict | None = None) -> list[str]:
    gene_name = v.get("gene_name") or (gene_record or {}).get("name", "")
    sys_key = system_key_for(gene_name, system_symbol)
    transcript_hgvs = v.get("hgvs_transcript") or ""
    g37 = v.get("hgvs_genomic_grch37") or ""
    g38 = v.get("hgvs_genomic_grch38") or ""
    chr37 = (v.get("grch37_chr") or "").lstrip("chr")
    chr38 = (v.get("grch38_chr") or "").lstrip("chr")
    pos37 = v.get("grch37_pos")
    pos38 = v.get("grch38_pos")
    ref37, alt37 = v.get("grch37_ref") or "", v.get("grch37_alt") or ""
    ref38, alt38 = v.get("grch38_ref") or "", v.get("grch38_alt") or ""

    # Strand: prefer external strand_map (gene_name -> '+'/'-'), then gene_record, then default '+'
    strand = (strand_map or {}).get(gene_name) or (gene_record or {}).get("strand") or "+"

    vcf_type = classify_type(ref37 or ref38, alt37 or alt38)

    def base_at(ref: str) -> str:
        if not ref:
            return ""
        # For indels, the bloodAGENT col 12 is typically the inserted/deleted base only
        if vcf_type == "ins":
            return alt37[len(ref37):] if ref37 and alt37 else ""
        if vcf_type == "del":
            return ref37[len(alt37):] if ref37 and alt37 else ""
        return ref

    # Reference base columns mirror the VCF orientation (NOT transcript strand)
    # because CISBTAnno::isVcfAlleleAnIsbtVariant compares them directly to VCF allele
    # strings, which are always +-strand per the VCF spec. Verified against existing
    # master: ABO (- strand) row has VCF_ref='C' and Reference base='C'.
    ref_base_hg19 = base_at(ref37)
    ref_base_hg38 = base_at(ref38)

    expected_NC_hg19 = None  # TODO: derive from gene.json hg19 chrom -> NC mapping table
    nc_hg19 = extract_NC(g37)
    is_concordant_hg19 = "TRUE" if not expected_NC_hg19 or nc_hg19 == expected_NC_hg19 else "FALSE"

    expected_NC_hg38 = None
    nc_hg38 = extract_NC(g38)
    is_concordant_hg38 = "TRUE" if not expected_NC_hg38 or nc_hg38 == expected_NC_hg38 else "FALSE"

    ref_concordant = "TRUE" if (ref37 and ref38 and ref37 == ref38) else ("FALSE" if (ref37 and ref38) else "")

    start37_0 = (pos37 - 1) if isinstance(pos37, int) else ""
    end37_1 = (pos37 + len(ref37) - 1) if isinstance(pos37, int) and ref37 else ""
    start38_0 = (pos38 - 1) if isinstance(pos38, int) else ""
    end38_1 = (pos38 + len(ref38) - 1) if isinstance(pos38, int) and ref38 else ""

    browser_hg19 = f"chr{chr37}:{start37_0}-{end37_1}" if chr37 and start37_0 != "" else ""

    sh = short_hgvs(transcript_hgvs)
    # Re-prepend '!' for hand-curated high-impact variants.
    # HIGH_IMPACT_MARKERS keys are (system_key, sh) — system_key is what
    # the existing master uses (ABO/SLC14A1/A4GALT/RHCE), matching system_key_for().
    if (sys_key, sh) in HIGH_IMPACT_MARKERS:
        sh = "!" + sh

    return [
        sys_key,
        transcript_hgvs,
        sh,
        g37,
        chr37, str(pos37 or ""), f"chr{chr37}" if chr37 else "",
        str(start37_0), str(end37_1), strand,
        browser_hg19,
        ref_base_hg19, ref_base_hg19, is_concordant_hg19,
        f"chr{chr38}" if chr38 else "",
        str(start38_0), str(end38_1), strand,
        ref_base_hg38, is_concordant_hg38, ref_concordant,
        str(pos37 or ""), ref37, alt37,
        str(pos38 or ""), ref38, alt38,
        vcf_type,
    ]


def load_strand_map(old_master_path: Path | None) -> dict:
    """Build {gene_name -> '+'/'-'} from:
       1. Existing master (data/config/variation_annotation.dat) — system/gene col + strand col
       2. SYM_TO_GENE alias (P1PK -> A4GALT, etc) — propagates old system-keyed strand to gene name
       3. EXTRA_STRAND hardcoded for V15 new systems (12 systems not in old master)
    """
    SYM_TO_GENE = {
        "P1PK":"A4GALT", "LU":"BCAM", "LE":"FUT3", "FY":"ACKR1", "JK":"SLC14A1",
        "DI":"SLC4A1", "YT":"ACHE", "SC":"ERMAP", "DO":"ART4", "CO":"AQP1",
        "LW":"ICAM4", "H":"FUT1", "GE":"GYPC", "CROM":"CD55", "KN":"CR1",
        "IN":"CD44", "OK":"BSG", "RAPH":"CD151", "JMH":"SEMA7A", "I":"GCNT2",
        "GLOB":"B3GALNT1", "GIL":"AQP3", "FORS":"GBGT1", "JR":"ABCG2",
        "LAN":"ABCB6", "VEL":"SMIM1", "AUG":"SLC29A1",
    }
    EXTRA = {
        "C4A": "+", "C4B": "+", "PRNP": "+", "B4GALNT2": "-", "SLC44A2": "+",
        "ABCC4": "-", "EMP3": "+", "PIGG": "+", "ABCC1": "+", "PIEZO1": "-",
        "CD36": "-", "ATP11C": "-", "MAL": "-", "MNS_HYBRID": "-",
        "RSRP1": "+", "CD99": "-",
    }
    out = dict(EXTRA)
    if old_master_path and old_master_path.exists():
        with old_master_path.open() as f:
            header = next(f).rstrip("\n").split("\t")
            try:
                sys_col = header.index("system/gene")
                strand_col = header.index("strand (hg19)")
            except ValueError:
                return out
            for line in f:
                if line.startswith("#"): continue
                p = line.rstrip("\n").split("\t")
                if len(p) <= max(sys_col, strand_col): continue
                key, strand = p[sys_col], p[strand_col]
                if strand in ("+", "-") and key:
                    out.setdefault(key, strand)
                    if key in SYM_TO_GENE:
                        out.setdefault(SYM_TO_GENE[key], strand)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw-dir", default="raw")
    ap.add_argument("--out-dir", default="derived")
    ap.add_argument("--allele-detail-dir", default="raw/alleles")
    ap.add_argument("--old-master", default="data/config/variation_annotation.dat",
                    help="Existing bloodAGENT master used to backport strand column.")
    args = ap.parse_args()

    raw = Path(args.raw_dir)
    out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)

    systems = {s["id"]: s["symbol"] for s in json.loads((raw / "system.json").read_text())}
    genes_by_id = {g["id"]: g for g in json.loads((raw / "gene.json").read_text())}
    genes_by_name = {g["name"]: g for g in json.loads((raw / "gene.json").read_text())}
    variants_by_id = {v["id"]: v for v in json.loads((raw / "variant.json").read_text())}
    strand_map = load_strand_map(Path(args.old_master))
    print(f"strand map: {len(strand_map)} gene(s) loaded", file=sys.stderr)

    rows = []
    missing = []
    unparseable = []
    allele_to_variants_rows = []

    allele_detail_dir = Path(args.allele_detail_dir)
    allele_files = sorted(allele_detail_dir.glob("*.json"))
    print(f"Reading {len(allele_files)} allele detail files...", file=sys.stderr)

    for fp in allele_files:
        try:
            a = json.loads(fp.read_text())
        except Exception:
            unparseable.append((fp.name, "json_parse_error"))
            continue
        if a.get("deleted") or a.get("obsolete"):
            continue
        isbt_allele = a.get("isbt_allele") or ""
        sys_sym = (a.get("system") or {}).get("symbol", "")
        gene_name = (a.get("gene") or {}).get("name", "")
        gene_record = genes_by_name.get(gene_name)
        if not isbt_allele:
            unparseable.append((fp.name, "no_isbt_allele"))
            continue
        variants = a.get("variants") or []
        if not variants:
            unparseable.append((isbt_allele, "no_variants"))
            continue
        var_ids = []
        for v in variants:
            var_ids.append(v["id"])
            # Sanity check: both hg19 + hg38 coords required
            if not (v.get("grch37_pos") and v.get("grch38_pos") and v.get("grch37_ref") and v.get("grch38_ref")):
                missing.append((isbt_allele, v["id"], v.get("hgvs_transcript") or "", v.get("hgvs_genomic_grch37") or "", v.get("hgvs_genomic_grch38") or ""))
                continue
            # Skip variants without a transcript-level HGVS — these are typically
            # large SVs (multi-kb deletions) that bloodAGENT detects via coverage,
            # not via the ISBT variant table (would crash CISBTAnno indexing).
            if not v.get("hgvs_transcript"):
                missing.append((isbt_allele, v["id"], "<no hgvs_transcript>",
                                v.get("hgvs_genomic_grch37") or "",
                                v.get("hgvs_genomic_grch38") or ""))
                continue
            rows.append(build_row(v, gene_record, sys_sym, strand_map))
        allele_to_variants_rows.append((a["id"], isbt_allele, sys_sym, gene_name, ",".join(map(str, var_ids))))

    # Deduplicate rows by (gene, transcript-short)
    seen = set(); deduped = []
    for r in rows:
        key = (r[0], r[2], r[24], r[25], r[26])  # gene + short_hgvs + hg38 pos/ref/alt
        if key in seen:
            continue
        seen.add(key); deduped.append(r)

    # Sort by system/gene then hg19 pos
    deduped.sort(key=lambda r: (r[0], int(r[5]) if r[5].isdigit() else 0))

    out_dat = out / "variation_annotation.v15.dat"
    with out_dat.open("w") as f:
        f.write("\t".join(COLUMNS) + "\n")
        for r in deduped:
            f.write("\t".join(sanitize(c) for c in r) + "\n")

    with (out / "allele_to_variants.tsv").open("w") as f:
        f.write("allele_id\tisbt_allele\tsystem\tgene\tvariant_ids\n")
        for r in allele_to_variants_rows:
            f.write("\t".join(map(str, r)) + "\n")

    with (out / "missing_coords_report.tsv").open("w") as f:
        f.write("isbt_allele\tvariant_id\thgvs_transcript\thgvs_genomic_grch37\thgvs_genomic_grch38\n")
        for r in missing:
            f.write("\t".join(map(str, r)) + "\n")

    with (out / "unparseable_alleles.tsv").open("w") as f:
        f.write("isbt_allele_or_file\treason\n")
        for r in unparseable:
            f.write("\t".join(map(str, r)) + "\n")

    print(f"OK: wrote {len(deduped)} variant rows to {out_dat}", file=sys.stderr)
    print(f"  - alleles with no variants: {len(unparseable)}", file=sys.stderr)
    print(f"  - variants missing hg19/hg38 coords: {len(missing)}", file=sys.stderr)


if __name__ == "__main__":
    main()
