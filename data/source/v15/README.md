# bloodAGENT × ISBT V15 — data integration workspace

This directory holds everything needed to align bloodAGENT's `.dat` knowledge base
with the ISBT Blood Group Database **V15** snapshot (released **2026-05-01**).

```
data/source/v15/
├── README.md                   # this file
├── SYSTEM_MAPPING.md           # ISBT 001-048 ↔ bloodAGENT system/gene mapping
├── raw/                        # raw API dumps (gitignored when large)
│   ├── system.json             # 55 entries (48 BGS + 2 TF + 3 collection + 2 series)
│   ├── gene.json               # 57 genes incl. hg19/hg38 coords
│   ├── antigen.json            # 398 antigens
│   ├── allele.json             # 2053 alleles (bulk — NO variants[])
│   ├── variant.json            # 1865 variants incl. hg19/hg38 VCF coords
│   ├── phenotype.json          # 18664 phenotype rows
│   ├── publication.json        # 854 citations
│   ├── genbank.json            # 846 reference sequences
│   ├── release.json            # release index (v1..v15)
│   ├── release_15.json         # V15 release object + updatedAlleleIds (1448 in V15)
│   ├── allele_ids.txt          # one allele id per line, drives the per-allele dump
│   └── alleles/<id>.json       # per-allele detail WITH variants[] (one file per allele)
├── tools/
│   ├── fetch_isbt_v15.sh       # one-shot script to (re)fetch raw/ from the public API
│   ├── build_variation_annotation.py
│   ├── build_gt2pt.py
│   ├── diff_against_current.py
│   └── validate_dat.py         # mimics CISBTAnno::readAnnotation+generateIndex checks
└── derived/                    # outputs of build_* scripts (created on first run)
    ├── variation_annotation.v15.dat
    ├── genotype_to_phenotype_annotation.v15.dat
    ├── allele_to_variants.tsv
    ├── missing_coords_report.tsv
    ├── unparseable_alleles.tsv
    ├── allele_diff.tsv
    └── system_summary.tsv
```

## ISBT V15 Blood Group Database — REST API reference

We discovered the database is backed by a **public NestJS REST API**. No auth
needed. Throttled (~roughly 1 req/sec sustained).

Useful endpoints:

| Method | Path | Purpose |
| --- | --- | --- |
| GET | `/api/release` | List all releases (v1..v15) |
| GET | `/api/release/15` | V15 release object: counts + `updatedAlleleIds` |
| GET | `/api/system` | All 55 system entries |
| GET | `/api/system/<id>` | Single system + nested description |
| GET | `/api/gene` | All 57 genes incl. transcript & hg19/hg38 coords |
| GET | `/api/antigen` | All 398 antigens |
| GET | `/api/allele` | All 2053 alleles (bulk — minimal fields, NO `variants[]`) |
| GET | `/api/allele/<id>` | **Single allele with embedded `variants[]` + `phenotypes[]` + `publications[]`** |
| GET | `/api/variant` | All 1865 variants (HGVS + hg19/hg38 VCF coords) |
| GET | `/api/phenotype` | All 18664 phenotype rows |
| GET | `/api/publication` | All 854 citations |
| GET | `/api/genbank` | All 846 reference sequences |

Note: query parameters like `?include=variants` are rejected with HTTP 400
(strict whitelist validation). The bulk `/api/allele` endpoint does NOT include
the variants[] array — you must call `/api/allele/<id>` for each allele to get
the allele→variant linkage. There are 2053 alleles; expect ~1 hour to dump
politely (the per-IP throttler is aggressive).

## End-to-end pipeline

```bash
# 1) Fetch raw V15 data from the public API (takes 30-60 minutes due to throttling)
data/source/v15/tools/fetch_isbt_v15.sh

# 2) Build bloodAGENT-format master files
python3 data/source/v15/tools/build_variation_annotation.py \
  --raw-dir data/source/v15/raw \
  --allele-detail-dir data/source/v15/raw/alleles \
  --out-dir data/source/v15/derived

python3 data/source/v15/tools/build_gt2pt.py \
  --raw-dir data/source/v15/raw \
  --allele-detail-dir data/source/v15/raw/alleles \
  --out-dir data/source/v15/derived

# 3) Validate the new master files won't crash C++ at runtime
python3 data/source/v15/tools/validate_dat.py \
  data/source/v15/derived/variation_annotation.v15.dat \
  data/source/v15/derived/genotype_to_phenotype_annotation.v15.dat

# 4) Diff vs current bloodAGENT data
python3 data/source/v15/tools/diff_against_current.py \
  --current data/config/HGDP/genotype_to_phenotype_annotation_HGDP.dat \
  --v15 data/source/v15/derived/genotype_to_phenotype_annotation.v15.dat \
  --out-dir data/source/v15/derived
```

## Mapping decisions baked into build_*.py

- **`system/gene` column** in bloodAGENT keeps using **gene names** (ABO, RHD,
  RHCE, GYPA, ...) because the C++ `CISBTAnno::generateIndex()` uses this string
  as a hash key. Renaming to ISBT system symbols (RH, MNS, ...) would break the
  existing genotype-to-phenotype lookups in `CIsbtGt2Pt`.

- **`MySystemKey` column** in gt2pt is now populated with the **ISBT V15 system
  symbol** (ABO, MNS, RH, ...) so downstream consumers can group by system.

- **`Phenotype_PDF_Table` column** is kept (legacy column name) but now stores
  the V15 `isbt_phenotype` value — the PDF era is over but the column persists
  for backward compatibility with anything that reads it.

- **`incidence` column** is populated from `gnomad_all` when present, otherwise
  blank. The PDF era used clinical-frequency strings ("26.50%") which are not
  in the V15 API — gnomAD is the closest objective replacement.

- **hg19 coordinates**: All V15 variants ship with both hg19 and hg38 in the
  API. New systems with hg38-only entries are flagged in `missing_coords_report.tsv`
  and need UCSC liftOver before they can be added (bloodAGENT's `CISBTAnno`
  crashes on `stoi("")` if hg19 coord is empty).

- **`Reference base (hg19)` / `Reference base (hg38)` columns** are populated by
  taking the variant's `grch3{7,8}_ref`, applying strand-aware revcomp when
  `gene.strand == '-'`, and stripping the VCF anchor base for indels (so the
  column carries only the inserted/deleted bases on transcript strand). This
  matches the convention in the existing master `variation_annotation.dat`.

## What's still manual

1. **Pipeline-specific overlay files** for the RHCE 109bp insertion (CMR/ICA,
   Dragen, ONT/sniffles, PacBio/GATK, PacBio/pbsv). These contain the same
   variant represented differently per VCF caller; they are NOT regenerated
   from V15 and remain hand-curated.

2. **exonic_annotation.{hg19,hg38}.BGStarget.txt** — needs new system exons
   (PIEZO1, CD36, MAL, ATP11C, PIGG, B4GALNT2, SLC44A2, etc.) appended. Build
   from `gene.json` columns `hg{19,38}_{start,end}` + `exon_count`. Script
   not yet written (deferred to phase 4).

3. **`detect_RHCplusminus.py` and `count_GYP.py`** — secondary-analysis helpers
   that hard-code allele names. Verify against V15 RH and MNS allele lists for
   any rename impact.

4. **gnomAD-derived incidence** is not directly comparable to the PDF-era
   clinical-frequency numbers. For ABO/RHD reference alleles the PDF said
   ~91% / ~38% (clinical), gnomAD gives much lower numbers (per-variant AF).
   Consider keeping the old `incidence` for reference alleles and using gnomAD
   only for novel V15 alleles.

## Verifying the C++ build (optional, deferred)

The C++ project uses NetBeans-generated Makefiles and htslib's autotools build.
Native macOS build is painful; the cleanest path is the existing Dockerfile.
Since this upgrade is **data-only**, we deliberately skip the C++ rebuild during
data preparation and validate via `validate_dat.py` instead.

Once the .dat changes are merged, a regression run against the four bundled
test samples (`data/testdata/HGDP00001` etc.) inside the rebuilt Singularity
image is the authoritative check.
