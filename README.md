# bloodAGENT (ISBT V15 fork)

**A fork of [ikmb/bloodAGENT](https://github.com/ikmb/bloodAGENT)** with the reference allele tables refreshed against the new **ISBT Blood Group Database V15** (May 2026 release). The C++ engine is unchanged; only the data layer under `data/config/` and a new toolchain under `data/source/v15/` have been added.

| | |
| --- | --- |
| Upstream | https://github.com/ikmb/bloodAGENT (Wittig lab, IKMB Kiel) |
| This fork | https://github.com/gangchen/bloodAGENT |
| Reference data | **ISBT Blood Group Database V15** (applied 2026-05-01) |
| Coverage | 48 blood group systems · 57 genes · 397 antigens · 2 036 alleles · 1 828 variants |
| License | BSD 2-Clause (unchanged from upstream) |

---

## Table of Contents

- [What this fork changes](#what-this-fork-changes)
- [Quick start (Docker)](#quick-start-docker)
- [Project origin](#project-origin)
- [Reference data: ISBT V15](#reference-data-isbt-v15)
- [V15 regression results](#v15-regression-results)
- [Regenerating the data from a future ISBT release](#regenerating-the-data-from-a-future-isbt-release)
- [Introduction (from upstream)](#introduction-from-upstream)
- [Key Features](#key-features)
- [System Requirements](#system-requirements)
- [Installation](#installation)
- [Input Data Format](#input-data-format)
  - [Testdata](#testdata)
- [Cosine Similarity Scoring](#cosine-similarity-scoring)
- [Running bloodAGENT](#running-bloodagent)
- [Output Format](#output-format)
- [Custom Secondary Analysis (RHCE)](#how-to-run-custom-secondary-analysis-scripts)
- [Special Case: RHD](#special-case-rhd)
- [Limitations](#limitations)
- [Licensing](#licensing)

---

## What this fork changes

Compared to the upstream commit (`fee47fd`, Jan 2026):

1. **All allele tables refreshed to ISBT V15** (May 2026 snapshot). 12 brand-new blood group systems gain first-time bloodAGENT coverage (CH_RG, KANNO, SID, CTL2, PEL, MAM, EMM, ABCC1, ER, CD36, ATP11C, MAL); all existing systems get added/renamed alleles. See [`data/source/v15/derived/baseline_analysis.md`](data/source/v15/derived/baseline_analysis.md) for the per-system delta.
2. **All six pipeline configs regenerated**: CMR, Dragen, HGDP, Microarray, ONT, PacBio (×2 callers). HGDP and Microarray keep their pipeline-specific subsetting.
3. **`exonic_annotation.hg{19,38}.BGStarget.txt`** extended with the 12 new gene rows.
4. **Reproducible regeneration toolchain** under [`data/source/v15/tools/`](data/source/v15/tools/) — fetches the V15 snapshot from the public ISBT REST API, rebuilds every `.dat`, validates against the C++ index loader, and writes the in-place swap with one-touch backup files.
5. **End-to-end regression harness** that builds the Docker image, runs the four bundled test samples, and emits a biological-equivalence verdict per system in [`data/source/v15/derived/regression_report.md`](data/source/v15/derived/regression_report.md).
6. **`Dockerfile.cn`** — single-stage build with Aliyun apt mirror + a patch dropping the NetBeans Makefile's hardcoded `-m64` flag so the binary compiles natively on Apple Silicon (linux/arm64). Use this if you are on a CN network and/or an ARM Mac.

The upstream C++ source is **untouched**; this is a pure data refresh + tooling addition.

---

## Quick start (Docker)

The fastest way to run bloodAGENT V15 on macOS / Linux:

```sh
# Clone with submodules
git clone --recurse-submodules https://github.com/gangchen/bloodAGENT.git
cd bloodAGENT

# Build the image — use Dockerfile.cn on Apple Silicon or CN networks
docker build -f Dockerfile.cn -t bloodagent:v15 .          # ~5–8 min on M-series Mac
# docker build              -t bloodagent:v15 .            # original upstream Dockerfile (x86_64 Linux)

# Smoke test
docker run --rm bloodagent:v15 --help

# Run the four bundled samples and write V15 JSONs to data/source/v15/derived/regression/
data/source/v15/tools/run_regression.sh bloodagent:v15

# Compare against the pre-V15 baselines bundled in data/testdata/
python3 data/source/v15/tools/analyze_regression_v2.py
```

---

## Project origin

bloodAGENT was created by the **IKMB Kiel** group (Wittig lab) as part of the work on resolving blood-group alleles from NGS / TGS data, with the original publication describing the cosine-similarity scoring approach and the HGDP benchmark. The upstream repository at https://github.com/ikmb/bloodAGENT remains the canonical source for the C++ engine.

This fork was started in May 2026 with one goal: keep the reference data layer aligned with the **ISBT Blood Group Database**, which replaced the per-system PDF allele tables in November 2025. Because the upstream `.dat` files were a PDF-era snapshot (2019–2024 vintage), they had drifted from current ISBT consensus by ~770 alleles by the time of V15. This fork picks up that drift and provides the toolchain to keep picking it up at each future ISBT release.

All credit for the algorithm, C++ implementation, and original benchmarking belongs to the upstream authors. This fork is data-only.

---

## Reference data: ISBT V15

| Metric | Value |
| --- | --- |
| Release | v15 (applied 2026-05-01, covers 2026-04 changes) |
| Source | https://blooddatabase.isbtweb.org/ (public REST API at `/api/`) |
| Blood group systems | 48 |
| Genes | 57 |
| Antigens | 397 |
| Effective alleles | 2 036 |
| Variants | 1 828 |
| Updated this release | 1 448 alleles touched |

**Newly added systems** (none of these existed in upstream bloodAGENT data): ER (044, *PIEZO1*), CD36 (045), ATP11C (046), MAL (047), plus older systems the upstream had not covered: CH_RG (017, *C4A/C4B*), KANNO (037, *PRNP*), SID (038, *B4GALNT2*), CTL2 (039, *SLC44A2*), PEL (040, *ABCC4*), MAM (041, *EMP3*), EMM (042, *PIGG*), ABCC1 (043). 82 new alleles across these 12 systems.

**Note on the `Phenotype_PDF_Table` column** in `genotype_to_phenotype_annotation_*.dat`: the column name is retained for backward compatibility with any downstream consumer, but it now holds the V15 `isbt_phenotype` value (the per-system PDF tables it used to reference were archived by ISBT in November 2025). Treat it as a phenotype label column going forward.

---

## V15 regression results

Full report in [`data/source/v15/derived/regression_report.md`](data/source/v15/derived/regression_report.md). Headline numbers from the four bundled test samples (HGDP00001/00003/00005 + NA24143):

| | Result |
| --- | --- |
| Schema-load failures on the new `.dat` files | **0** |
| Crashes / segfaults across 4 samples | **0** |
| **Biological equivalence vs pre-V15** (allele-level) | **120 / 123 system calls preserved (97.6 %)** |
| Phenotype calls equivalent after canonicalisation | 90 / 123 (73.2 %) |
| Genuine V15 reclassifications (need clinical review) | **3** — P1PK on HGDP00001/3, GLOB on HGDP00003 |

The report classifies every deviation into one of four categories:

- **A-renaming** — V15 renamed the identifier but the underlying biology is identical (e.g. upstream `GLOB*02` ≡ V15 `GLOB*01.02`; upstream `secretor` ≡ V15 `FUT2*01`). 0 prediction change.
- **B-obsoleted** — V15 retired an upstream allele identifier (`KLF1*BGM12` is `obsolete:true` in V15). V15 calls a phenotypically-equivalent successor.
- **C-reclassification** — V15 changed which DNA variants define an allele based on new evidence. The flagship case is **P1PK**: V15 redefined `A4GALT*02` (the P2 reference) as requiring the deep-intronic regulatory SNP `c.-188+3010G>T`, not the exonic `c.109A>G` (Met37Val). A sample typed only on exonic data can no longer be confidently called P1 vs P2 — that's a real consequence of the ISBT working party adopting newer evidence.
- **D-expansion** — V15's larger allele table means more alleles share the same SNV signature when the input VCF can't discriminate. Top systems on HGDP00001: RHD +98, XK +52, KEL +49, FUT2 +43, JK +40 extra alleles tied at top score. This is a precision effect, not a correctness regression — tune `--scoreRange` or use higher-resolution sequencing.

The full per-system verdict table sits in [`data/source/v15/derived/regression_report.md`](data/source/v15/derived/regression_report.md). The recommendation for clinical reports generated by this fork is: surface the V15 allele name as primary, append `(formerly: <upstream name>)` for systems with renaming, and flag P1PK calls as needing intronic-SNP or serology confirmation when c.109A>G is the only signal.

---

## Regenerating the data from a future ISBT release

When ISBT publishes V16 (or any later snapshot), the full regeneration takes about an hour mostly waiting on the rate-limited API:

```sh
# 1. Pull the V<N> snapshot (about 1 hour with the throttler).
data/source/v15/tools/fetch_isbt_v15.sh

# 2. Rebuild the master variation_annotation + master gt2pt
python3 data/source/v15/tools/build_variation_annotation.py \
    --raw-dir data/source/v15/raw \
    --allele-detail-dir data/source/v15/raw/alleles \
    --out-dir data/source/v15/derived \
    --old-master data/config/variation_annotation.dat

python3 data/source/v15/tools/build_gt2pt.py \
    --raw-dir data/source/v15/raw \
    --allele-detail-dir data/source/v15/raw/alleles \
    --out-dir data/source/v15/derived

# 3. Fork per-pipeline gt2pt files (HGDP, Microarray)
python3 data/source/v15/tools/build_pipeline_gt2pt.py

# 4. Generate exonic_annotation supplement for any new systems
python3 data/source/v15/tools/build_exonic_supplement.py

# 5. Validate everything against the C++ index loader semantics
python3 data/source/v15/tools/validate_dat.py \
    data/source/v15/derived/variation_annotation.v15.dat \
    data/source/v15/derived/genotype_to_phenotype_annotation.v15.dat

# 6. Swap the new files into data/config/ (creates *.pre-v15.bak backups)
data/source/v15/tools/apply_to_config.sh --apply

# 7. Rebuild the container and run regression
docker build -f Dockerfile.cn -t bloodagent:v16 .
data/source/v15/tools/run_regression.sh bloodagent:v16
python3 data/source/v15/tools/analyze_regression_v2.py
```

To roll back: `for f in $(find data/config -name '*.pre-v15.bak'); do mv "$f" "${f%.pre-v15.bak}"; done`.

---

## Introduction (from upstream)

**bloodAGENT** (Blood Antigen GENo Typer) is an open-source software tool designed for the determination of blood group alleles based on genetic markers. By analyzing genomic data from Next-Generation Sequencing (NGS) and Third-Generation Sequencing (TGS), bloodAGENT resolves blood group alleles and provides insights into genomic variations.

## Key Features
- **High accuracy** in allele determination under typical conditions.
- **Modular and flexible architecture**, allowing for future adaptations.
- **Uses cosine similarity scoring** to determine the best haplotype match.
- **Supports VCF and BigWig file formats** for variant and coverage data.
- **Open-source** for transparency and community collaboration.
- **(this fork)** Reference data pinned to **ISBT V15** (May 2026), with a reproducible upgrade path.

## System Requirements

**Supported platforms:** Compatible with Windows, macOS, and Linux through the Singularity / Docker image.

Native build dependencies (Linux only — see Apple Silicon caveat below):
- GCC or Clang compiler
- `libhts` library: https://github.com/samtools/htslib
- `libBigWig` library: https://github.com/dpryan79/libBigWig
- Python 3 (for parsing output files + running the V15 toolchain)
- https://github.com/mirror/tclap (CLI argument parsing)
- https://github.com/nlohmann/json (JSON output)

> **Apple Silicon (M-series Mac):** the upstream NetBeans Makefile passes `-m64` which g++ on ARM64 doesn't recognise. Use `Dockerfile.cn` (or apply the same `sed -i 's/=-m64$/=/' nbproject/Makefile-Release.mk` patch before `make`) to build on ARM. Native ARM build runs ~5–8 minutes; cross-build for `linux/amd64` via Docker Desktop emulation is ~3–5× slower.

## Installation

### Docker (recommended)
```sh
git clone --recurse-submodules https://github.com/gangchen/bloodAGENT.git
cd bloodAGENT
docker build -f Dockerfile.cn -t bloodagent:v15 .   # apply CN mirror + ARM64 patch
# or: docker build -t bloodagent:v15 .              # upstream Dockerfile, x86_64 Linux only
```

### Native (Linux x86_64)
```sh
sudo apt-get update && sudo apt-get install -y \
    g++ make zlib1g-dev libbz2-dev git liblzma-dev libcurl4-openssl-dev

git clone --recurse-submodules https://github.com/gangchen/bloodAGENT.git
cd bloodAGENT
cd external/htslib && make
cd ../libBigWig && make
cd ../..
make CONF=Release

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PWD/external/htslib:$PWD/external/libBigWig
./dist/Release/GNU-Linux/bloodAGENT --help
```

## Input Data Format

bloodAGENT requires two main input files:
- **VCF files**: genomic variants (compatible with hg19 and hg38).
- **BigWig files**: coverage data for SNVs that are not present in the VCF.

Configuration files (already provisioned in this fork at V15):
- `data/config/exonic_annotation.${build}.BGStarget.txt` — transcript annotation for blood-group targets (one per genome build).
- `data/config/<PIPELINE>/variation_annotation_<TAG>.dat` — variant annotation per pipeline.
- `data/config/<PIPELINE>/genotype_to_phenotype_annotation_<TAG>.dat` — genotype→phenotype mapping per pipeline.

Per-pipeline differences are primarily about how the RHCE *02 109 bp insertion is represented in the VCF (Illumina/ICA vs GATK vs sniffles vs pbsv); the V15 refresh preserves these handler files unchanged.

### Pipeline settings

- **HGDP** — original HGDP project pipeline. RHCE alleles are excluded from the gt2pt and detection is coverage-based via `-k/--trick`.
- **CMR** — Illumina ICA cloud constellation-mapped reads.
- **Dragen** — Illumina Dragen platform.
- **PacBio** — pbmm2 + GATK (or DeepVariant) for SNVs, pbsv for SVs.
- **ONT** — Oxford Nanopore + minimap2 + clair3 + sniffles (for the 109 bp insertion).
- **Microarray** — Affy / Illumina array. Subset of the gt2pt restricted to SNV / short-indel single-variant alleles.

Pipeline-selection details: [`data/config/`](data/config/).

### Variant Phasing

WhatsHap is recommended for long-read data. SHAPEIT5 is recommended for short-read HGDP-class data — but **SHAPEIT5 does not set the VCF Phasing-ID field**, which bloodAGENT requires. The `append_phasingID.py` script in this repo writes it for you.

Workflow for SHAPEIT5: phase a multi-sample VCF → split back into single-sample VCFs → run `append_phasingID.py` on each → feed to bloodAGENT.

### Testdata
Located under `./data/testdata/` — HGDP00001, HGDP00003, HGDP00005, and the GIAB NA24143. The bundled `*.phased.json` files are pre-V15 baseline outputs kept for the regression harness; the V15 outputs go to `data/source/v15/derived/regression/` (gitignored, reproducible). The complete HGDP benchmark dataset used in the upstream publication is available at https://www.internationalgenome.org/data-portal/data-collection/hgdp.

## Cosine Similarity Scoring

bloodAGENT uses **cosine similarity** to measure the similarity between observed haplotypes and reference haplotypes from ISBT. The score ranges from **0 to 2**:
- **1 per haplotype** is the theoretical maximum.
- **2** is the best possible match for a diploid genome.

Scores between different blood groups or individuals are **not directly comparable** — different systems have different numbers of relevant SNPs. A score of 1.9 vs 1.8 is meaningful only within the same system+allele pair.

> **V15 caveat:** the larger V15 allele table can produce more tied alleles at the top score on a given input VCF (e.g. RHD ties up to 99 alleles per haplotype on the bundled HGDP00001 sample). This is a precision effect of V15 having alleles whose discriminating variants aren't in the input VCF; tune `--scoreRange` or use higher-resolution sequencing if you need a single-allele call.

## Running bloodAGENT

### Job Type: Phenotype Analysis

Native binary:
```sh
bloodAGENT --job phenotype \
  --target ./data/config/exonic_annotation.hg38.BGStarget.txt \
  --variants ./data/config/HGDP/variation_annotation_HGDP.dat \
  --gt2pt ./data/config/HGDP/genotype_to_phenotype_annotation_HGDP.dat \
  --vcf ./data/testdata/HGDP00001/HGDP00001.phased.vcf.gz \
  --bigwig ./data/testdata/HGDP00001/HGDP00001.BGStarget.bw \
  --coverage 12 --verbose 2 --scoreRange 1 \
  --out HGDP00001.json \
  --build hg38 -k --id "HGDP00001"
```

Docker:
```sh
docker run --rm -v "$PWD":/work -w /work bloodagent:v15 \
  --job phenotype \
  --target /work/data/config/exonic_annotation.hg38.BGStarget.txt \
  --variants /work/data/config/HGDP/variation_annotation_HGDP.dat \
  --gt2pt /work/data/config/HGDP/genotype_to_phenotype_annotation_HGDP.dat \
  --vcf /work/data/testdata/HGDP00001/HGDP00001.phased.vcf.gz \
  --bigwig /work/data/testdata/HGDP00001/HGDP00001.BGStarget.bw \
  --coverage 12 --verbose 2 --scoreRange 1 \
  --out /work/HGDP00001.json \
  --build hg38 -k --id "HGDP00001"
```

Singularity (if you build a `.sif` from the docker image):
```sh
singularity exec bloodagent.sif /runtime/bloodAGENT --job phenotype <...same flags...>
```

### Job Type: Simulated Data Generation
```sh
bloodAGENT --job vcf \
  --variants ./data/config/variation_annotation.dat \
  --gt2pt ./data/config/genotype_to_phenotype_annotation.dat \
  -a "ABO*A1.01" -b "ABO*O.01.01" \
  --phased --dropout 1 --crack 5
```

### Command-Line Parameters

`--job phenotype`:

| Short | Long | Description | Type | Required | Default |
|---|---|---|---|---|---|
| `-j` | `--job phenotype` | Run phenotype determination. | String | Yes | – |
| `-t` | `--target <file>` | Transcript annotation for blood-group targets. | File | Yes | – |
| `-s` | `--variants <file>` | Variant annotation (V15 in this fork). | File | Yes | – |
| `-g` | `--gt2pt <file>` | Genotype→phenotype mapping (V15 in this fork). | File | Yes | – |
| `-v` | `--vcf <file>` | Phased VCF (comma-separate multiple files). | File | Yes | – |
| `-b` | `--bigwig <file>` | Coverage BigWig (or BAM). | File | No | – |
| `-c` | `--coverage <int>` | Minimum coverage for a reliable call. | Integer | No | `10` |
| `-d` | `--verbose <int>` | Verbosity 0–3. | Integer | No | `1` |
| `-r` | `--scoreRange <float>` | Score range below the top to report. | Float | No | – |
| `-o` | `--out <file>` | JSON output path. | File | No | `bloodAGENT.json` |
| `-u` | `--build <hg19\|hg38>` | Genome reference build. | String | Yes | – |
| `-k` | `--trick` | Use coverage-based RhD typing instead of variant-based. | Flag | No | `false` |
| `-f` | `--id <string>` | Sample identifier. | String | No | `unknown` |

`--job vcf`:

| Short | Long | Description | Type | Required | Default |
|---|---|---|---|---|---|
| `-j` | `--job vcf` | Generate a simulated VCF. | String | Yes | – |
| `-s` | `--variants <file>` | Variant annotation. | File | Yes | – |
| `-g` | `--gt2pt <file>` | gt→pt mapping. | File | Yes | – |
| `-a` | `--alleleA <string>` | First allele. | String | Yes | – |
| `-b` | `--alleleB <string>` | Second allele. | String | Yes | – |
| `-p` | `--phased` | Output phased haplotypes. | Flag | No | `false` |
| `-o` | `--dropout <float>` | SNP dropout probability. | Float | No | – |
| `-x` | `--crack <float>` | Haplotype-breakage probability at heterozygous sites. | Float | No | – |

## Output Format

`--job phenotype` writes JSON. The helper `deepBlood_values.py` extracts a tab-delimited summary. `--job vcf` writes VCF data lines to stdout.

### JSON structure

- **genome** — reference build (e.g. `hg38`).
- **sample_id** — value passed via `--id`.
- **version** — bloodAGENT binary version.
- **parameters** — full command-line snapshot.
- **loci** — blood-group systems.
  - **`<system>`**
    - **calls** — list of genotype/phenotype determinations:
      - **alleles** — detected V15 allele names per haplotype (plus per-allele coverage / mismatch issues).
      - **haplotypes** — genotype calls.
      - **phenotypes** — predicted phenotype(s) in V15 notation.
      - **score** — cosine similarity (0 if coverage-failed variants exist; otherwise up to 2).
      - **weak_score** — score ignoring coverage-failed variants.
      - **coverage_failed_variants** — variants with insufficient depth.
      - **mean_coverage** — coverage stats per region.
      - **relevant_variations** — every ISBT variant for the system (chrom, position, ref/alt, high-impact flag, depth, etc.).

> **V15 notation note:** phenotype strings now follow the ISBT-V15 convention (`Co(a+)` → `CO:1 or Co(a+)`). The biological meaning is unchanged; the dual notation is V15's canonical format. The regression harness includes a phenotype canonicaliser that handles this.

## How to Run Custom Secondary Analysis Scripts

This is the upstream-documented strategy for **RHCE antigens**, retained verbatim under V15.

### RHCE Antigen Detection Strategy

- **RHCE-Cc Antigen** — coverage of *RHCE* exon 2: absence ⇒ C-negative.
- **RHCE-Ee Antigen** — tagging SNP `676G>C`.

Both antigens are reported separately in the JSON output. The HGDP gt2pt at V15 still uses the upstream `RHCE*c / RHCE*C / RHCE*e / RHCE*E` pseudo-allele names so existing downstream consumers keep working.

### Step-by-step
1. Run `detect_RHCplusminus.py` on the sample → produces a `*.RHC.vcf`.
2. Pass it as a second `--vcf` to bloodAGENT (comma-separated):

```sh
--vcf sample.phased.vcf.gz,sample.RHC.vcf
```

## Special Case: RHD

Secondary analysis tools struggle to call RHD-specific variants reliably — even the complete `RHD*01N.01` deletion is often missed. bloodAGENT therefore offers a **coverage-based RHD detection** mode via `-k / --trick`:

- Homozygous `RHD*01N.01` deletion → **RhD negative**
- Any other allele combination → **RhD positive**

We **strongly recommend always enabling `-k`** unless you have an exceptionally reliable secondary-analysis pipeline.

## Limitations

- **Dropout**: missing variants degrade accuracy; ~50 % accuracy at a 50 % dropout rate.
- **Phasing**: marginal on most systems but important for KEL / ABO / Duffy.
- **Paralogous regions**: RHCE, GYPA/B/E, C4A/C4B, etc. can be mis-aligned or mis-called by upstream pipelines.
- **V15-specific**: P1PK typing on exonic-only data is now ambiguous (V15 moved the P2-defining variant to the deep intron `c.-188+3010G>T`). 12 new V15 systems have empty `cdsStart/exonStarts/exonEnds` in the supplementary exonic-annotation rows — pure SNV typing works, but coverage-based detection on those new systems requires UCSC refGene backfill. The TODO list is at [`data/source/v15/derived/exonic_annotation_TODO.tsv`](data/source/v15/derived/exonic_annotation_TODO.tsv).

## Licensing

BSD 2-Clause (unchanged from upstream). Third-party licenses are in `Third_Party_Licenses.md`; the docker / singularity image's `/licenses` directory contains per-component licenses for the bundled libraries.

Source code of this fork: https://github.com/gangchen/bloodAGENT.
Upstream source code: https://github.com/ikmb/bloodAGENT.
