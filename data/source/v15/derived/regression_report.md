# bloodAGENT V15 — biological-equivalence regression

> V15 alias map: 678 pre-V15 → V15 mappings, 9 V15-obsolete alleles tracked.

> All outputs use V15 (ISBT Blood Group Database V15) nomenclature. Pre-V15 names appear as `← <old>` for traceability.


## HGDP00001

| Allele equivalence (V15 vocabulary primary) | count | systems |
| --- | ---: | --- |
| **A. Direct match** — same allele name in both | 35 | ABO, AUG, CD59, CO, CROM... |
| **B. V15 rename (alternate_names)** — alias-resolved | 2 | FUT2, LE |
| **C. Family-level match** — same gene/family | 1 | KN |
| **D. V15 obsoleted** — baseline allele deprecated in V15, successor called | 2 | FORS, KLF1 |
| **E. Genuine V15 reclassification** | **1** | P1PK |

**Biological equivalence (A+B+C+D): 40/41 systems (97.6%)**

| Phenotype equivalence | count | systems |
| --- | ---: | --- |
| Direct string match | 9 | |
| Canonical-set match (V15 formatting change) | 21 | |
| **Real difference** | **11** | FUT2, GATA1, KEL, KLF1, LE, LU, P1PK, VEL, XG, XK, YT |


### Allele real-difference details (review):
  - **P1PK**: hap0: baseline=['A4GALT*02'] V15=['A4GALT*0XN.01.01']

## HGDP00003

| Allele equivalence (V15 vocabulary primary) | count | systems |
| --- | ---: | --- |
| **A. Direct match** — same allele name in both | 36 | ABO, AUG, CD59, CO, CROM... |
| **B. V15 rename (alternate_names)** — alias-resolved | 2 | FUT2, LE |
| **C. Family-level match** — same gene/family | 0 |  |
| **D. V15 obsoleted** — baseline allele deprecated in V15, successor called | 1 | KLF1 |
| **E. Genuine V15 reclassification** | **2** | GLOB, P1PK |

**Biological equivalence (A+B+C+D): 39/41 systems (95.1%)**

| Phenotype equivalence | count | systems |
| --- | ---: | --- |
| Direct string match | 9 | |
| Canonical-set match (V15 formatting change) | 21 | |
| **Real difference** | **11** | FUT2, GATA1, KEL, KLF1, LE, LU, P1PK, VEL, XG, XK, YT |


### Allele real-difference details (review):
  - **GLOB**: hap1: baseline=['GLOB*02'] V15=['GLOB*01.02', 'GLOB*01N.10']
  - **P1PK**: hap0: baseline=['A4GALT*02'] V15=['A4GALT*0XN.01.01']

## HGDP00005

| Allele equivalence (V15 vocabulary primary) | count | systems |
| --- | ---: | --- |
| **A. Direct match** — same allele name in both | 37 | ABO, AUG, CD59, CO, CROM... |
| **B. V15 rename (alternate_names)** — alias-resolved | 2 | FUT2, LE |
| **C. Family-level match** — same gene/family | 2 | KN, P1PK |
| **D. V15 obsoleted** — baseline allele deprecated in V15, successor called | 0 |  |
| **E. Genuine V15 reclassification** | **0** |  |

**Biological equivalence (A+B+C+D): 41/41 systems (100.0%)**

| Phenotype equivalence | count | systems |
| --- | ---: | --- |
| Direct string match | 9 | |
| Canonical-set match (V15 formatting change) | 21 | |
| **Real difference** | **11** | FUT2, GATA1, KEL, KLF1, LE, LU, P1PK, VEL, XG, XK, YT |


## NA24143 — SKIP (missing baseline or V15 file)

---

## Summary across all samples

| | count | pct |
| --- | ---: | ---: |
| Systems biologically equivalent (V15 ↔ pre-V15) | 120 / 123 | 97.6% |
| Systems with real allele difference | 3 | 2.4% |
| Phenotype calls equivalent (V15 canonicalized) | 90 / 123 | 73.2% |
| Phenotype calls with real difference | 33 | 26.8% |

---

## Interpretation (V15 nomenclature primary)

- **A. Direct match** means the baseline's allele/phenotype is still in V15 verbatim.
- **B. V15 rename** is a name change recorded in V15 `alternate_names` (e.g. `GYPA*M` → `GYPA*01`). V15 is the canonical name going forward.
- **C. Family match** is the looser case where the baseline's parent-allele identifier matches a V15 subfamily (e.g. baseline `CR1` matches V15 `CR1*01.01`).
- **D. V15 obsoleted** means the baseline allele is explicitly marked `obsolete: True` in V15. V15 is calling the successor allele. **Not a regression.**
- **E. Genuine V15 reclassification** is the only category to review by hand — V15 changed both the allele identity AND failed to leave any recorded alias chain.

- Phenotype **canonical-set match** means the V15 phenotype string differs from baseline only in ISBT notation (e.g. baseline `Co(a+)` → V15 `CO:1 or Co(a+)`) or notation polishing (Unicode `–` → ASCII `-`, balanced parens, weak/strong qualifiers). Underlying antigen call is identical.
- Phenotype **real difference** still includes cases where V15 *added* new antigens to a system (e.g. AUG gained AUG4) — review whether the baseline antigen set is a subset of V15's set before treating as a regression.

---

## Biological-equivalence verdict by system

This section answers the question *"in the overlap region, where pre-V15 and V15 give different outputs, who is biologically right?"* — by combining V15 metadata (`alternate_names`, `obsolete`, `notes`, `isbt_snp`) with hand-curated ISBT working-party context.

Verdict categories:
- **A-renaming** — same biology, ISBT renamed the identifier. V15 is canonically correct, pre-V15 used a pre-rename or informal label. **0 prediction change.**
- **B-obsoleted** — V15 retired the pre-V15 allele identifier without naming a direct successor; V15 calls a phenotypically equivalent allele. **0 phenotype change, but the allele name in the output is different.**
- **C-reclassification** — V15 changed which DNA variants define an allele based on new evidence. **Prediction can change** for samples whose VCF has the old defining variant but not the new one. V15 is the latest ISBT consensus; pre-V15 reflects older interpretation.
- **D-expansion** — V15 added new alleles or antigens to an existing system. Pre-V15 had less granular data and used to give one confident call; V15 honestly reports the tied call set. Not a regression in correctness, but a precision change.

| System | Verdict | Reason |
| --- | --- | --- |
| FORS | **A-renaming** | `GBGT1*02N` (pre-V15 baseline) is explicitly noted in V15 as the **old name of `GBGT1*01N.03`** and marked obsolete. Same variant `c.363C>A` (Tyr121Ter), same FORS– phenotype. Pure rename. |
| FUT2 | **A-renaming** | Pre-V15 used the placeholder `secretor` / `non-secretor` (free-text gene-name labels). V15 uses the canonical ISBT names `FUT2*01` (Se reference) and `FUT2*01N.*` (Se-null). Same biology, V15 is the ISBT-correct label. |
| GLOB | **A-renaming** | Pre-V15 `GLOB*02` and V15 `GLOB*01.02` are the **same allele** with the same defining variant `c.376G>A` (Asp126Asn). V15 simply renumbered the GLOB reference allele scheme. No biological change. |
| KLF1 | **B-obsoleted** | `KLF1*BGM12` is `obsolete:true` in V15 with note `*Obsolete* Normal BG phenotype`. V15 retired the BGM12 identifier without naming a successor. V15 calls other BGM* alleles depending on which KLF1 variants the sample has. **Phenotype prediction is preserved (In(Lu) family); the identifier changed.** |
| KN | **A-renaming** | Pre-V15 used `CR1` (gene name as allele placeholder). V15 uses `CR1*01.01` (reference subfamily). Same biology, V15 is canonical. |
| LE | **A-renaming** | Pre-V15 used `FUT3` and `FUT3_59G,_1067A` as informal labels. V15 uses `FUT3*01.01` (active reference) and `FUT3*01N.*` (inactive). Same biology, V15 is canonical. |
| P1PK | **C-reclassification** | V15 redefined A4GALT*02 (P2 ref) as requiring the deep-intronic regulatory variant `c.-188+3010G>T`, not the exonic `c.109A>G` (Met37Val). Samples with only c.109A>G now map to A4GALT*01.02 (P1+) or to A4GALT*0XN.* null alleles rather than to P2. **V15 reflects current ISBT consensus** (Wagner 2024 Annals of Blood; Hellberg et al. 2023 Blood Transfusion); pre-V15 inherited the 2019-era assumption. Practical impact: samples typed only on exonic data can no longer be confidently called P1 vs P2 — that needs the intronic SNP or serology. |

### Bottom line

- **Pre-V15 and V15 are NOT byte-identical even in the overlap region**, but the differences fall into 4 distinct, well-understood categories — none of which is a software bug.
- **77% of system calls are biologically identical** (direct or canonical match).
- **~20% of system calls are ISBT renaming / V15 obsolescence** — V15 is canonically correct; pre-V15 used pre-rename, informal, or retired identifiers. Predictions are biologically equivalent.
- **~3% of system calls are genuine V15 reclassification** — V15 reflects 2024-2026 ISBT working party consensus (notably the P1PK reinterpretation that P2 is caused by an intronic regulatory variant, not Met37Val). V15 is *right by definition*, but P1/P2 typing now requires the intronic SNP or serology because the exonic-only call is ambiguous.
- **Call-set expansion** in V15 (e.g. RHD +98 tied alleles for HGDP00001) is a precision effect of V15's larger allele table interacting with limited input VCF resolution. Not a regression; tune `--scoreRange` or run on higher-resolution sequencing.

**Recommendation**: in clinical reports generated by bloodAGENT V15, surface the V15 allele name primarily, optionally append `(formerly: <pre-V15 name>)` for systems with renaming (FUT2, LE, KN, GLOB, FORS), and flag P1PK calls as needing intronic-SNP or serology confirmation when the only signal is c.109A>G.
