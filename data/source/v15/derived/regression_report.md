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
