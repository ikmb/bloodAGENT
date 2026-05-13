# ISBT V15 系统 ↔ bloodAGENT 数据 — 系统名映射

> 数据源：`raw/system.json`（来自 `GET https://blooddatabase.isbtweb.org/api/system`，55 条目，含 48 个血型系统 + 2 转录因子 + 3 collection + 2 series）。
> 当前数据：`data/config/variation_annotation.dat` 与 6 个 pipeline 的 .dat。

## ISBT V15 全部 48 个血型系统（按 isbt_number 排序）

| ISBT# | DB id | symbol | name | bloodAGENT 当前是否覆盖 | 备注 |
|---|---|---|---|---|---|
| 1 | 1 | ABO | ABO | ✅ ABO | |
| 2 | 2 | MNS | MNS | ⚠️ 拆为 GYPA + GYPB | 缺 GYPE / MNS_HYBRID |
| 3 | 3 | P1PK | P1PK | ✅ P1PK | gene = A4GALT |
| 4 | 4 | RH | Rh | ✅ RHD + RHCE | bloodAGENT 又把 RHCE 在 gt2pt 拆为 RHC/RHE |
| 5 | 5 | LU | Lutheran | ✅ LU | gene = BCAM |
| 6 | 6 | KEL | Kell | ✅ KEL | |
| 7 | 7 | LE | Lewis | ✅ LE | gene = FUT3 |
| 8 | 8 | FY | Duffy | ✅ FY | gene = ACKR1 |
| 9 | 9 | JK | Kidd | ✅ JK | gene = SLC14A1 |
| 10 | 10 | DI | Diego | ✅ DI | gene = SLC4A1 |
| 11 | 11 | YT | Yt | ✅ YT | gene = ACHE |
| 12 | 12 | XG | Xg | ✅ XG | |
| 13 | 13 | SC | Scianna | ✅ SC | gene = ERMAP |
| 14 | 14 | DO | Dombrock | ✅ DO | gene = ART4 |
| 15 | 15 | CO | Colton | ✅ CO | gene = AQP1 |
| 16 | 16 | LW | Landsteiner-Wiener | ✅ LW | gene = ICAM4 |
| 17 | 17 | CH_RG | Chido/Rodgers | ❌ **缺** | 多基因 C4A/C4B (复杂 CNV) |
| 18 | 18 | H | H | ✅ H | gene = FUT1（bloodAGENT 还有独立 FUT2 条目） |
| 19 | 19 | XK | Kx | ✅ XK | |
| 20 | 20 | GE | Gerbich | ✅ GE | gene = GYPC |
| 21 | 21 | CROM | Cromer | ✅ CROM | gene = CD55 |
| 22 | 22 | KN | Knops | ✅ KN | gene = CR1 |
| 23 | 23 | IN | Indian | ✅ IN | gene = CD44 |
| 24 | 24 | OK | Ok | ✅ OK | gene = BSG |
| 25 | 25 | RAPH | Raph | ✅ RAPH | gene = CD151 |
| 26 | 26 | JMH | John Milton Hagen | ✅ JMH | gene = SEMA7A |
| 27 | 27 | I | I | ✅ I | gene = GCNT2 |
| 28 | 28 | GLOB | Globoside | ✅ GLOB | gene = B3GALNT1 |
| 29 | 29 | GIL | Gill | ✅ GIL | gene = AQP3 |
| 30 | 30 | RHAG | Rh-associated glycoprotein | ✅ RHAG | |
| 31 | 31 | FORS | FORS | ✅ FORS | gene = GBGT1 |
| 32 | 32 | JR | JR | ✅ JR | gene = ABCG2 |
| 33 | 33 | LAN | LAN | ✅ LAN | gene = ABCB6 |
| 34 | 34 | VEL | Vel | ✅ VEL | gene = SMIM1 |
| 35 | 35 | CD59 | CD59 | ✅ CD59 | |
| 36 | 36 | AUG | Augustine | ✅ AUG | gene = SLC29A1 |
| 37 | 37 | KANNO | Kanno | ❌ **缺** | gene = PRNP |
| 38 | 38 | SID | SID | ❌ **缺** | gene = B4GALNT2 |
| 39 | 39 | CTL2 | CTL2 | ❌ **缺** | gene = SLC44A2 |
| 40 | 40 | PEL | PEL | ❌ **缺** | gene = ABCC4 |
| 41 | 41 | MAM | MAM | ❌ **缺** | gene = EMP3 |
| 42 | 42 | EMM | EMM | ❌ **缺** | gene = PIGG |
| 43 | 43 | ABCC1 | ABCC1 | ❌ **缺** | gene = ABCC1 |
| 44 | 44 | ER | Er | ❌ **缺** | gene = PIEZO1 |
| 45 | 45 | CD36 | CD36 | ❌ **缺** | |
| 46 | 46 | ATP11C | ATP11C | ❌ **缺** | |
| 47 | 47 | MAL | MAL | ❌ **缺** | |
| 48 | 48 | PIGZ | PIGZ | ❌ **缺** | gene = PIGZ |

**结论：bloodAGENT 缺 13 个血型系统**：CH_RG (017)、KANNO (037)、SID (038)、CTL2 (039)、PEL (040)、MAM (041)、EMM (042)、ABCC1 (043)、ER (044)、CD36 (045)、ATP11C (046)、MAL (047)、PIGZ (048)。

## V15 中非系统类（不计入 48 之列）

| ISBT# | DB id | symbol | category | 备注 |
|---|---|---|---|---|
| 101 | 46 | KLF1 | Transcription Factor | bloodAGENT 已收（影响 LU 等系统） |
| 102 | 45 | GATA1 | Transcription Factor | bloodAGENT 已收 |
| 207 | 61 | 207 | Collection (Ii) | 未收 |
| 210 | 62 | 210 | Collection (Le) | 未收 |
| 213 | 60 | 213 | Collection (MN CHO) | 未收 |
| 700 | 58 | 700 | Series (Low Prevalence Antigens) | 未收 |
| 901 | 57 | 901 | Series (High Prevalence Antigens) | 未收 |

## bloodAGENT `system/gene` 列与 ISBT V15 的命名差异

| bloodAGENT 写法 | ISBT V15 system symbol | 说明 |
|---|---|---|
| `GYPA` | `MNS` (gene GYPA) | bloodAGENT 用 gene 名 |
| `GYPB` | `MNS` (gene GYPB) | bloodAGENT 用 gene 名；缺 GYPE |
| `FUT2` | `H` (gene FUT2 = Secretor) | bloodAGENT 单独列出，ISBT 归 H 系统 |
| `RHC` (仅 gt2pt) | `RH` (gene RHCE) | bloodAGENT 把 RHCE 拆为 RHC + RHE 伪系统记 phenotype |
| `RHE` (仅 gt2pt) | `RH` (gene RHCE) | 同上 |
| `RHCE` (variation_annotation) | `RH` (gene RHCE) | 与 ISBT 一致 |

**升级策略建议**：保持 bloodAGENT 现有 `system/gene` 列内容为"gene name"（与代码 `CIsbtGt2Pt` 的索引键一致，破坏成本太高），新增系统按其主基因名（如 PIEZO1 而非 ER）作为 `system/gene` 值。V15 的 `gene.name` 字段即可直接复用。

## 13 个新系统的主基因（用于 bloodAGENT system/gene 列）

| ISBT V15 symbol | 主基因（拟用作 bloodAGENT `system/gene`） | 染色体 |
|---|---|---|
| CH_RG | C4A 或 C4B（CNV 复杂，可能延后） | chr6 |
| KANNO | PRNP | chr20 |
| SID | B4GALNT2 | chr17 |
| CTL2 | SLC44A2 | chr19 |
| PEL | ABCC4 | chr13 |
| MAM | EMP3 | chr19 |
| EMM | PIGG | chr4 |
| ABCC1 | ABCC1 | chr16 |
| ER | PIEZO1 | chr16 |
| CD36 | CD36 | chr7 |
| ATP11C | ATP11C | chrX |
| MAL | MAL | chr2 |
| PIGZ | PIGZ | chr3 |

> 染色体取自 `raw/gene.json`，待数据拉取脚本汇总到 `derived/genes_for_new_systems.tsv`。
