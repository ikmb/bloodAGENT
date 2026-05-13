# bloodAGENT × ISBT V15 — 数据基线分析报告（最终版，全量数据）

> 数据基于 V15 release 全量拉取（2053 个等位基因详情已全部抓取）。
> Release: v15 applied 2026-05-01

## 1. ISBT V15 release 概览

| 指标 | V15 |
| --- | --- |
| Release | v15 (2026-05-01) |
| 覆盖周期 | 2026-04-01 ~ 2026-04-30 |
| 系统总数 | 55 (BGS: 48 + TF: 2 + Series: 2 + Collection: 3) |
| 基因数 | 57 |
| 抗原数 | 397 |
| 等位基因数（含 deleted/obsolete） | 2053 (实际有效 2036) |
| 变异数 | 1828 |
| 本期更新的等位基因 | 1448 |

## 2. ISBT release 历史

| Release | 日期 | 系统数 | 等位基因 | 抗原 |
| --- | --- | ---: | ---: | ---: |
| v2 | 2025-08-22 | 55 | 2007 | 396 |
| v3 | 2025-09-01 | 55 | 2007 | 396 |
| v8 | 2025-10-01 | 55 | 2007 | 397 |
| v9 | 2025-11-03 | 55 | 1965 | 397 |
| v10 | 2025-12-01 | 55 | 1971 | 397 |
| v11 | 2026-01-01 | 55 | 1971 | 397 |
| v12 | 2026-02-02 | 55 | 2000 | 397 |
| v13 | 2026-03-05 | 55 | 2002 | 397 |
| v14 | 2026-04-02 | 55 | 2002 | 397 |
| v15 | 2026-05-01 | 55 | 2004 | 397 |

## 3. 每系统对比：V15 vs bloodAGENT 当前 (HGDP)

| ISBT# | symbol | 主基因 | V15 alleles | bloodAGENT 当前 | delta | 备注 |
| --- | --- | --- | ---: | ---: | ---: | --- |
| 001 | ABO | ABO | 207 | 202 | +5 | +5 新增 |
| 002 | MNS | GYPA、GYPB、MNS_HYBRID | 106 | 30 | +76 | +76 新增 |
| 003 | P1PK | A4GALT | 43 | 38 | +5 | +5 新增 |
| 004 | RH | RHCE、RHD | 633 | 369 | +264 | +264 新增 |
| 005 | LU | BCAM | 37 | 26 | +11 | +11 新增 |
| 006 | KEL | KEL | 129 | 73 | +56 | +56 新增 |
| 007 | LE | FUT3 | 62 | 52 | +10 | +10 新增 |
| 008 | FY | ACKR1 | 32 | 24 | +8 | +8 新增 |
| 009 | JK | SLC14A1 | 81 | 31 | +50 | +50 新增 |
| 010 | DI | SLC4A1 | 25 | 21 | +4 | +4 新增 |
| 011 | YT | ACHE | 5 | 2 | +3 | +3 新增 |
| 012 | XG | CD99、XG | 8 | 2 | +6 | +6 新增 |
| 013 | SC | ERMAP | 10 | 8 | +2 | +2 新增 |
| 014 | DO | ART4 | 24 | 15 | +9 | +9 新增 |
| 015 | CO | AQP1 | 11 | 9 | +2 | +2 新增 |
| 016 | LW | ICAM4 | 6 | 3 | +3 | +3 新增 |
| 017 | CH_RG | C4A、C4B | 14 | 0 | +14 | **完全缺失** |
| 018 | H | FUT1、FUT2 | 122 | 58 | +64 | +64 新增 |
| 019 | XK | XK | 68 | 25 | +43 | +43 新增 |
| 020 | GE | GYPC | 21 | 9 | +12 | +12 新增 |
| 021 | CROM | CD55 | 24 | 18 | +6 | +6 新增 |
| 022 | KN | CR1 | 11 | 7 | +4 | +4 新增 |
| 023 | IN | CD44 | 6 | 5 | +1 | +1 新增 |
| 024 | OK | BSG | 4 | 4 | +0 | 对齐 |
| 025 | RAPH | CD151 | 5 | 4 | +1 | +1 新增 |
| 026 | JMH | SEMA7A | 8 | 6 | +2 | +2 新增 |
| 027 | I | GCNT2 | 14 | 11 | +3 | +3 新增 |
| 028 | GLOB | B3GALNT1 | 15 | 14 | +1 | +1 新增 |
| 029 | GIL | AQP3 | 3 | 2 | +1 | +1 新增 |
| 030 | RHAG | RHAG | 48 | 32 | +16 | +16 新增 |
| 031 | FORS | GBGT1 | 6 | 5 | +1 | +1 新增 |
| 032 | JR | ABCG2 | 38 | 29 | +9 | +9 新增 |
| 033 | LAN | ABCB6 | 47 | 42 | +5 | +5 新增 |
| 034 | VEL | SMIM1 | 7 | 4 | +3 | +3 新增 |
| 035 | CD59 | CD59 | 8 | 4 | +4 | +4 新增 |
| 036 | AUG | SLC29A1 | 5 | 3 | +2 | +2 新增 |
| 037 | KANNO | PRNP | 2 | 0 | +2 | **完全缺失** |
| 038 | SID | B4GALNT2 | 5 | 0 | +5 | **完全缺失** |
| 039 | CTL2 | SLC44A2 | 5 | 0 | +5 | **完全缺失** |
| 040 | PEL | ABCC4 | 5 | 0 | +5 | **完全缺失** |
| 041 | MAM | EMP3 | 6 | 0 | +6 | **完全缺失** |
| 042 | EMM | PIGG | 9 | 0 | +9 | **完全缺失** |
| 043 | ABCC1 | ABCC1 | 2 | 0 | +2 | **完全缺失** |
| 044 | ER | PIEZO1 | 7 | 0 | +7 | **完全缺失** |
| 045 | CD36 | CD36 | 23 | 0 | +23 | **完全缺失** |
| 046 | ATP11C | ATP11C | 2 | 0 | +2 | **完全缺失** |
| 047 | MAL | MAL | 2 | 0 | +2 | **完全缺失** |
| 048 | PIGZ | - | 0 | 0 | +0 |  |

**合计：V15 = 2036 等位基因 / bloodAGENT HGDP = 1267 等位基因 / delta = +769**

## 4. V15 中 bloodAGENT 完全缺失的系统（12 个）

| ISBT# | symbol | 主基因 | V15 等位基因数 |
| --- | --- | --- | ---: |
| 017 | CH_RG | C4A、C4B | 14 |
| 037 | KANNO | PRNP | 2 |
| 038 | SID | B4GALNT2 | 5 |
| 039 | CTL2 | SLC44A2 | 5 |
| 040 | PEL | ABCC4 | 5 |
| 041 | MAM | EMP3 | 6 |
| 042 | EMM | PIGG | 9 |
| 043 | ABCC1 | ABCC1 | 2 |
| 044 | ER | PIEZO1 | 7 |
| 045 | CD36 | CD36 | 23 |
| 046 | ATP11C | ATP11C | 2 |
| 047 | MAL | MAL | 2 |

**12 个缺失系统，合计 82 个新等位基因需要新增。**

## 5. 等位基因级别 diff 汇总

- **添加** (V15 有但 HGDP 没有的 allele)：887
- **撤销** (HGDP 有但 V15 没有的 allele)：118
- **表型变更** (同名 allele 但表型字符串不同)：443
- 完整列表：`data/source/v15/derived/allele_diff.tsv`

## 6. 数据获取通路（已打通并跑通）

- REST API base：`https://blooddatabase.isbtweb.org/api/` （NestJS，无需认证）
- 关键端点：
  - `GET /api/release/15` — V15 元信息 + `updatedAlleleIds` (1448 个本期更新)
  - `GET /api/system|gene|antigen|allele|variant|phenotype` — 列表（bulk allele **不含** variants[]）
  - `GET /api/allele/<id>` — **核心**：单 allele 详情含 `variants[]` 嵌套数组
- 限流：实测约 1 req/s 持续，HTTP 429 后退避。完整 2053 个 allele 详情 ~1 小时。
- 一次性脚本：`data/source/v15/tools/fetch_isbt_v15.sh`

## 7. 转换规则（已实现于 build_*.py）

| bloodAGENT 列 | V15 字段来源 | 处理 |
| --- | --- | --- |
| system/gene (vanno col 1) | variant.gene_name | 直接复用 |
| Transcript annotation (vanno col 2) | variant.hgvs_transcript | 直接 |
| Transcript annotation short (vanno col 3) | 取 hgvs_transcript 中 `:c.` 之后部分 | 自动 strip `c.` 前缀 |
| Chr/Pos/Coord-VCF (vanno) | variant.grch37/38_chr/pos | 'chr' 前缀按列区分 |
| Ref/Alt (vanno) | variant.grch37/38_ref/alt | 直接 |
| strand (vanno col 10/18) | gene.strand（V15 未暴露） | **TODO**：从现有 master 表 backport 或 RefSeq 查询 |
| Reference base (vanno col 12/19) | 取 VCF ref，负链 revcomp，ins/del 去锚定碱基 | C++ 约定 |
| TYPE (vanno col 28) | 由 ref/alt 长度推断 | SNV/ins/del/delins |
| MySystemKey (gt2pt col 2) | system.symbol | RH/MNS 拆为 gene 名以兼容现有 |
| Phenotype_PDF_Table/Phenotype (gt2pt col 5/6/7) | allele.isbt_phenotype | HTML 已清洗，空白合并 |
| base_change (gt2pt col 8) | join(variant.hgvs_transcript 短形式) | 空格分隔 |
| acid_change (gt2pt col 9) | join(variant.hgvs_predicted_protein 短形式) | |
| incidence (gt2pt col 10) | max(variant.gnomad_all) × 100% | gnomAD 替代 PDF 时代临床频率 |

## 8. 已知 / 待手工处理的问题

1. **gene.strand 信息**：V15 `/api/gene` 不返回链方向，新系统的 strand 需要 backport。已在 build_variation_annotation.py 标记为 TODO。
2. **大 SV 等位基因**（26 个）：无 `hgvs_transcript`（多 kb 缺失），已主动跳过并落入 `missing_coords_report.tsv`。bloodAGENT 通过 coverage 检测这类 SV（如 RHD 全基因缺失），不收录是正确做法。
3. **无 variants[] 的等位基因**（113 个）：在 V15 API 中只有 isbt_allele 名称但 variants 数组为空。包含参考等位基因（如 ABO*A1.01）与一些占位条目。这些不进入 vanno 表但需要进入 gt2pt 表（已正确处理）。
4. **RHCE 109bp 插入**：5 个 pipeline-specific overlay 文件（CMR/Dragen/ONT/PacBio-GATK/PacBio-pbsv）保留手工维护，**不**由 V15 数据再生。
5. **exonic_annotation.{hg19,hg38}.BGStarget.txt**：12 个新系统需追加外显子区间，从 `gene.json` 衍生（待实现）。
6. **CH_RG (system 017)**：基因 C4A+C4B 高度同源 + CNV 多态，短读 NGS 难解析；建议初版跳过。
7. **incidence 列语义变更**：PDF 时代为临床频率，V15 是 gnomAD AF，两者不可比。**建议**：对参考等位基因保留旧 incidence。

## 9. 输出产物清单

| 文件 | 行数 | 状态 |
| --- | ---: | --- |
| `data/source/v15/raw/*.json` | - | ✅ V15 全量原始数据 |
| `data/source/v15/raw/alleles/*.json` | 2053 | ✅ 每 allele 详情（含 variants[]） |
| `data/source/v15/derived/variation_annotation.v15.dat` | 1767+1 | ✅ schema 合规（CISBTAnno 可直接读） |
| `data/source/v15/derived/genotype_to_phenotype_annotation.v15.dat` | 2036+1 | ✅ schema 合规 |
| `data/source/v15/derived/allele_diff.tsv` | 1448 | ✅ 887 added + 443 phen_changed + 118 retired |
| `data/source/v15/derived/system_summary.tsv` | 49 | ✅ per-system 增量 |
| `data/source/v15/derived/missing_coords_report.tsv` | 26 | ✅ 需 liftOver 的 SV |
| `data/source/v15/derived/unparseable_alleles.tsv` | 113 | ✅ 无 variants[] 的等位基因 |

## 10. 下一步（实施阶段，未启动）

- [ ] 用 `variation_annotation.v15.dat` 替换 `data/config/variation_annotation.dat`，保留 strand 等列从老表 backport 缺失字段
- [ ] 用新 gt2pt 替换 6 个 pipeline 的 `genotype_to_phenotype_annotation_*.dat`（RHCE 在 HGDP 仍要保留删除逻辑）
- [ ] 5 个 pipeline-specific variation_annotation overlay 保留不动（手工维护的 RHCE 109bp 插入）
- [ ] 12 个新系统补 exonic_annotation 区间（从 gene.json 自动衍生）
- [ ] 4 个测试样本 (HGDP00001/00003/00005/NA24143) 跑回归
- [ ] 重建 Singularity 镜像 + 更新 README 中 ISBT 版本声明
- [ ] 视情况向 ikmb/bloodAGENT 提 PR
