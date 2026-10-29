# 01b_format_hayles_phenotypes.py — Hayles 2013 表型表

当前流程的第一步。以 Hayles 2013 为主表构建基因级表型表。Kim 2010 只提供系统名和 Paralog 列。Kim 表型数据不合并。

## 流程位置

`00 → 01b → 03 → 04 → 05 → {02, 06, 07}`

`01b` 的输出 `Hayles2013_phenotypes_annotated.xlsx` 是 `03` 的输入。

## 逻辑

1. 读取 Kim 的 `Systemic ID` 和 `Paralog`。
2. 用当前 PomBase 注释更新 Kim 主 ID 和 Paralog ID。
3. 生成 Kim paralog pre-merge，带 `Kim_Note`。
4. 读取 Hayles 的 4 列，改成清晰列名。
5. 清洗 Hayles 表型文本：去空白、去句尾句号、把 `25, 32` / `25 32` / `32, 25` 统一为 `25,32`。
6. 更新 Hayles ID，生成 Hayles pre-merge，带 `Hayles_Note`。
7. 以 Hayles pre-merge 为主表，左连接 Kim paralog pre-merge。
8. 按更新后的 `Systematic ID` 加 `Gene name`、`Gene product`、`Gene type`。`Gene name` 为空时用 `Systematic ID`。

## Hayles 列名

| 原列名 | 新列名 |
|---|---|
| `Deletion mutant phenotype description` | `Phenotype_description` |
| `Phenotypic classification used for analysis` | `Morphological_classification` |
| `Gene dispensability. This study` | `Gene_dispensability` |

## 输出

| 文件 | 行数 | 列数 |
|---|---:|---:|
| `data/1_formatted/Kim2010_paralog_premerge.xlsx` | 4,913 | 4 |
| `data/1_formatted/Hayles2013_phenotypes_premerge.xlsx` | 4,843 | 6 |
| `data/1_formatted/Hayles2013_phenotypes_annotated.xlsx` | 4,843 | 12 |

## 下游使用的列

| 脚本 | 使用的列 |
|---|---|
| `03_group_genes.py` | `Systematic ID`, `Phenotype_description` |
| `04_categorize_phenotypes.py` | `Phenotype_description`, `Morphological_classification`, `Gene_dispensability`, `Basic phenotype` |
| `02_keyword_profile.py` | `Phenotype_description` |
| `06_plot_dr_distribution.py` | `Systematic ID`, `Sub_category`, `Category` |
| `07_plot_dr_distribution_pure.py` | `Systematic ID`, `Sub_category`, `Morphological_classification` |

## 共享模块

- `src/id_mapping.py` — PomBase ID 更新、Paralog 处理、note 合并。
- `src/phenotype.py` — 列名常量和表型文本清洗。
- `src/dataset_merge.py` — Kim/Hayles 原始表读取和重复 ID 合并。
- `src/pipeline_utils.py` — 日志和 DIT-HAP 工具。

## 运行

```bash
mamba run -n bioinformatics python scripts/01b_format_hayles_phenotypes.py
```

## 说明

- Kim 和 Hayles 的表型数据不合并。
- 不使用 Hayles 的 `Paralogues` 列。
- `SPAC27D7.10c` 折叠进 `SPAC27D7.09c`。`Kim_ID` 写两个 ID。
- `SPAC823.02` 更新为 `SPAC823.17`。
- `SPAC1556.06.1` 更新为 `SPAC1556.06`。`note` 仍是 `NO_CHANGE`（已知问题）。
- `SPBC8E4.02c` 保留旧 ID。`Gene name`、`Gene product`、`Gene type` 来自 synonym `SPNCRNA.9001`。
- 5 个 `RECLASSIFIED` 行：4 个 pseudogene，1 个 lncRNA。Hayles 表中有 3 个 pseudogene 和 1 个 lncRNA。
- 旧 `01_format_and_update_ids.py` 已移到 `tmp/01_format_and_update_ids/`。它生成 25 列旧表和 changelog，不再进入当前流程。
