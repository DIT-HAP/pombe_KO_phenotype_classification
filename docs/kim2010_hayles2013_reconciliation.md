# Kim 2010 vs Hayles 2013 deletion datasets: reconciliation and data-quality notes

> **Archive note (2026-10-29).** This note describes the archived
> `01_format_and_update_ids.py` union, now in `tmp/01_format_and_update_ids/`.
> The live pipeline starts at `01_format_hayles_phenotypes.py` and does not
> build the union. The per-gene facts below still apply. Columns named
> `Hayles2013_*` belong to the archived schema; the live `01` table calls the
> study call `Gene_dispensability`.

This note records the reconciliation behind the question: **why do some genes
have a 2010 dispensability result but no 2013 result?** It is the analytical
companion to the machine-readable union produced by the archived
`01_format_and_update_ids.py`; the former `00b_compare_kim2010_hayles2013.py`
script has been retired because that union now carries its per-gene results.

The direct answer is that **11 constructed Kim 2010 strains are absent from
the 4,843-gene Hayles 2013 screen**. The paper does not name them or give a
per-gene reason; they are the arithmetic remainder of the paper's own
statement that 4,825 of Kim's strains were used.

The 2010 input is the **full Supplementary Table 1**: 4,914 candidate genes,
of which 4,836 were constructed and have a V/E dispensability call. The 2013
input has 4,843 genes. The reconciliation is exact:

```
4,836 constructed Kim strains
-   11 not included in the 2013 screen
+   18 Hayles additions
= 4,843
```

Hayles 2013, Methods 5.1:

> The 4843 deletion strain collection used for this analysis consists of 4825
> strains described by Kim et al. [7] plus 18 additional gene deletion mutants
> as shown in the electronic supplementary material 1, table S2.

`4,836 - 11 = 4,825`, so the 11 are exactly the difference between Kim's
constructed collection and the 4,825 Kim strains used by Hayles.

The loaders and reconciliation primitives live in the shared module
`src/dataset_merge.py`. The archived `01_format_and_update_ids.py` uses them
to build the production union table (`Presence`,
`Dispensability_consistency`, `Dispensability_discrepancy`,
`Dispensability_change_class`); the sections below record the interpretation
that is not encoded in those columns.

## 1. What each file contains

| File | Rows | What it is |
| --- | ---: | --- |
| `data/references/kim2010/kim2010_supp_table1.xls` | 4,914 | Kim 2010 full **Supplementary Table 1**; 4,836 rows have `Constructed? = yes`. |
| `data/references/hayles2013/hayles2013_supp_table3.xlsx` | 4,843 | Hayles 2013 main supplementary table. |

The 2010 table has two dispensability columns, defined in its cell comments:

- `Deletion dataset` — essentiality determined by the Kim project's own
  experiments (primary deletion screen).
- `Analysis dataset` — essentiality determined by both publication and the
  Kim experiments (the reconciled call).

The 2013 workbook also copies both columns (`Kim et al 2010 Deletion dataset`,
`Kim et al 2010 Analysis dataset`) and adds its own `Gene dispensability. This
study`. The copies match Kim's full table exactly for all 4,825 shared
constructed genes. The final union keeps Kim's originals only — the two 2013
copies are dropped — and the 2013 study call is written as
`Hayles2013_Gene_dispensability` (the live `01` table calls it
`Gene_dispensability`).

### Reconciliation roles

| Role | Count | Meaning |
| --- | ---: | --- |
| Shared, constructed | 4,825 | Kim-constructed and screened in 2013 |
| Shared, not constructed | 18 | Kim candidates added by Hayles (Table S2) |
| Only 2010, constructed | **11** | Constructed by Kim but absent from 2013 |
| Only 2010, not constructed | 60 | Kim candidates never constructed or screened |
| Only 2013 | 0 | Every 2013 gene is in Kim's 4,914-row candidate list |

The 71 Kim-only IDs are therefore not 71 "missing results": only the 11
constructed strains have a 2010 V/E call. The other 60 are design candidates
that were never constructed, so they never had a 2010 dispensability result.

At the **gene level**, the archived `01_format_and_update_ids.py` updates IDs before merging and folds
`SPAC27D7.10c` into `SPAC27D7.09c` (PomBase treats `.10c` as a synonym of
`.09c`). The union's `Presence == 2010_only` therefore has **10 rows**, not
11; the `.10c` historical record is kept in the `.09c` row's
`Original Systematic ID`. The 11-strain statement above remains true at the
strain level.

### Kim's own three-stage call

Kim's Table S1 reports three calls per constructed gene, and the archived
`01` keeps all of them in the Kim intermediate:

| Column | Meaning |
| --- | --- |
| `Kim2010_Deletion_dataset` | Kim's own knockout/primary-screen call |
| `Kim2010_Published_phenotype` | the literature call (`viable` / `essential` / `nodata`) |
| `Kim2010_Analysis_dataset` | Kim's final reconciled call |

Kim compared its own data with published data. In the current tables 1,209
genes have a V/E published call: **1,163 agree** with Kim's own deletion
call, **40 disagree**, and 6 have no Kim deletion call. The 40 disagreements
are listed in **Kim Supplementary Table 3**, which is only published inside
the supplementary PDF; it is hand-transcribed into
`data/references/kim2010/kim2010_supp_table3.tsv` (tracked in git) and
checked against Table S1 at runtime by `validate_kim_table_s3()`.

The stated criteria for the 40 disagreements:

| `criteria_for_judgment` | Genes | Meaning |
| --- | ---: | --- |
| `Disagreement` | 19 | Kim kept the published call in its final analysis |
| `More of CDS deleted in this study` | 13 | Kim used its own call |
| `Essential in our assay conditions` | 4 | Kim used its own call (condition-dependent) |
| `Different interpretation due to slow growth or dying cells` | 2 | Kim used its own call |
| `Published data not confirmed` | 2 | Kim used its own call |

The footnotes give the condition details for `atm1` (essentiality depends on
KO%), `vma1`/`vma3` (media pH), `met9` (methionine supplementation) and
`tol1` (sulfite supplementation). They are carried as
`Kim2010_Judgment_criteria_detail` in the Kim intermediate.

The final judgment is not simply "always use our data": for the 19
`Disagreement` genes Kim kept the published call in `Analysis dataset`, which
is exactly why some `Deletion dataset` and `Analysis dataset` values differ.

## 2. The 11 constructed strains with a 2010 result but no 2013 row

| ID | Gene (2010) | Product (2010) | KO% | 2010 deletion | 2010 analysis | Published phenotype | Family / paralog note | Likely category |
| --- | --- | --- | ---: | --- | --- | --- | --- | --- |
| `SPAC19B12.13` | `cox112` | Fusion Cox1102/Rsm2202 | 76 | V | V | nodata | Fusion gene; paralog `SPAC1420.04c` (`cox1101`) is in 2013 | gene fusion / paralog |
| `SPAC1F3.02c` | `mkh1` | MEK kinase (MEKK) Mkh1 | 97 | V | V | viable | No paralog | unknown |
| `SPAC23A1.10` | `ef1a-b` | Translation elongation factor EF-1α Ef1a-b | 100 | V | V | nodata | EF-1α duplicate; `ef1a-a` (`SPCC794.09c`) and `ef1a-c` (`SPBC839.15c`) are in 2013 | duplicated family |
| `SPAC25B8.12c` | — | Nucleotide-sugar phosphatase (predicted) | 96 | V | V | nodata | No paralog | unknown |
| `SPAC27D7.10c` | — | But2 family protein | 94 | V | V | nodata | Neighboring `.09c`/`.11c` are in 2013 | temperature-sensitive remake (`ts` flag) |
| `SPAC977.07c` | `pfl6` | Glycoprotein / flocculin Pfl6 | 100 | V | V | nodata | Four paralogs | paralogous family |
| `SPCC16A11.02` | `utp13` | U3 snoRNP-associated protein Utp13 | 79 | E | E | nodata | Essential; no paralog | unknown |
| `SPCC16C4.16c` | `cbc3` | Nuclear cap-binding complex variant subunit Cbc3 | 85 | V | V | nodata | No paralog | unknown |
| `SPCC736.05` | `wtf7` | wtf element Wtf7 | 48 | V | V | nodata | 13 paralogs | repetitive family |
| `SPCC737.04` | — | UPF0300 family protein 6 | 97 | V | V | nodata | Five paralogs | paralogous family |
| `SPCC970.11c` | `wtf9` | wtf element Wtf9 | 63 | V | V | nodata | Many paralogs | repetitive family |

In the gene-level union, `SPAC27D7.10c` is collapsed into `SPAC27D7.09c`, so
this table has 11 strain records but the union has 10 `2010_only` rows. The
collapsed row keeps both original IDs.

**The 2013 paper does not give a per-gene reason for excluding these 11.**
It states only that 4,825 Kim strains were used. The per-gene notes above are
annotations, not documented exclusions:

- One strain (`SPAC27D7.10c`) is explicitly marked in Kim's own table as a
  heterozygous diploid containing a `ts` mutation "being remade".
- Six are members of large duplicated/repetitive families (`cox112`,
  `ef1a-b`, `pfl6`, `wtf7`, `wtf9`, `SPCC737.04`); for each, other family
  members are present in 2013. This is consistent with strain-validation
  difficulty, but it is not stated in the paper.
- `mkh1`, `utp13`, `SPAC25B8.12c` and `SPCC16C4.16c` have no obvious
  biological reason in the available records. For these, "strain not
  available / not included in the screened library" is the most likely
  explanation, but it is not documented.

If the exact reason matters, the next place to look is the strain-level
library QC record (the 2010 supplementary validation tables `bPCR`, `sPCR` and
`GS`), not the 2013 paper. Absence from those validation sheets does not
cleanly separate the 11: 1,606 of the 4,836 constructed genes are absent from
all three sheets, so validation-sheet absence alone is not diagnostic.

## 3. The 60 not-constructed candidates absent from 2013

Kim's 4,914-row table includes 78 genes with `Constructed? = no`. Of these,
18 were constructed/added by Hayles (next section) and 60 are absent from
2013. The 60 have neither a `Deletion dataset` nor an `Analysis dataset`
call, so they are **not** "2010 experimental results missing from 2013"; they
are candidate genes that were never made into strains. Some carry a
`Published phenotype` annotation from the literature, but none has a Kim
knockout V/E call. In the archived 01 union they are the rows with
`Presence` = `neither`.

## 4. The 18 Hayles additions (Table S2)

Hayles Table S2 is titled "Dispensability of 18 genes not previously included
in the set of genome-wide gene deletions". All 18 are Kim `Constructed? = no`
candidates. The table has two groups.

### 4.1 Six newly constructed deletions

| ID | Gene | 2013 call | % deleted |
| --- | --- | --- | ---: |
| `SPCC320.04c` | — | V | 99 |
| `SPCC4B3.16` | `tip41` | V | 83 |
| `SPCPJ732.03` | `meu15` | V | 100 |
| `SPAC22E12.16c` | `pik1` | E | 94 |
| `SPAC11E3.12` | — | V | 100 |
| `SPAC750.08c` | — | V | 100 |

### 4.2 Twelve genes not analysed by Kim because they were <40% deleted

| ID | Gene | 2013 call | % deleted |
| --- | --- | --- | ---: |
| `SPAC26A3.08` | `smb1` | E | 6 |
| `SPAC3G9.03` | `rpl231` | V | 9 |
| `SPAC824.06` | `tim14` | E | 5 |
| `SPAC977.06` | — | V | 8 |
| `SPBC1347.01c` | `rev1` | V | 2 |
| `SPBC16E9.01c` | `php4` | V | 4 |
| `SPBC17A3.02` | — | V | 34 |
| `SPBC405.01` | `ade1` | V | 3 |
| `SPCC1235.01` | — | V | 4 |
| `SPCC1682.16` | `rpt4` | V | 5 |
| `SPCC24B10.06` | — | V | 32 |
| `SPAC9E9.04` | — | V | 40 |

The last two (`SPBC17A3.02`, `SPCC24B10.06`) carry a `V` in Kim's
`Deletion dataset` column (and in the 2013 copy of it) but no analysis-dataset
value. The
other ten have neither marker. This is why the 2013 workbook has 16 blank
deletion markers but 18 blank analysis markers. Their 2010 V/E calls are
missing by design, so they must not be counted as "2013 has no 2010 result".

## 5. What the 2013 paper actually corrected

The 2013 supplement has two relevant tables:

- **Table S2** — the 18 additions above.
- **Table S3** — "Gene dispensability changes from Kim et al 2010": 21 genes,
  of which 9 had the deletion remade and 12 were re-analysed.

Table S3's closing note:

> 3 gene deletion strains have been reconstructed as they had a different
> phenotype to previously published data. We find that 2 are still different
> to published data and one is now the same as previously published. Deletion
> mutants of 6 other genes that were non-essential (V) were remade and these
> genes are now essential (E). 12 gene deletion mutants were re-analysed and
> found to have a different dispensability, 7 of these now have the same
> phenotype as previously published deletion mutants.

Table S3 is carried into the union as three per-gene columns (blank outside
its 21 genes): `Hayles2013_Strain_status` (`Deletion remade` for 9,
`Re-analysed` for 12), `Hayles2013_2010_dispensability_change` (the
Kim-2010-analysis → 2013-study change, e.g. `V -> E`, `E -> V`, `same (E)`),
and `Hayles2013_TableS3_comment` (the verbatim supplement comment, including
`same as pub` / `confirmed different to pub`).

### 5.1 Change taxonomy

For the 4,825 shared constructed genes, the 2010 `Analysis dataset` and the
2013 `This study` call agree for 4,796 and differ for 29
(12 E→V, 17 V→E). The 2013 table also contains the 18 additions, so 54 rows
have an incomplete or changed record:

| Class | Count | Meaning |
| --- | ---: | --- |
| `hayles2013_addition` | 18 | No 2010 call; Table S2 additions. |
| `documented_change` | 14 | 2013 study differs from Kim's **analysis** call and is listed in Table S3. |
| `documented_discrepancy_in_kim2010_columns` | 7 | 2013 study agrees with Kim's **analysis** call but differs from Kim's **deletion** call; listed in Table S3 because Table S3 shows both columns. |
| `matches_kim2010_deletion_not_analysis` | 11 | 2013 study differs from Kim's analysis call but equals Kim's deletion call. Not a new correction; the 2013 screen sided with Kim's primary call. |
| `undocumented_change` | 4 | 2013 study differs from **both** Kim columns and is **not** in Table S3. |

So Table S3's 21 rows decompose as 14 changes against Kim's final (analysis)
call plus 7 rows where the 2013 call simply confirms Kim's analysis call
against Kim's primary deletion call.

Kim's own `Deletion dataset` vs `Analysis dataset` columns disagree for 21
rows (19 true V↔E revisions + 2 partial calls); in the union all 21 are
folded into the documented classes above (14 `documented_change` + 7
`documented_discrepancy_in_kim2010_columns`). The generic
`kim2010_deletion_analysis_discrepancy` class is never emitted on the default
run; it appears only when the Hayles docx is unavailable. These are the source
of the 11 "reversion" rows below.

### 5.2 The 4 undocumented V/E changes

These four are present in the 2013 main table and in the 2013 phenotype
category sheets, but absent from Table S3:

| ID | Gene | Kim deletion | Kim analysis | 2013 study | 2013 phenotype | Note |
| --- | --- | --- | --- | --- | --- | --- |
| `SPAC17C9.13c` | `cut8` | E | E | V | long high penetrance | 2013 description says "VIABLE slightly misshapen cells at 25, ESSENTIAL germinated spores long, often divide at 32"; the published phenotype and `cut8Δ` phenotypes are essential/inviable. The 2013 V call is the outlier. |
| `SPAC25B8.06c` | `dia4` | V | V | E | skittle | Mitochondrial serine-tRNA ligase; 2013 category sheet says E. Likely genuine correction, but undocumented. |
| `SPAC2F3.06c` | `kap104` (2010: `kap14`) | V | V | E | misshapen essential | Karyopherin/importin β; 2010 published phenotype was viable; 2013 category sheet says E. Likely genuine correction, but undocumented. |
| `SPAC823.12` | `vps11` | V | V | E | misshapen essential | 2013 category sheet says E. Likely genuine correction, but undocumented. |

PomBase 2026 annotates all four as **"Gene Deletion Viability: Depends on
conditions"**, which is why the 2010 and 2013 calls can differ by condition.
The current phenotype evidence mostly supports the 2013 call:

- `kap104Δ`: "inviable after spore germination, multiple cell divisions,
  abnormal morphology"; "inviable vegetative cell population"; loss of
  viability at high and low temperature.
- `vps11Δ`: "inviable after spore germination, multiple cell divisions,
  abnormal morphology".
- `dia4Δ`: "inviable vegetative cell population"; "inviable tapered
  vegetative cell".
- `cut8Δ`: "inviable after spore germination, single or double cell division"
  and "inviable after spore germination, without cell division, with elongated
  germ tube"; the literature describes Cut8 as "essential for anaphase"
  (Tatebe et al., Curr Biol 2000, PMID 11084332). **This one supports the 2010
  E call, not the 2013 V call.**

So three of the four are likely genuine corrections missing from Table S3
(`kap104`, `vps11`, `dia4`), while `cut8` is more likely an error in the 2013
main table or a condition-specific call. Either way, the four are internally
consistent within the 2013 main table and category sheets; the inconsistency
is between the 2013 main table and Table S3. The paper does not resolve it.

PomBase pages: <https://www.pombase.org/gene/SPAC2F3.06c>,
<https://www.pombase.org/gene/SPAC823.12>,
<https://www.pombase.org/gene/SPAC25B8.06c>,
<https://www.pombase.org/gene/SPAC17C9.13c>.

### 5.3 The 11 genes where 2013 matches Kim's deletion call, not the analysis call

`rho1`, `hus5`, `spo20`, `bot1`, `msp1`, `rpt6`, `ste11`, `tlg1`, `sec31`,
`ssl1`, `bgs4`.

For these, the 2013 visual screen agrees with Kim's primary
`Deletion dataset` call rather than Kim's `Analysis dataset` call. Examples:
`rho1` (deletion V, analysis E, 2013 V, WT), `hus5` (V/E/V), `ste11`
(E/V/E). These are not new corrections to Kim's final call; they are
reversions to the primary call. They are worth flagging when using the 2010
analysis dataset as a gold standard.

### 5.4 Do strains outside Table S3 reproduce Kim?

Only the 21 Table S3 strains were remade (9) or re-analysed (12). The other
4,822 Hayles rows comprise 4,804 genes with a Kim call and the 18 Table S2
additions. Their 2013 calls break down as:

| `Hayles2013_vs_Kim2010` | Rows |
| --- | ---: |
| `same as Kim analysis` | 4,789 |
| `matches Kim deletion (not analysis)` | 11 |
| `differs from both Kim calls` | 4 |
| `no Kim call` | 18 |

So **15 of 4,804** strains with a Kim call do not reproduce Kim's final
`Analysis dataset` call, and **4 of 4,804** do not reproduce Kim's own
`Deletion dataset` call either. The answer is therefore **no, not exactly**
— but the exceptions are small and mostly inherited from Kim's own
literature reconciliation.

The 11 "matches deletion" rows are not new Hayles disagreements: Kim's own
Table S1 already had `Deletion dataset != Analysis dataset` for these,
because Kim had used the published call in its final analysis. The 2013
screen agrees with Kim's primary experiment. The 4 genuine differences are
the undocumented changes in §5.2 (`cut8`, `dia4`, `kap104`, `vps11`).

At the whole-table level, including Table S3, the comparison column in the
Hayles intermediate reads:

| `Hayles2013_vs_Kim2010` | Rows |
| --- | ---: |
| `same as Kim analysis` | 4,796 |
| `matches Kim deletion (not analysis)` | 12 |
| `differs from both Kim calls` | 17 |
| `no Kim call` | 18 |

The whole-table "differs from both" count adds the 13 Table S3 changes to the
4 undocumented changes; the twelfth "matches deletion" row is `cti1`
(`SPCC1739.07`), a remade strain that returned to Kim's deletion call.

## 6. Data-quality issues that affect merges

The most important pitfalls, all handled by the loader or recorded here:

1. **Call encoding.** Use `V`/`E` consistently; the loader also accepts
   `viable`/`essential` defensively. Strip stray whitespace before comparing.
2. **Constructed status.** 78 of 4,914 Kim rows are `Constructed? = no`; they
   have no full analysis call. Do not treat them as missing 2013 data. The 18
   that appear in 2013 are the Table S2 additions.
3. **Whitespace.** The 2013 table has stray whitespace in 9
   `Gene dispensability. This study` values (renamed
   `Hayles2013_Gene_dispensability` in the output), 58 phenotype classifications,
   170 phenotype descriptions, 4 gene-product descriptions, 1 paralogue
   field, and 1 cell-length value. Kim's table also has a few trailing spaces
   in primer/gene-function fields.
4. **Excel date corruption in gene names.** Three 2013 `Gene name` cells are
   stored as dates (`mmm-yy` format). `load_hayles2013()` reconstructs them:
   - `SPAC1F3.10c` → `2001-10-01` → `Oct-1` (`oct1`)
   - `SPAC5D6.05` → `2011-09-01` → `Sep-11` (`med18`, synonyms
     `pmc6`, `sep11`)
   - `SPBC4C3.12` → `2001-09-01` → `Sep-1` (`sep1`)
5. **Gene-name drift.** 294 shared genes have different gene names between
   the two tables (e.g. `kap14` → `kap104`, `rps152` → `rps1502`,
   `mug31` → `pom34`, `rad5` → `rad50`, `spo2` → `spo20`). Match on
   systematic IDs only.
6. **Blank Kim-2010 markers in 2013.** 16 rows have no deletion marker and 18
   have no analysis marker; they are the Table S2 additions.
7. **Cell comments.** The 2013 workbook carries the definitions of the two
   Kim columns (E1/F1) and a strain-mix-up note at C3885: "w3582-1(b3669,
   SPBC29A3.21) and w3582 (SPBP23A10.16) are completely different genes."
   `pandas.read_excel` drops comments silently; they were read with
   `openpyxl` during the audit. The 2010 `.xls` comments (D1/E1/L1/S1) define
   `Deletion dataset`, `Analysis dataset`, the construction methods and the
   abbreviation key; `xlrd` cannot read them, so they are recorded as a known
   limitation.
8. **Truncated ID in Table S3.** The `alp11` row is written as `SPAC13D6.0`;
   the full ID is `SPAC13D6.05`. `parse_table_s3()` repairs it explicitly.
9. **Allele-style ID.** `SPAC1556.06.1` (`meu1`/`meu1-1`) does not match the
   standard `SP[A-Z0-9]+.NN[c]` pattern. Keep it; do not drop it as malformed.
10. **Header typo.** The 2010 locus column is `Systemic ID`, not
    `Systematic ID`.
11. **Stray value.** The 2013 workbook has an unnamed trailing column with a
    single `x`. Ignore the column.
12. **Non-ASCII.** One phenotype description contains `36°C`; write UTF-8.
13. **Paralog columns.** The 2013 `Paralogues` column is essentially a copy
    of Kim's `Paralog` (same pipe-separated IDs, order, and even duplicate
    tokens), with the sentinel `none` → `NONE` and only three content
    differences: `.10c` is removed from the `SPAC27D7.09c`/`.11c` lists
    (annotation update), and `SPAC212.01c` carries the Kim paralog list of
    the dropped `SPCC737.04` row (a likely copy-paste error). Both `01` and
    the archived `01` therefore use Kim's `Paralog`, update its IDs/names,
    and drop `Paralogues`.

## 7. Reproducing the reconciliation

The machine-readable reconciliation is the union table written by the
archived `01`:

```python
import pandas as pd

union = pd.read_excel(
    "tmp/01_format_and_update_ids/data/1_formatted/"
    "Kim2010_Hayles2013_merged_formatted.xlsx"
)
union[union["Presence"] == "2010_only"]        # 10 gene-level rows; .10c is folded into .09c
union[union["Presence"] == "neither"]          # the 60 never-constructed candidates
union[union["Presence"] == "2013_only"]        # the 18 Table S2 additions
union[union["Dispensability_consistency"] == "inconsistent"]   # the 29 call changes
union[union["Dispensability_change_class"] == "undocumented_change"]  # the 4
```

The former `00b` script's 13 generated files (12 TSVs and a summary) were
removed with it. Their content is reconstructible from the union with the
filters above, except for the 316-row data-quality ledger; its important
counts and pitfalls are summarised in section 6.

Reading the Kim `.xls` requires `xlrd` (already installed in the
`bioinformatics` environment).

### Optional reference file

`data/references/` is gitignored. For the Table S3 `documented_change` /
`documented_discrepancy_in_kim2010_columns` classes, the archived `01` reads the Hayles
Supplementary Material 1 (`hayles2013/hayles2013_supp.docx`); without it 17
rows degrade to `changed_table_s3_not_available`, `cti1` joins
`matches_kim2010_deletion_not_analysis`, and the 7 documented discrepancies
become `kim2010_deletion_analysis_discrepancy`. It is available from the PMC
supplementary archive:

<https://www.ebi.ac.uk/europepmc/webservices/rest/PMC3866870/supplementaryFiles>

## 8. References

- Kim DU, Hayles J, Kim D, Wood V, Park HO, Won M, Yoo HS, et al.
  *Analysis of a genome-wide set of gene deletions in the fission yeast
  Schizosaccharomyces pombe.* Nature Biotechnology 28:617–623 (2010).
  DOI: [10.1038/nbt.1628](https://doi.org/10.1038/nbt.1628).
  PMID: 20473289; PMC3962850.
- Hayles J, Wood V, Jeffery D, Hoe KL, Kim DU, Park HO, Salas-Pino S,
  Heichinger C, Nurse P. *A genome-wide resource of cell cycle and cell shape
  genes of fission yeast.* Open Biology 3:130053 (2013).
  DOI: [10.1098/rsob.130053](https://doi.org/10.1098/rsob.130053).
  PMID: 23697806; PMC3866870.
- PomBase release 2026-06-01, used for current gene names/products and
  phenotype evidence (`data/references/pombase-2026-06-01_gene_IDs_names_products.tsv`).
