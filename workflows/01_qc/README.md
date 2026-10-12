# 01_qc: Cell quality control

Removes low-quality cells from a Seurat object. It runs in two steps with a
human decision in between: first look at the QC metrics, then choose the
cutoffs, then apply them.

## Start here

**You need:** one Seurat object saved as `.rds` in `data/raw/`.

**You edit only two files**, both in `config/`:

| File | What goes in it | When |
|---|---|---|
| `params_override.R` | What you already know about your dataset (grouping columns, cells to exclude up front, species gene patterns) | Before step 1 |
| `decisions.R` | The QC cutoffs you choose after looking at step 1's plots | Between step 1 and step 2 |

Create them by copying the examples:
`params_override_example.R` → `params_override.R`, `decisions_example.R` → `decisions.R`.

> **Do not edit files in `pipeline/` or `shared/` to analyse a new dataset.**
> Dataset settings belong in `params_override.R`; your choices belong in `decisions.R`.

## Run

Open `scRNAseq-tools.Rproj` in RStudio first (so paths resolve), then:

1. **Characterize QC metrics**
   ```r
   source(here::here("workflows", "01_qc", "pipeline", "01_qc_prefilter.R"))
   ```
   Look at the plots and tables in `results/01_qc/01_qc_prefilter/`
   (violin plots, scatter plot, summary tables, MAD-suggested cutoffs).

2. **Record your cutoffs** in `config/decisions.R`, for example:
   ```r
   params$thresholds <- list(
     nFeature_RNA = list(min = 600, max = 3500),
     nCount_RNA   = list(max = 1000000),
     percent.mt   = list(max = 50)
   )
   ```
   Cells are kept when `value > min` and `value < max`.

3. **Apply the cutoffs and assess the result**
   ```r
   source(here::here("workflows", "01_qc", "pipeline", "02_qc_postfilter.R"))
   ```

If a setting name in `params_override.R` is misspelled or not a valid
setting, the scripts stop and list it, instead of silently ignoring it.
Valid names are the ones in `config/config.R`. Set single entries
(`params$retained_vs_excluded$padj_cutoff <- 0.01`) rather than replacing a
whole list.

## Optional sections

Switched on/off in `params_override.R` (defaults in `config/config.R`):

| Setting | What it adds |
|---|---|
| `steps$mad_estimation` | Suggested cutoffs from the median absolute deviation (advisory only) |
| `steps$retained_vs_excluded` | Compares excluded with retained cells (DE, GO/KEGG, and an Excel workbook with every table) to check that the cutoffs removed low-quality cells rather than a distinct population. Uses human gene annotation: switch off for other species. KEGG needs an internet connection |
| `steps$spatial_qc` | QC metrics and retained/excluded cells by plate/well position (plate-based protocols only) |

## Output

| Where | What |
|---|---|
| `results/01_qc/01_qc_prefilter/` | QC metrics before filtering; cell counts; MAD-suggested cutoffs |
| `results/01_qc/02_qc_postfilter/` | QC metrics after filtering; cells retained per group; `qc_de_enrichment.xlsx` and the excluded vs retained plots; spatial QC |
| `data/processed/01_qc_filtered.rds` | QC-filtered object: **the input to `02_population_refinement`** |
| `data/processed/*_manifest.txt` | Next to each checkpoint: the settings used, its input and the input's checksum, cell count, date, repo version and checksum |

## Record what you did

Keep a copy of `config/decisions.R` and `config/params_override.R`, and note
the version in the repo's `VERSION` file. Together they reproduce your run.
