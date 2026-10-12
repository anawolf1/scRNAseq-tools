# scRNAseq-tools

Modular R workflows for single-cell RNA-seq of sorted cells: quality
control, clustering and population refinement, functional program scores,
and differential abundance / expression analysis. Built with Seurat, miloR
and miloDE.

Developed and tested on plate-based data (FLASH-seq) from human samples.
Droplet-based data would need additional preprocessing before `01_qc`.
Non-human data needs the species settings changed: the gene patterns in `01`
and `02` (`percent_features`, `clustering$exclude_gene_patterns`), and the
GO/KEGG steps switched off (`01`: `retained_vs_excluded`, `02`:
`exclusion_check`); `03` maps gene sets to human symbols.

## Structure

```
scRNAseq-tools/
├── data/
│   ├── raw/          your Seurat object (.rds) - never committed
│   ├── processed/    checkpoints passed between workflows - never committed
│   └── gene_sets/    curated gene sets (.gmt) for functional scores
├── shared/functions/ functions used by more than one workflow
├── workflows/
│   ├── 01_qc/
│   ├── 02_population_refinement/
│   ├── 03_functional_scores/
│   └── 04_neighborhood_analysis/
└── results/          figures and tables - created when you run, never committed
```

Every workflow has the same layout:

| Folder / file | Contents |
|---|---|
| `README.md` | What the workflow does and how to run it - **start here** |
| `setup.R` | Loads functions, settings and packages (sourced by the scripts) |
| `config/config.R` | Every setting, with generic defaults |
| `config/params_override_example.R` | Template for your dataset's settings |
| `config/decisions_example.R` | Template for the choices you make after looking at results |
| `pipeline/` | The scripts you run, in numbered order |
| `functions/`, `tools/` | Workflow-specific functions; helper scripts for choosing parameters |

**To analyze a dataset you only edit two files per workflow:** copy
`params_override_example.R` → `params_override.R` and `decisions_example.R` →
`decisions.R`, and fill them in. Both names are gitignored, so your data-specific
values stay local. A misspelled setting stops the scripts instead of being ignored.

## Workflows

| Workflow | Does | Reads | Writes |
|---|---|---|---|
| `01_qc` | Metadata exclusions, QC metrics, cutoffs, checks of what was removed | `data/raw/*.rds` | `01_qc_filtered.rds` |
| `02_population_refinement` | Clustering in rounds, characterization, removal of populations that do not belong | `01_qc_filtered.rds` | `02_population_refined.rds` |
| `03_functional_scores` | Module scores for curated programs; cluster and timepoint comparisons on patient means | `02_population_refined.rds` | `03_functional_scores.rds` |
| `04_neighborhood_analysis` | Differential abundance (miloR) and neighborhood DE (miloDE) between timepoints | `02_population_refined.rds` (+ `03_functional_scores.rds`) | `04_milo_da.rds`, `04_milo_de_*.rds` |

Workflows communicate only through files in `data/processed/`, so each can be
run on its own:

- **As a chain.** Run in order; each workflow picks up the previous one's
  checkpoint by default. `03` and `04` both start from `02_population_refined.rds`
  and do not depend on each other: without `03_functional_scores.rds`, `04`
  runs but skips the program-score plots.
- **From your own object.** Every workflow's input is a setting (`01`:
  `params$paths$raw_data_file`; `02`: round 1's `input` in `decisions.R`;
  `03`, `04`: `params$input`), which accepts a file path as well as a
  checkpoint name. An already-clustered Seurat object can go straight to `03`
  or `04`, provided it has what that workflow needs (see *Input requirements*
  below).

Re-running a workflow does not update later checkpoints or results. Nothing
warns you at that point; instead, a workflow stops when it would combine a
checkpoint with a newer version of the input it was made from (e.g. `04`
with a `03_functional_scores.rds` made before `02` was re-run). Re-run the
workflows downstream of the one you changed.

Each checkpoint has a `_manifest.txt` next to it: the script and settings that
produced it, its input and the input's checksum, date, repository version and
checksum, and for Seurat objects the cell and gene counts, assays, reductions
and metadata columns.

### Input requirements

What a workflow needs when you give it your own object instead of the
previous workflow's checkpoint. `03` and `04` check their required columns
and reductions and stop with a message if one is missing.

| Workflow | Required | Optional |
|---|---|---|
| `01_qc` | Raw counts (RNA assay) | Columns named in `params$group_by_vars`, `cell_count_vars` or `exclusions` (all empty by default); plate position columns (only with `steps$spatial_qc`, off by default) |
| `02_population_refinement` | Counts | QC metric columns from `01` (`steps$qc_metrics`, `clustering$vars_to_regress`); columns named in `group_by_vars`, `cell_count_vars` or a round's exclusion steps; protein columns for `marker_check_cols` |
| `03_functional_scores` | Patient and timepoint columns (`params$study_design`); `cluster_label`; gene sets in `data/gene_sets/`, with the species of each set in `params$gene_sets$species` | Columns used by the cell sets in `decisions.R` |
| `04_neighborhood_analysis` | Patient and timepoint columns; `cluster_label`; `pca` and `umap` reductions; the timepoints named in the contrasts | `03_functional_scores.rds` (program scores per neighborhood); QC metric columns (per-neighborhood QC, used if present) |

## Requirements

- R 4.5.1 (the version recorded in `renv.lock`; other versions are untested)
- RStudio recommended (paths resolve from `scRNAseq-tools.Rproj`)

## Getting started

1. Clone the repository and open `scRNAseq-tools.Rproj` in RStudio (all paths
   resolve from the project root).
2. Install the pinned package versions:
   ```r
   renv::restore()
   ```
3. Put your Seurat object in `data/raw/` (one `.rds`). If you will run `03`
   or `04`, its metadata needs a patient and a timepoint column (`01` and
   `02` carry them through but do not need them). To start from a later
   workflow instead, see *Input requirements*.
4. Follow `workflows/01_qc/README.md`, or the README of the workflow you start
   from.

`library_help()` lists every package the repository uses and which workflow
needs it; `palette_help()` lists the color palettes and how to change them.

## Reproducibility

- Package versions are pinned in `renv.lock`.
- Every choice made after looking at results (cutoffs, clustering parameters,
  exclusions, labels, contrasts) is recorded in a workflow's `decisions.R`,
  the audit trail of the analysis. `decisions.R` is gitignored, so keep your
  own copy under version control (e.g. a private repository) or backed up.
- Workflows communicate only through checkpoints, and each checkpoint records
  the checksum of the input it was made from. A workflow stops when a
  checkpoint it reads was made from an earlier version of its input;
  population refinement also stops when a round's settings changed since it
  was computed.
- To add a package: add it to `library_registry` in `shared/functions/libraries.R`,
  run `tools/write_description.R`, then `renv::snapshot()`.
- The repository version is in `VERSION`.

## Citation

If you use these workflows, please cite the methods they are built on:

- Seurat: Hao et al. (2024), *Nature Biotechnology* 42, 293–304
- miloR: Dann et al. (2022), *Nature Biotechnology* 40, 245–253
- miloDE: Missarova et al. (2024), *Genome Biology* 25, 189

Run `citation("<package>")` in R for the package references.

## License

MIT, see [`LICENSE`](LICENSE).
