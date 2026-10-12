# 02_population_refinement: Clustering and population refinement

Clusters the QC-filtered cells, characterizes each clustering, and removes
populations that do not belong (e.g. low-quality cells or contaminating lineages), in **rounds**.
The final round is the clustering later workflows use.

## Start here

**You need:** the QC workflow's output, `data/processed/01_qc_filtered.rds`.

**You edit only two files**, both in `config/` (copy them from the `*_example.R` files):

| File | What goes in it |
|---|---|
| `params_override.R` | What you know about your dataset: grouping columns, variables to regress, genes to exclude, signature genes, colors |
| `decisions.R` | Your rounds: clustering parameters, exclusions, cluster labels, and which round is final |

> **Do not edit files in `pipeline/`, `functions/` or `shared/` to analyse a new dataset.**

Open `scRNAseq-tools.Rproj` in RStudio before running anything, so paths resolve.

| Script | What it does |
|---|---|
| `pipeline/01_population_refinement.R` | Runs every round in `decisions.R` |
| `tools/select_clustering_parameters.R` | Screens dims × resolution (× k_param × prune_snn) for a round |
| `tools/inspect_clustering_candidates.R` | Compares a few candidate clusterings in depth; saves cell lists |
| `tools/compare_rounds.R` | Compares rounds run as alternatives |

---

## How rounds work

A **round** is one pass of:

```
load input ─▶ steps (remove cells/genes) ─▶ preprocess ─▶ cluster ─▶ UMAP ─▶ characterize ─▶ save
```

Each round is an entry in `config/decisions.R`. Four fields are required:
`input`, `algorithm`, `dims`, `resolution`. The others are optional:
`label`, `rationale`, `k_param`, `prune_snn` (defaults in `config.R`),
`exclude_genes`, `batch_correction`, `steps`, `cluster_labels`,
`pairwise_de`. Fields ending in `sweep` record how parameters were chosen and
are not read. Any other field stops the pipeline, so a typo is not ignored.

**Rounds chain through `input`.** Round 1 reads `"qc_filtered"`; a later
round reads an earlier round's result by its name (`input = "round_1"`).
Every round is preprocessed and clustered from scratch: only the cells and
metadata of the input are used, not its PCA or clusters.

**A round finds populations; the next round removes them.** A round's steps
run on its input, before clustering. So round 1 clusters and characterizes
all cells; you look at it and label its clusters; round 2 removes what you
decided (e.g. the cells round 1 labeled `"Low quality"`) and clusters the
rest. Each round's figures show exactly the cells it clustered, and every
round's checkpoint keeps all the cells it clustered, so alternative rounds
(remove all, half, a third) can all start from `input = "round_1"`.

**Rounds are a record.** You append rounds; you do not rewrite past ones.
Together with git history, `decisions.R` is the audit trail of your choices.

### Removing cells: `steps`

A step removes cells (or genes). Three kinds:

| Step | Removes |
|---|---|
| `list(type = "exclude_metadata", column = "cluster_label", values = "Low quality")` | cells whose metadata column has one of these values (any column: cluster, patient, plate…) |
| `list(type = "exclude_cells", cells_file = "low_quality_half.txt")` | specific cells, from a file saved with `save_cell_list()` (or `cell_ids = c(...)` inline) |
| `list(type = "exclude_genes", genes = c("XIST"))` | genes (`pattern = TRUE` for one regular expression) |

Steps select on the input's metadata: the previous round's `cluster_label`
or `seurat_clusters`, or any other column. A step field the step type does not
use stops the round, so a typo is not ignored. To remove the
`params$clustering$exclude_gene_patterns` genes, set the round's
`exclude_genes = TRUE` instead of writing a gene step.

**`expect_n_cells`**: add it to cluster-based steps. Cluster numbers only stay
the same while every input and setting does; if a step would remove a different
number of cells than when you decided on it, the round stops instead of
removing other cells.

**Titles.** Every figure and table of a round is titled with the round's
`label`, and its files start with the round's name. Write the label as you
want it in the figures (e.g. `label = "Final clustering"`).

### Re-running

- A round that is already computed and unchanged is **skipped**.
- If a round's settings (including the clustering and preprocessing defaults
  in `config.R`/`params_override.R`), **its input file** or a `cells_file`
  changed (e.g. QC re-run with new cutoffs), the script **stops** and names
  what changed, because reusing it would mix old and new results. Set
  `params$population_refinement$force_rerun <- TRUE` in `params_override.R`
  to recompute (every later round then recomputes too).
- Changing a round's `label`, `rationale` or `cluster_labels` does not change
  that round: labels are applied from `decisions.R` whenever they are used (the
  final checkpoint, rounds and tools that read the round). But rounds that
  build on it stop, because their steps may select on those labels.

### A typical history

```r
round_1 = list(label = "Lineage screening", input = "qc_filtered", algorithm = 1,
               dims = 10, resolution = 0.6, steps = list(),
               cluster_labels = c("0" = "Memory", "1" = "Effector", "2" = "Low quality", "3" = "Naive")),

round_2 = list(label = "Final clustering", input = "round_1", algorithm = 4,
               dims = 12, resolution = 0.5, exclude_genes = TRUE,
               steps = list(list(type = "exclude_metadata", column = "cluster_label",
                                 values = "Low quality", expect_n_cells = 58)),
               cluster_labels = c("1" = "Memory", "2" = "Effector", "3" = "Naive"))

params$population_refinement$final_round <- "round_2"
```

1. Screen parameters for round 1 (`select_clustering_parameters.R` on `"qc_filtered"`), record the round, run the pipeline.
2. Look at round 1's results; fill in `cluster_labels`.
3. Write round 2 with its step, screen its parameters with
   `sweep_for_round <- "round_2"` (the tool then removes the low-quality cluster first,
   as the round will), fill in `dims`/`resolution`, run the pipeline.
4. Label round 2's clusters and set `final_round`.

### Comparing alternative exclusions

Say round 1 has a low-quality cluster and you want to compare removing **all**
of it, **half**, or **a third**. Each alternative is simply a round from the same
input:

1. **Find the parts.** Subcluster it: run
   `select_clustering_parameters.R` and then `inspect_clustering_candidates.R`
   with `sweep_checkpoint <- "round_1"` (`inspect_checkpoint` in the inspect
   tool), `subset_column <- "cluster_label"`, `subset_values <- "Low quality"`.
   Pick the subclusters that make up "half" and "a third", and save them from
   the inspect tool:
   ```r
   save_cell_list(Seurat::WhichCells(inspection_results[["dims8_res0.4"]]$obj, idents = "1"),
                  "low_quality_half", data_processed_dir)
   ```
2. **Write one round per alternative**, all with `input = "round_1"`:
   ```r
   round_2_full  = list(label = "Remove all of the low-quality cluster", input = "round_1", ...,
                        steps = list(list(type = "exclude_metadata", column = "cluster_label",
                                          values = "Low quality"))),
   round_2_half  = list(label = "Remove half of the low-quality cluster", input = "round_1", ...,
                        steps = list(list(type = "exclude_cells", cells_file = "low_quality_half.txt"))),
   round_2_third = list(label = "Remove a third of the low-quality cluster", input = "round_1", ...,
                        steps = list(list(type = "exclude_cells", cells_file = "low_quality_third.txt")))
   ```
3. **Use the same parameters for every alternative**, so that only the removed
   cells differ between them. Check them for each alternative: run
   `select_clustering_parameters.R` with `sweep_for_round <- "round_2_half"`, and
   so on, and confirm the chosen dims/resolution are still in a stable region.
   If they are not for one alternative, report that rather than changing its
   parameters.
4. **Run the pipeline.** Each alternative gets the full characterization, DEG
   Venns and its own **exclusion check** (excluded vs retained cells, cells
   retained per group).
5. **Compare them** with `tools/compare_rounds.R`: summary table, how clusters
   correspond between alternatives, and per population, whether its DEGs depend
   on what was removed.
6. **Keep one**: note why in its `rationale`, set it as `final_round` (or continue
   from it with a round 3). The other alternatives stay in `decisions.R` as the
   record of what you compared.

---

## Optional sections

Switched off in `params_override.R` with `steps$<name> <- FALSE` (all on by default, see `config/config.R`):

| Setting | What it adds to each round |
|---|---|
| `composition_stability` | Cluster sizes, percentages and per-cluster stability (Jaccard, under subsampling) |
| `qc_metrics` | QC metrics by cluster |
| `group_umaps` | UMAPs colored by each grouping column |
| `cell_counts` | Cell counts per `cell_count_vars` entry, and cells per cluster and grouping column |
| `marker_characterization` | Marker-column UMAPs and tables (`parameter_selection$marker_check_cols`), signature-gene violins |
| `differential_expression` | One-vs-rest DE per cluster: volcanoes, tables, `_de_results.csv` and a workbook |
| `deg_overlap` | Venns of genes up (and separately down) in 2-4 clusters vs the rest, so shared genes change the same way |
| `exclusion_check` | For rounds that remove cells: excluded vs retained DE (table, volcano), GO/KEGG of the genes up in each group, a workbook with every result, and cells retained per group. GO/KEGG use human annotation; KEGG needs an internet connection |

The exclusion check uses `run_de_enrichment()` (shared): the enrichment
background is the genes detected in at least `exclusion_check$min_pct` of
either group. Its p-values come from a cell-level Wilcoxon test, which
treats cells as independent: they describe the cells compared, not
replicate-level evidence.

## DE results per round

Thresholds are in `params$de` (`min_pct`, `log2fc_cutoff`, `padj_cutoff`);
`params$de$exclude_gene_patterns` removes genes from every DE result.

- **One-vs-rest** (`differential_expression`): `<round>_de_one_vs_rest.xlsx`, one sheet per
  cluster, plus volcanoes and top-gene tables.
- **Pairwise**, if a round has `pairwise_de` in `decisions.R`: DE between every pair of
  clusters, with volcanoes and `pairwise_de/<round>_pairwise_de.xlsx` (a summary sheet with
  DEGs per pair, then one sheet per pair). Leave clusters out by label:
  ```r
  pairwise_de = list(exclude = c("Cytotoxic"))   # or pairwise_de = TRUE for all clusters
  ```
  It runs from the saved checkpoint, so it can be added to a round that is already
  computed (delete the workbook, or set `force_rerun`, to redo it).

## Output

| Where | What |
|---|---|
| `results/02_population_refinement/<round>_<label>/` | Each round's characterization, DE workbooks, Venns, `pairwise_de/`, `exclusion_check/` (with `<round>_excluded_vs_retained_de_enrichment.xlsx`), round summary |
| `results/02_population_refinement/exploration/` | Parameter screening, candidate inspection, round comparisons |
| `data/processed/02_population_refinement_<round>.rds` | Checkpoint per round (+ `_manifest.txt`; `_config.rds`: the settings used to detect a changed round) |
| `data/processed/cell_lists/` | Cell lists for `exclude_cells` steps |
| `data/processed/02_population_refined.rds` | **Final clustering with a `cluster_label` column: input to later workflows** (+ `_manifest.txt`) |

## Record what you did

Keep a copy of `config/decisions.R` and `config/params_override.R`, and note
the version in the repo's `VERSION` file. Together they reproduce your run.
