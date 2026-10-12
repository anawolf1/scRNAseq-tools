# 04_neighborhood_analysis: differential abundance and expression in neighborhoods

Tests which cell states change in abundance between timepoints (miloR), and
which genes change within cell states (miloDE). Both work on
**neighborhoods**: small, overlapping groups of transcriptionally similar
cells, so continuous states are captured without relying on cluster borders.

## Start here

**You need:**
- population refinement's output, `data/processed/02_population_refined.rds`
- optionally, `data/processed/03_functional_scores.rds` (programs per neighborhood)

**You edit only two files**, both in `config/` (copy them from the `*_example.R` files):

| File | What goes in it |
|---|---|
| `params_override.R` | Study description, colors, genes hidden from marker tables |
| `decisions.R` | The analyses (which cells build each graph, which contrasts), and the neighborhood parameters you chose |

> **Do not edit files in `pipeline/`, `functions/` or `shared/` to analyze a new dataset.**

## Workflow

Open `scRNAseq-tools.Rproj` in RStudio first, then:

1. **Set the analyses** in `decisions.R` (see [Analysis design](#analysis-design)).
2. **Choose neighborhood parameters** for the primary analysis:
   ```r
   source(here::here("workflows", "04_neighborhood_analysis", "tools", "select_milo_parameters.R"))
   ```
   Part A screens `k` × `d` × `prop` for DA neighborhoods (target: on average ≥
   `params$da$size_multiplier` cells per sample compared), Part B screens `k` for DE
   neighborhoods, Part C screens the Leiden resolution of the neighborhood groups.
   Record the choices in `decisions.R` (`params$da$k`, `params$da$d`, `params$da$prop`,
   `params$de$k`, `params$da$nhood_groups$resolution`); the tool prints the lines to paste.
3. **Run the analysis:**
   ```r
   source(here::here("workflows", "04_neighborhood_analysis", "pipeline", "01_milo_da_de.R"))
   ```

---

## How Milo works — where each number comes from

1. **Graph.** Every cell is linked to its `k` most similar cells, measured in the
   first `d` principal components (the same PCA the clustering used — not the UMAP).
2. **Neighborhoods.** A random share (`prop`) of cells is sampled and refined into
   *index cells*; each neighborhood is an index cell plus its neighbors. A cell can
   belong to several neighborhoods, so they overlap.
3. **Counting.** For each neighborhood and each **sample** (patient × timepoint), the
   number of its cells is counted. This table — neighborhoods × samples — is all the
   test sees.
4. **Test.** Per neighborhood, a negative-binomial model compares the counts of the
   two timepoints, accounting for how many cells each sample contributed in total
   (TMM normalization, `params$da$norm_method`). **The replicates are the samples, not
   the cells**: a neighborhood full of one patient's cells is one observation, however
   many cells it holds.
5. **SpatialFDR.** Neighborhoods overlap, so testing them all would count densely
   sampled regions many times. SpatialFDR is a BH correction that down-weights
   overlapping neighborhoods; significant = SpatialFDR < `params$da$spatial_fdr`.
6. **logFC** is log2 of the (normalized) abundance ratio, test vs reference. It is a
   change in *relative composition*: "enriched in FU" means that state is a larger share
   of FU samples' cells than of SPR samples' cells.

**Neighborhood groups (from all neighborhoods).** Neighborhoods are linked by the
cells they share (Jaccard index) and split into communities with Leiden — like clustering
cells, but applied to neighborhoods. **No DA result is used**, so the groups can then be
tested and described without circularity. Each cell is assigned to the group most of its
neighborhoods belong to, so groups never double-count a cell. The resolution is set in
`params$da$nhood_groups`; nearby resolutions are reported to show how much the grouping
depends on it. Each group gets a **group-level abundance test** (cells per sample, the
same designs), per-sample abundance, composition, maps and a row in the program heatmap.
The same test on the **Seurat clusters** is reported as a comparator: what clustering
alone would conclude.

**Meta-neighborhoods** (kept for comparison). Significant neighborhoods with the same
direction that share ≥ `params$da$markers$min_overlap` of their cells, merged. Because
they are selected by significance, they are descriptive only.

**Cluster label of a neighborhood** = its majority cluster if it holds ≥
`params$da$mixed_threshold` (50 %) of the cells, otherwise "Mixed". An annotation rule, not
a statistical call.

**miloDE.** Separate, larger neighborhoods, built per contrast from its cells
(`params$de$neighborhoods = "milode"`, the default), or the DA neighborhoods themselves
(`"da"`). In each, cells are summed per sample (pseudobulk) and genes are compared
between the two timepoints with edgeR; a gene is DE if its within-neighborhood FDR <
`params$de$alpha`.

## Analysis design

| Analysis | Cells building the graph | Contrasts | Question |
|---|---|---|---|
| **Primary** | BL + 3M + FU + SPR | 3M vs SPR, FU vs SPR, 3M vs FU | Which cell states differ? The same cells the clustering and cluster labels were built on; a neighborhood is the same cells in every contrast. BL shapes the landscape but is not tested (2 patients) |
| **Sensitivity: without BL** | 3M + FU + SPR | the same | Do conclusions depend on BL shaping the landscape? |
| **Sensitivity: pairwise** | only the contrast's two timepoints | each contrast | Do conclusions depend on the other timepoints shaping the neighborhoods? |

Both sensitivity analyses change the **graph construction** (which cells define the
neighborhoods), not the tests: a test always counts only its two timepoints' samples.

The primary analysis is a **common landscape**, not a longitudinal analysis: pairing is
handled in each test's design. **3M vs FU** (the same patients) is tested
**paired**: `~ patient + condition` on patients sampled at both, so each patient is
compared with themself. **3M vs SPR** and **FU vs SPR** (different patients) are
unpaired, `~ condition`. miloDE uses the same designs.

Neighborhood numbers differ between analyses, so they are compared **per cell**: a
cell is "enriched" if it lies in a significant enriched neighborhood, and so on. The
robustness table reports, per contrast, how many of the primary's DA cells get the same
call in each other analysis, and whether any flip direction.

Only the primary analysis gets the full characterization (diagnostics, neighborhood
groups, meta-neighborhoods, markers, miloDE). Other analyses report DA results only.

## Three levels of resolution

| Level | Answers |
|---|---|
| **Neighborhoods** | *Where exactly* in the landscape does abundance change? |
| **Neighborhood groups** | Which coherent regions differ? A summary of the neighborhoods, with fewer tests |
| **Seurat clusters** (comparator) | What would the conventional cluster-level analysis conclude? |

Group-level tests reduce the number of hypotheses; they do not automatically have more
power — aggregating neighborhoods can also dilute a local change.

## Reading the results

| Output (per contrast, primary) | What it shows |
|---|---|
| `da_graph_*` | Neighborhoods at their index cell's UMAP position; size = cells, edges = shared cells; colored if significant; significant ones numbered |
| `da_beeswarm_*` | Every neighborhood's logFC, by majority cluster |
| `da_variability_*` | Why a neighborhood is (not) significant: logFC vs the share of its abundance from a single sample, and vs within-condition variation |
| `da_extreme_nhoods_abundance_*` | For the neighborhoods with the largest \|logFC\|, each sample's % of cells inside — exactly what the test compares |
| `da_sig_nhood_cells_umap_*` | Which cells make up each significant neighborhood |
| `da_qc_vs_logfc_*`, `nhood_qc_vs_programs` | Mean QC metrics per neighborhood against logFC and against program scores |
| `markers/`, `meta_nhoods/` | Meta-neighborhoods: cells on UMAP, markers, results in every contrast, per-sample abundance |

- **Large logFCs (e.g. > 3)** come from neighborhoods with almost no cells in one group:
  trust the direction, not the number.
- **Marker DE of meta-neighborhoods is descriptive**: the cells were selected in the same
  expression space being tested. The within-sample check shows whether a marker is high
  in the group's cells across samples, or only in one.
- **miloDE**: `pct_extreme_logFC` flags genes present in one group only (|logFC| >
  `params$de$extreme_logfc`),
  which usually reflects differences between individuals, not cell state.

## Frequently asked questions

**The UMAP looks SPR-enriched in a region, but no neighborhood there is significant. Why?**
The UMAP shows cells; the test compares samples. Typical reasons, each visible in the outputs:
- **Normalization.** SPR samples may simply contribute more cells overall. Milo compares
  each sample's *share* of cells, so a region can look denser in SPR without SPR being
  enriched there (`sample_table`).
- **One or two samples.** If most of the region's SPR cells come from one patient, it is one
  high value among four, and the test weighs it as such (`da_variability_*`: high
  top-sample share; `da_extreme_nhoods_abundance_*`: one outlying point).
- **Variation within groups.** With 4–5 samples per group, a shift must be consistent
  across samples to be detected; large within-group spread hides it.
- **UMAP distorts density.** Neighborhoods are built in PCA space; UMAP compresses some
  regions and expands others, so apparent density is not cell count.
- **Few replicates limit the smallest attainable p-value.** Non-significant is "not
  detected with these samples", not "no difference".

To look at any region directly, add its neighborhood numbers to
`params$da$markers$nhood_sets`: you get its cells, markers and per-sample abundance.

**Why are there not more DE neighborhoods?** miloDE pseudobulks per sample *within*
each neighborhood: a sample with fewer than `params$de$min_n_cells_per_sample` cells there
leaves the test, and with ~1,300 cells across ~13 samples many neighborhoods keep only a
few samples per group. Neighborhoods are also homogeneous by construction, so
within-state differences are subtler than between-cluster ones, and each neighborhood
is corrected across thousands of genes. The neighborhood summary lists, per
neighborhood, whether it was tested and how many genes were DE.

**Should cells be batch-corrected before building the graph?** Only for technical
batches that cut across the compared groups (e.g. sequencing runs or plates containing
samples of every timepoint). Do not correct by patient: SPR patients are a different
group from the longitudinal patients, so removing patient differences would also remove
the SPR-vs-longitudinal differences being tested. Check first whether batches and
timepoints are crossed (a cross-table of plate × timepoint); if they are confounded,
correction cannot separate them.

## Optional sections

Switched off in `params_override.R` with `steps$<name> <- FALSE` (all on by default):

| Setting | What it adds |
|---|---|
| `milo_da` | miloR differential abundance for every analysis (needed by the next three) |
| `nhood_program_scores` | Functional-program scores per neighborhood: heatmap of significant ones, programs vs logFC across all (descriptive — neighborhoods overlap) |
| `da_nhood_markers` | Meta-neighborhoods and `params$da$markers$nhood_sets`: marker DE, volcanoes, within-sample marker check |
| `nhood_groups` | Groups of all neighborhoods (Leiden). One figure per question: **what are they** (net-program heatmaps of neighborhoods and of groups), **where** (UMAP and neighborhood graph by group), **what are they made of** (cluster and timepoint composition), **do they change** (group- vs cluster-level abundance dot plot, per-sample abundance), **what distinguishes them** (one-vs-rest DE: volcanoes, top-gene tables, workbook) |
| `milo_de` | miloDE per primary contrast |

## Output

| Where | What |
|---|---|
| `results/04_neighborhood_analysis/01_milo_da/<analysis>/` | DA per analysis and contrast; primary: diagnostics, `nhood_groups/`, `meta_nhoods/`, `chosen_nhoods/` |
| `results/04_neighborhood_analysis/01_milo_da/robustness/` | Cell-level DA calls compared between analyses |
| `results/04_neighborhood_analysis/02_milo_de/` | miloDE per primary contrast |
| `results/04_neighborhood_analysis/exploration/` | Parameter screening |
| `data/processed/04_milo_da.rds` (primary), `04_milo_da_<analysis>.rds`, `04_milo_de_<contrast>.rds` | Milo objects and results |
