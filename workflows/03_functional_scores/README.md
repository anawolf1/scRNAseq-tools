# 03_functional_scores: Functional program scores

Scores every cell for curated functional programs (e.g. exhaustion, memory)
from gene sets, then compares the scores between clusters and timepoints.
Timepoint effects are tested on patient means: longitudinally with a mixed
model, and against a separately sampled group with Wilcoxon tests.

## Start here

**You need:**
- population refinement's output, `data/processed/02_population_refined.rds`
- your gene sets as `.gmt` files in `data/gene_sets/`

**You edit only two files**, both in `config/` (copy them from the `*_example.R` files):

| File | What goes in it |
|---|---|
| `params_override.R` | Your study design (patient/timepoint columns, which timepoints are repeated vs cross-sectional), the species of each gene set, colors, which statistics to run |
| `decisions.R` | Which cells each comparison uses |

> **Do not edit files in `pipeline/`, `functions/` or `shared/` to analyse a new dataset.**

## Run

Open `scRNAseq-tools.Rproj` in RStudio first, then:

```r
source(here::here("workflows", "03_functional_scores", "pipeline", "01_functional_scores.R"))
```

## Choosing cells per comparison

Scores are computed **once, for every cell**. You then decide which cells each
comparison uses, in `decisions.R`, in two parts.

**1. Define cell sets.** Any metadata column can define one:

```r
params$cell_sets <- list(
  all_cells = list(label = "All cells"),
  non_naive = list(label   = "Naive cells excluded",
                   exclude = list(column = "Naive", values = "Naive"))
)
```

`keep` keeps only cells with the given values, `exclude` removes them. Rules
work on any column: a protein-gated population (like a `Naive` column from
index sorting), a cluster label (`column = "cluster_label"`), a patient.
A value that does not occur in the column stops the run and lists the values
that do, so a typo cannot match nothing. Check a column with `table(data_obj$Naive)`.

> **Naive cells vs the naive cluster.** `exclude = list(column = "Naive", ...)`
> removes cells gated as naive, wherever they clustered.
> `exclude = list(column = "cluster_label", values = "Naive")` removes the whole
> naive *cluster*, including any non-naive cells in it and missing naive cells
> elsewhere. They are different analyses: say which one you mean.

**2. Assign a cell set to each comparison:**

```r
params$analysis_cells <- list(
  cluster_summaries        = "all_cells",   # clusters compared on every cell
  program_umaps            = "all_cells",
  timepoint_summaries      = "non_naive",
  trajectory_plot          = "non_naive",
  longitudinal_analysis    = "non_naive",
  cross_sectional_analysis = "non_naive"
)
```

Give several names to run a comparison on each, e.g. as a sensitivity check:
`longitudinal_analysis = c("non_naive", "all_cells")`.

Because every cell keeps the same score in every set, differences between sets
come only from which cells are compared.

## Comparisons

Switched on/off in `params_override.R` with `steps$<name>` (defaults in `config/config.R`):

| Setting | Output | Unit |
|---|---|---|
| `cluster_summaries` | Score table + violins by cluster | cells (descriptive) |
| `program_umaps` | UMAP colored by each program (z-scored) | cells (descriptive) |
| `timepoint_summaries` | Score table + violins by timepoint | cells (descriptive) |
| `trajectory_plot` | Per-patient score lines across repeated timepoints, cross-sectional group alongside | patient means |
| `longitudinal_analysis` (off by default) | Mixed model `score ~ timepoint + (1 \| patient)` on the repeated timepoints: overall effect (likelihood-ratio test, BH across programs) and pairwise contrasts (Kenward-Roger df by default, BH within program); needs `statistics$min_patients_lmm` patients (default 3) | patient means |
| `cross_sectional_analysis` (off by default) | Unpaired Wilcoxon rank-sum, each repeated timepoint vs each cross-sectional group, or the pairs in `statistics$cross_sectional_comparisons` (BH across programs within each comparison); two repeated timepoints cannot be a pair, since they are the same patients | patient means |

**Statistics never use individual cells** as replicates: cells from one patient
are not independent. Every test uses one mean per patient × timepoint;
`params$statistics$min_cells_per_sample` leaves out samples with too few cells
(the `cells_per_sample` table shows what each value would remove).

With few patients, report non-significant results as "not detected", not as
"no difference": e.g. with 4 vs 4 patients the smallest possible Wilcoxon
p-value is 0.029, before any correction.

## Gene sets

Every set in the `.gmt` files needs its species in `params$gene_sets$species`
(`"human"` or `"mouse"`). Human symbols are updated to current names; mouse sets
are mapped to human orthologs. The coverage table shows, per program, how many
genes survived each step and which are missing from the data.

For paired programs (`<name>_UP` and `<name>_DN`), a net score `<name>_NET` = UP − DN
is added (`params$scoring$net_scores`). A shift that raises both halves together
(e.g. samples detecting more genes overall) cancels in the net score.

## Output

| Where | What |
|---|---|
| `results/03_functional_scores/` | Gene-set coverage table |
| `results/03_functional_scores/<cell set>/` | The comparisons that use this cell set, cells per sample, workbook with every table |
| `data/processed/03_functional_scores.rds` | **Every cell with its scores: input to later workflows** (+ `_manifest.txt`) |
| `data/processed/03_term2gene_used.rds` | The resolved gene sets actually scored (+ `_manifest.txt`) |

## Record what you did

Keep a copy of `config/decisions.R` and `config/params_override.R`, the
`.gmt` files, and note the version in the repo's `VERSION` file. Together
they reproduce your run.
