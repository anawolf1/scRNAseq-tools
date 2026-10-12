# ------------------------------------------------------------------------------
# enrichment.R
#
# Differential expression between two cell groups and GO/KEGG/curated
# gene-set enrichment of the resulting gene lists.
#
# Functions
# run_de_enrichment()           DE between two cell groups, with GO/KEGG/curated enrichment per direction
# enrich_gene_list()            GO and KEGG enrichment of one gene list
# enrich_curated_sets()         over-representation test against a curated gene-set collection
# read_gmt()                    .gmt gene-set file(s) -> term/gene table
# ------------------------------------------------------------------------------

#' DE between two cell groups, with GO/KEGG/curated enrichment per direction
#'
#' Runs FindMarkers() (Wilcoxon) between two values of `group_col`, or takes
#' precomputed DE results, keeps the genes passing `log2fc_cutoff` and
#' `padj_cutoff`, and runs enrich_gene_list() and, with `curated_term2gene`,
#' enrich_curated_sets() on the genes up in each group. Saves a table of the
#' top DE genes, the enrichment dot plots, and one Excel workbook with every
#' tested gene, the reported genes and every enrichment result
#' (<file_prefix>_de_enrichment.xlsx).
#'
#' FindMarkers() only tests genes detected in at least `min_pct` of either
#' group and with |log2FC| >= `log2fc_cutoff`; genes below that could not
#' be reported anyway. Seurat's adjusted p-value is corrected for every gene
#' in the object, so this filter does not change it.
#'
#' The enrichment background (`universe`) defaults to the genes that could
#' have been reported: those detected in at least `min_pct` of either group,
#' minus `exclude_gene_patterns`. GO, KEGG and the curated sets all use it.
#'
#' The Wilcoxon test treats cells as independent observations: its p-values
#' describe differences between the cells compared, not replicate-level
#' evidence. When cells come from several patients, patient structure can
#' still matter; for conditions compared across patients, use a
#' replicate-aware (pseudobulk) method.
#'
#' @param obj Seurat object containing both groups. Not needed when
#'   `de_results` is given.
#' @param group_col Metadata column holding the two groups. Not needed when
#'   `de_results` is given.
#' @param ident_1,ident_2 The two groups compared (values of `group_col`).
#' @param min_pct Minimum fraction of cells in either group expressing a
#'   gene for it to be tested. Default 0.1.
#' @param log2fc_cutoff,padj_cutoff A gene is reported and enriched when
#'   |log2FC| > `log2fc_cutoff` and adjusted p-value < `padj_cutoff`.
#'   Defaults 0.5 and 0.05.
#' @param exclude_gene_patterns Regex patterns of genes to drop (e.g.
#'   unannotated Ensembl IDs). Default "^ENSG".
#' @param label_1,label_2 Group labels for titles and file names. Default
#'   `ident_1`/`ident_2`.
#' @param curated_term2gene Optional term/gene table from read_gmt(); runs
#'   enrich_curated_sets() for each direction. Default NULL.
#' @param de_results Optional precomputed DE results for this comparison
#'   (genes as row names or a `gene` column, plus `avg_log2FC` and
#'   `p_val_adj`). Skips FindMarkers(); `obj`, `group_col`, `min_pct` and
#'   `exclude_gene_patterns` are then not used. The two cutoffs still apply.
#' @param universe Optional background gene symbols for the enrichment tests.
#'   NULL (default): the detected genes described above when FindMarkers()
#'   runs here; each database's default with `de_results`.
#' @param n_terms Number of top terms in each enrichment dot plot. Default 10.
#' @param output_dir Directory for the saved table, plots and workbook.
#' @param file_prefix Prefix of every saved file name.
#'
#' @return Invisibly, a list: `de_tested` (every tested gene), `de_results`
#'   (reported genes, sorted by adjusted p-value), `universe`, and
#'   `enrich_1`/`enrich_2`, each a list with `go`, `kegg` and `curated` (NULL
#'   where nothing was tested or enriched).
#'
#' @examples
#' \dontrun{
#' run_de_enrichment(
#'   obj, group_col = "exclusion_status", ident_1 = "excluded", ident_2 = "retained",
#'   label_1 = "Excluded", label_2 = "Retained",
#'   output_dir = fig_dir, file_prefix = "n3_subcluster"
#' )
#' }
run_de_enrichment <- function(obj = NULL,
                              group_col = NULL,
                              ident_1,
                              ident_2,
                              min_pct = 0.1,
                              log2fc_cutoff = 0.5,
                              padj_cutoff = 0.05,
                              exclude_gene_patterns = c("^ENSG"),
                              label_1 = ident_1,
                              label_2 = ident_2,
                              curated_term2gene = NULL,
                              de_results = NULL,
                              universe = NULL,
                              n_terms = 10,
                              output_dir,
                              file_prefix) {
  
  if (identical(ident_1, ident_2)) {
    stop("run_de_enrichment(): `ident_1` and `ident_2` must be different.", call. = FALSE)
  }
  
  if (is.null(de_results)) {
    
    if (is.null(obj) || is.null(group_col)) {
      stop("run_de_enrichment(): give `obj` and `group_col`, or `de_results`.", call. = FALSE)
    }
    validate_required_cols(obj[[]], group_col, "run_de_enrichment")
    group_values <- obj[[group_col, drop = TRUE]]
    if (!ident_1 %in% group_values || !ident_2 %in% group_values) {
      stop("run_de_enrichment(): both `ident_1` and `ident_2` must be present in `",
           group_col, "`.", call. = FALSE)
    }
    
    # the object may not be normalized yet (e.g. before QC filtering)
    obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
    Seurat::Idents(obj) <- group_values
    
    de_markers <- tryCatch(
      Seurat::FindMarkers(obj, ident.1 = ident_1, ident.2 = ident_2, test.use = "wilcox",
                          min.pct = min_pct, logfc.threshold = log2fc_cutoff),
      error = function(e) {
        warning("run_de_enrichment(): FindMarkers() failed for ", label_1, " vs. ", label_2,
                ": ", conditionMessage(e), call. = FALSE)
        data.frame(p_val = numeric(0), avg_log2FC = numeric(0), pct.1 = numeric(0),
                   pct.2 = numeric(0), p_val_adj = numeric(0))
      }
    )
    
    # background: genes detected in at least min_pct of either group, the
    # same eligibility FindMarkers() applies
    if (is.null(universe)) {
      counts <- Seurat::GetAssayData(obj, layer = "counts")
      pct_1 <- Matrix::rowMeans(counts[, group_values == ident_1, drop = FALSE] > 0)
      pct_2 <- Matrix::rowMeans(counts[, group_values == ident_2, drop = FALSE] > 0)
      universe <- rownames(counts)[pmax(pct_1, pct_2) >= min_pct]
    }
    
    if (length(exclude_gene_patterns) > 0) {
      exclude_regex <- paste(exclude_gene_patterns, collapse = "|")
      de_markers <- de_markers[!grepl(exclude_regex, rownames(de_markers)), , drop = FALSE]
      universe <- universe[!grepl(exclude_regex, universe)]
    }
    
  } else {
    
    if (!is.data.frame(de_results)) {
      stop("run_de_enrichment(): `de_results` must be a data.frame.", call. = FALSE)
    }
    validate_required_cols(de_results, c("avg_log2FC", "p_val_adj"), "run_de_enrichment")
    de_markers <- de_results
    
    if ("gene" %in% colnames(de_markers)) {
      if (anyNA(de_markers$gene) || anyDuplicated(de_markers$gene)) {
        stop("run_de_enrichment(): `de_results$gene` has missing or duplicated genes; ",
             "pass the results of one comparison only.", call. = FALSE)
      }
      rownames(de_markers) <- de_markers$gene
      de_markers$gene <- NULL
    } else if (identical(rownames(de_markers), as.character(seq_len(nrow(de_markers))))) {
      stop("run_de_enrichment(): `de_results` needs gene names as row names or a 'gene' column.",
           call. = FALSE)
    }
  }
  
  de_tested <- de_markers[order(de_markers$p_val_adj, na.last = TRUE), , drop = FALSE]
  de_markers <- de_markers[which(de_markers$p_val_adj < padj_cutoff &
                                   abs(de_markers$avg_log2FC) > log2fc_cutoff), , drop = FALSE]
  de_markers <- de_markers[order(de_markers$p_val_adj, na.last = TRUE), , drop = FALSE]
  
  if (nrow(de_markers) > 0) {
    de_table <- gt_de(
      de_markers, title = paste(label_1, "vs.", label_2),
      subtitle = "Differentially expressed genes",
      pct1_label = paste0("% expr. (", label_1, ")"),
      pct2_label = paste0("% expr. (", label_2, ")")
    )
    print(de_table)
    save_gt(de_table, paste0(file_prefix, "_de_table.png"), output_dir)
  } else {
    message("run_de_enrichment(): no DE genes passed the cutoffs for ", label_1,
            " vs. ", label_2, "; no DE table.")
  }
  
  up_in_1 <- rownames(de_markers)[de_markers$avg_log2FC > 0]
  up_in_2 <- rownames(de_markers)[de_markers$avg_log2FC < 0]
  
  go_kegg_1 <- enrich_gene_list(up_in_1, label = label_1, output_dir = output_dir,
                                file_prefix = file_prefix, universe = universe, n_terms = n_terms)
  go_kegg_2 <- enrich_gene_list(up_in_2, label = label_2, output_dir = output_dir,
                                file_prefix = file_prefix, universe = universe, n_terms = n_terms)
  
  curated_1 <- NULL
  curated_2 <- NULL
  if (!is.null(curated_term2gene)) {
    curated_1 <- enrich_curated_sets(up_in_1, curated_term2gene, label = label_1, output_dir = output_dir,
                                     file_prefix = file_prefix, universe = universe, n_terms = n_terms)
    curated_2 <- enrich_curated_sets(up_in_2, curated_term2gene, label = label_2, output_dir = output_dir,
                                     file_prefix = file_prefix, universe = universe, n_terms = n_terms)
  }
  
  gene_table <- function(df) tibble::rownames_to_column(as.data.frame(df), "gene")
  sheets <- list("DE tested" = gene_table(de_tested), "DE reported" = gene_table(de_markers))
  enrichments <- list(
    GO      = list(go_kegg_1$go,   go_kegg_2$go),
    KEGG    = list(go_kegg_1$kegg, go_kegg_2$kegg),
    Curated = list(curated_1,      curated_2)
  )
  for (source in names(enrichments)) {
    for (i in 1:2) {
      result <- enrichments[[source]][[i]]
      if (!is.null(result)) {
        sheets[[paste(source, "up in", c(label_1, label_2)[i])]] <- as.data.frame(result)
      }
    }
  }
  save_xlsx(sheets, paste0(file_prefix, "_de_enrichment.xlsx"), output_dir)
  
  invisible(list(
    de_tested  = de_tested,
    de_results = de_markers,
    universe   = universe,
    enrich_1   = list(go = go_kegg_1$go, kegg = go_kegg_1$kegg, curated = curated_1),
    enrich_2   = list(go = go_kegg_2$go, kegg = go_kegg_2$kegg, curated = curated_2)
  ))
}
#' GO and KEGG enrichment of one gene list
#'
#' GO biological process (enrichGO()) on gene symbols and KEGG
#' (enrichKEGG(), after mapping to Entrez IDs), both BH-adjusted. Saves a
#' dot plot of the top terms of each. enrichKEGG() downloads KEGG data; if
#' that fails (e.g. offline), KEGG is skipped with a warning and GO is still
#' returned.
#'
#' @param genes Human gene symbols to test.
#' @param label Label for the plot title and file name.
#' @param output_dir Directory for the plots.
#' @param file_prefix Prefix of the saved file names.
#' @param universe Optional background gene symbols: the genes that could
#'   have been called DE (e.g. every gene FindMarkers() tested). NULL
#'   (default) uses every gene annotated in the database; pass it when
#'   comparing results between analyses that tested different genes.
#' @param qvalue_cutoff Significance threshold (p and q value). Default 0.05.
#' @param n_terms Number of top terms in each dot plot. Default 10.
#'
#' @return Invisibly, a list with `go` and `kegg` (enrichResult objects, or
#'   NULL when there were no genes, no enriched terms, or KEGG failed).
#'
#' @examples
#' \dontrun{
#' enrich_gene_list(naive_up_genes, label = "Naive",
#'                  output_dir = fig_dir, file_prefix = "whole_cluster_naive")
#' }
enrich_gene_list <- function(genes,
                             label,
                             output_dir,
                             file_prefix,
                             universe = NULL,
                             qvalue_cutoff = 0.05,
                             n_terms = 10) {
  
  if (length(genes) == 0) {
    message("[skip] no genes to test for: ", label)
    return(invisible(list(go = NULL, kegg = NULL)))
  }
  
  go_result <- clusterProfiler::enrichGO(
    gene = genes, universe = universe, OrgDb = org.Hs.eg.db::org.Hs.eg.db, keyType = "SYMBOL",
    ont = "BP", pAdjustMethod = "BH", pvalueCutoff = qvalue_cutoff, qvalueCutoff = qvalue_cutoff
  )

  if (is.null(go_result) || nrow(as.data.frame(go_result)) == 0) {
    message("[no enrichment] GO for: ", label)
    go_result <- NULL
  } else {
    go_result@result <- go_result@result[order(go_result@result$p.adjust, na.last = TRUE), , drop = FALSE]
    go_plot <- enrichplot::dotplot(go_result, showCategory = n_terms,
                                   title = paste("Top enriched GO biological processes -", label))
    print(go_plot)
    save_plot(go_plot, paste0(file_prefix, "_go_", label, ".png"), output_dir, width = 8, height = 6)
  }
  
  entrez_ids <- clusterProfiler::bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                                      OrgDb = org.Hs.eg.db::org.Hs.eg.db)
  
  kegg_result <- NULL
  if (nrow(entrez_ids) == 0) {
    message("[skip] no genes mapped to Entrez IDs for KEGG: ", label)
  } else {
    universe_entrez <- if (!is.null(universe)) {
      clusterProfiler::bitr(universe, fromType = "SYMBOL", toType = "ENTREZID",
                            OrgDb = org.Hs.eg.db::org.Hs.eg.db)$ENTREZID
    }
    kegg_result <- tryCatch(
      clusterProfiler::enrichKEGG(
        gene = entrez_ids$ENTREZID, universe = universe_entrez, organism = "hsa",
        pAdjustMethod = "BH", pvalueCutoff = qvalue_cutoff, qvalueCutoff = qvalue_cutoff
      ),
      error = function(e) {
        warning("enrich_gene_list(): KEGG enrichment failed for ", label, " (KEGG needs an internet ",
                "connection): ", conditionMessage(e), call. = FALSE)
        NULL
      }
    )
    if (is.null(kegg_result) || nrow(as.data.frame(kegg_result)) == 0) {
      message("[no enrichment] KEGG for: ", label)
      kegg_result <- NULL
    } else {
      kegg_result@result <- kegg_result@result[order(kegg_result@result$p.adjust, na.last = TRUE), , drop = FALSE]
      kegg_plot <- enrichplot::dotplot(kegg_result, showCategory = n_terms,
                                       title = paste("Top enriched KEGG pathways -", label))
      print(kegg_plot)
      save_plot(kegg_plot, paste0(file_prefix, "_kegg_", label, ".png"), output_dir, width = 8, height = 6)
    }
  }
  
  invisible(list(go = go_result, kegg = kegg_result))
}
#' Over-representation test against a curated gene-set collection
#'
#' clusterProfiler::enricher() with a user-supplied collection (e.g.
#' published exhaustion/memory signatures not in GO or KEGG). Saves a dot
#' plot of the top sets.
#'
#' @param genes Gene symbols to test.
#' @param term2gene Data.frame with `term` (set name) and `gene` columns,
#'   one row per set/gene pair, e.g. from read_gmt().
#' @param label Label for the plot title and file name.
#' @param output_dir Directory for the plot.
#' @param file_prefix Prefix of the saved file name.
#' @param universe Optional background gene symbols (see enrich_gene_list()).
#'   NULL (default) uses every gene in `term2gene`.
#' @param qvalue_cutoff Significance threshold (p and q value). Default 0.05.
#' @param n_terms Number of top sets in the dot plot. Default 10.
#'
#' @return Invisibly, the enrichResult, or NULL when there were no genes or
#'   no enriched sets.
#'
#' @examples
#' \dontrun{
#' curated_sets <- read_gmt(list.files("data/gene_sets", pattern = "\\.gmt$", full.names = TRUE))
#' enrich_curated_sets(fulike_up_genes, curated_sets, label = "FU-like",
#'                     output_dir = fig_dir, file_prefix = "whole_cluster_fulike")
#' }
enrich_curated_sets <- function(genes,
                                term2gene,
                                label,
                                output_dir,
                                file_prefix,
                                universe = NULL,
                                qvalue_cutoff = 0.05,
                                n_terms = 10) {
  
  validate_required_cols(term2gene, c("term", "gene"), "enrich_curated_sets")
  
  if (length(genes) == 0) {
    message("[skip] no genes to test for: ", label)
    return(invisible(NULL))
  }
  
  curated_result <- clusterProfiler::enricher(
    gene = genes, universe = universe, TERM2GENE = term2gene, pAdjustMethod = "BH",
    pvalueCutoff = qvalue_cutoff, qvalueCutoff = qvalue_cutoff
  )
  if (is.null(curated_result) || nrow(as.data.frame(curated_result)) == 0) {
    message("[no enrichment] curated gene sets for: ", label)
    return(invisible(NULL))
  }
  curated_result@result <- curated_result@result[order(curated_result@result$p.adjust, na.last = TRUE), ,
                                                 drop = FALSE]
  
  curated_plot <- enrichplot::dotplot(curated_result, showCategory = n_terms,
                                      title = paste("Top curated gene sets -", label))
  print(curated_plot)
  save_plot(curated_plot, paste0(file_prefix, "_curated_", label, ".png"), output_dir, width = 8, height = 6)
  
  invisible(curated_result)
}
#' .gmt gene-set file(s) -> term/gene table
#'
#' One gene set per line, tab-separated: name, description (ignored), genes.
#' Several files are combined into one table. Gene symbols are kept as they
#' are, with no species conversion; check each collection's species in its
#' source, since symbol case is not a reliable indicator.
#'
#' @param paths .gmt file paths.
#'
#' @return A data.frame with `term` and `gene`, one row per set/gene pair.
#'
#' @examples
#' \dontrun{
#' term2gene <- read_gmt(list.files("data/gene_sets", pattern = "\\.gmt$", full.names = TRUE))
#' }
read_gmt <- function(paths) {
  
  rows <- lapply(paths, function(path) {
    lines <- readLines(path, warn = FALSE)
    do.call(rbind, lapply(lines, function(line) {
      if (!nzchar(trimws(line))) return(NULL)
      fields <- strsplit(line, "\t", fixed = TRUE)[[1]]
      if (length(fields) < 3 || !nzchar(fields[1])) return(NULL)
      genes <- fields[-c(1, 2)]
      genes <- genes[nzchar(genes)]
      if (length(genes) == 0) return(NULL)
      data.frame(term = fields[1], gene = genes, stringsAsFactors = FALSE)
    }))
  })
  
  result <- do.call(rbind, rows)
  if (is.null(result)) {
    return(data.frame(term = character(0), gene = character(0), stringsAsFactors = FALSE))
  }
  rownames(result) <- NULL
  result
}