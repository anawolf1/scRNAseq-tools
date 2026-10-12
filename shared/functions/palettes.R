# ------------------------------------------------------------------------------
# palettes.R
#
# Registry of which color palette each context (a kind of plot or table)
# uses, so choosing a different one is a one-line override in
# params_override.R rather than an edit to a plotting function.
#
# Five kinds of context:
#   - "continuous": a viridis-family palette name, for scale_*_viridis_c(option = ...)
#   - "discrete":   a positional hex vector, or a named palette matched by value
#   - "diverging":  a list(low =, mid =, high =) hex triplet
#   - "sequential": a list(low =, high =) pair - like diverging, no midpoint
#   - "status":     a named hex vector for fixed semantic states
#
# Sections
#   1. Continuous palettes
#   2. Discrete palettes
#   3. Diverging palettes
#   4. Sequential palettes
#   5. Status palettes
#   6. Registry: default palette per context
#   7. Functions
#
# Usage:
#   (sourced automatically by shared/functions/load_all.R)
#   resolve_palette("qc_scatter")                                # "plasma"
#   resolve_palette("cluster")                                   # Okabe-Ito hex vector
#   resolve_palette("group_umap", n = 4)                         # 4 colors from carto_safe
#   resolve_palette("count_heatmap")                             # list(low=, mid=, high=)
#   resolve_palette("sweep_heatmap")                             # list(low=, high=)
#   resolve_palette("qc_status")                                 # c(retained=, excluded=) named vector
#   palette_help()                                               # print the whole registry
#   palette_help("cluster")                                      # filter by context/description
#   add_palette("lab_house", c("#1B9E77", "#D95F02", "#7570B3")) # add a custom discrete palette
#   add_status_palette("patients", c(P001 = "#E41A1C"))          # add a fixed color per value
#
# To use a different palette for a context, set it in params_override.R:
#   params$palettes$cluster    <- "dark2"
#   params$palettes$qc_scatter <- "magma"
#
# Functions
# resolve_palette()         the palette (name, colors or list) for a context
# resolve_group_colors()    named colors for the levels of a grouping column
# palette_help()            print the registry, optionally filtered
# add_palette()             register a custom discrete palette
# add_status_palette()      register a named palette (fixed color per value)
# resolve_cluster_colors()  stable colors for cluster labels
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 1. Continuous palettes
# ------------------------------------------------------------------------------
continuous_palette_catalog <- list(
  viridis = c(
    "#440154", "#414487", "#2A788E",
    "#22A884", "#7AD151", "#FDE725"
  ),
  magma = c(
    "#000004", "#3B0F70", "#8C2981",
    "#DE4968", "#FE9F6D", "#FCFDBF"
  ),
  plasma = c(
    "#0D0887", "#6A00A8", "#B12A90",
    "#E16462", "#FCA636", "#F0F921"
  ),
  inferno = c(
    "#000004", "#420A68", "#932667",
    "#DD513A", "#FCA50A", "#FCFFA4"
  ),
  cividis = c(
    "#00204D", "#31446B", "#666970",
    "#958F78", "#CBBA69", "#FFEA46"
  )
)

# ------------------------------------------------------------------------------
# 2. Discrete palettes
# ------------------------------------------------------------------------------
# "viridis_discrete" is generated on demand from cividis at whatever length
# is requested, so it has no fixed vector here.
#
#   okabe_ito:    Okabe & Ito (2008), designed for color-vision accessibility
#   dark2:        ColorBrewer Dark2, marked colorblind-safe by ColorBrewer
#   tol_muted:    Paul Tol's muted qualitative scheme
#   carto_safe:   CARTOColors "Safe" qualitative palette (12 colors)
#   colorbrewer1: ColorBrewer Set1 - NOT colorblind-safe (red/green, pale yellow)
discrete_palette_catalog <- list(
  okabe_ito = c(
    "#E69F00", "#56B4E9", "#009E73",
    "#F0E442", "#0072B2", "#D55E00",
    "#CC79A7"
  ),
  dark2 = c(
    "#1B9E77", "#D95F02", "#7570B3",
    "#E7298A", "#66A61E", "#E6AB02",
    "#A6761D", "#666666"
  ),
  tol_muted = c(
    "#CC6677", "#332288", "#DDCC77",
    "#117733", "#88CCEE", "#882255",
    "#44AA99", "#999933", "#AA4499"
  ),
  carto_safe = c(
    "#88CCEE", "#CC6677", "#DDCC77",
    "#117733", "#332288", "#AA4499",
    "#44AA99", "#999933", "#882255",
    "#661100", "#6699CC", "#888888"
  ),
  colorbrewer1 = c(
    "#E41A1C", "#377EB8", "#4DAF4A",
    "#984EA3", "#FF7F00", "#FFFF33",
    "#A65628", "#F781BF", "#999999"
  )
)

# ------------------------------------------------------------------------------
# 3. Diverging palettes
# ------------------------------------------------------------------------------
diverging_palette_catalog <- list(
  redblue = list(
    low = "#6495ED",
    mid = "#FFFFFF",
    high = "#B22222"
  )
)

# ------------------------------------------------------------------------------
# 4. Sequential palettes
# ------------------------------------------------------------------------------
# Like diverging, but a plain low-to-high gradient with no meaningful
# midpoint (e.g. marker/stability counts, which have no "bad below / good
# above" center the way a percentage does).
sequential_palette_catalog <- list(
  grey_red  = list(low = "grey95", high = "#B22222"),
  grey_blue = list(low = "grey95", high = "#0571B0")
)

# ------------------------------------------------------------------------------
# 5. Status palettes
# ------------------------------------------------------------------------------
status_palette_catalog <- list(
  cornflower_firebrick = c(retained = "#6495ED", excluded = "#B22222")
)

# ------------------------------------------------------------------------------
# 6. Registry: default palette per context
# ------------------------------------------------------------------------------
palette_registry <- tibble::tribble(
  ~context,                   ~kind,        ~default,                     ~description,
  "cluster",                  "discrete",   "okabe_ito",                  "Cluster ID color - UMAP, composition bars, and anywhere else clusters are colored; always the same mapping",
  "count_heatmap",            "diverging",  "redblue",                    "Count/percent heatmaps (plot_count_heatmap) - low/mid/high colors when diverging = TRUE",
  "count_heatmap_sequential", "sequential", "grey_red",                   "Count/percent heatmaps (plot_count_heatmap) - low/high colors when diverging = FALSE (the function's default mode)",
  "group_umap",               "discrete",   "carto_safe",                 "UMAP colored by a categorical variable with no dedicated context",
  "qc_scatter",               "continuous", "plasma",                     "QC scatter plots (nCount vs nFeature)",
  "sweep_heatmap",            "sequential", "grey_blue",                  "Parameter-sweep heatmaps (plot_sweep_tile)",
  "de_log2fc",                "diverging",  "redblue",                    "Differential-expression log2FC coloring - blue/white/firebrick centered at zero",
  "module_score_z",           "diverging",  "redblue",                    "Module-score z-score coloring - blue/white/firebrick centered at zero",
  "qc_status",                "status",     "cornflower_firebrick",       "Retained/excluded status (plot_plate_grid)",
  "qc_metric",                "continuous", "viridis",                    "Continuous QC metric coloring (plot_plate_grid)",
  "venn",                     "discrete",   "okabe_ito",                  "DEG-overlap Venn diagrams (compare_deg_venn) - one fill color per condition/set",
  "group_var1",               "discrete",   "colorbrewer1",               "First configured grouping variable",
  "group_var2",               "discrete",   "tol_muted",                  "Second configured grouping variable"
)

# ------------------------------------------------------------------------------
# 7. Functions
# ------------------------------------------------------------------------------

#' Resolve which palette to use for a given context
#'
#' Checks `params$palettes[[context]]` first (set in params_override.R),
#' falling back to `palette_registry`'s default if nothing was overridden.
#' A discrete context also accepts a status palette (a fixed color per
#' value, see add_status_palette()).
#'
#' @param context Character. Must match a `context` entry in palette_registry.
#' @param n Integer. Number of colors needed, for discrete contexts:
#'   required for "viridis_discrete"; a fixed palette with fewer colors is
#'   extended with generic colors (with a message). Ignored for
#'   "continuous", "diverging", "sequential" and "status" contexts.
#'
#' @return A viridis option string ("continuous"), a hex color vector
#'   ("discrete"), a list(low =, mid =, high =) ("diverging"), a
#'   list(low =, high =) ("sequential"), or a named hex vector ("status").
#'
#' @examples
#' \dontrun{
#' plot_qc_scatter(data, color_option = resolve_palette("qc_scatter"))
#' plot_umap_by_group(data, "Timepoint", group_colors = setNames(
#'   resolve_palette("group_umap", n = 4), c("BL", "3M", "FU", "SPR")))
#' plot_sweep_tile(sweep_results, "mean_ari_stability", "resolution", "dims", title = "ARI")
#' }
resolve_palette <- function(context,
                            n = NULL) {

  if (!is.null(n) &&
      (!is.numeric(n) || length(n) != 1 || !is.finite(n) ||
       n < 1 || n != as.integer(n))) {
    stop("resolve_palette(): `n` must be a single positive integer.", call. = FALSE)
  }

  match_row <- palette_registry[palette_registry$context == context, ]
  if (nrow(match_row) == 0) {
    stop("resolve_palette(): unknown palette context '", context,
         "'. Check palette_registry in shared/functions/palettes.R.", call. = FALSE)
  }

  # a user override in params_override.R always wins over the registry default
  requested <- if (!is.null(params$palettes[[context]])) {
    params$palettes[[context]]
  } else {
    match_row$default
  }

  if (match_row$kind == "continuous") {
    if (!requested %in% names(continuous_palette_catalog)) {
      stop(
        "resolve_palette(): unknown continuous palette '", requested, "' for context '", context, "'. ",
        "Available palettes: ", paste(names(continuous_palette_catalog), collapse = ", "),
        call. = FALSE
      )
    }
    return(requested)
  }

  if (match_row$kind == "diverging") {
    if (!requested %in% names(diverging_palette_catalog)) {
      stop(
        "resolve_palette(): unknown diverging palette '", requested, "' for context '", context, "'. ",
        "Available palettes: ", paste(names(diverging_palette_catalog), collapse = ", "),
        call. = FALSE
      )
    }
    return(diverging_palette_catalog[[requested]])
  }

  if (match_row$kind == "sequential") {
    if (!requested %in% names(sequential_palette_catalog)) {
      stop(
        "resolve_palette(): unknown sequential palette '", requested, "' for context '", context, "'. ",
        "Available palettes: ", paste(names(sequential_palette_catalog), collapse = ", "),
        call. = FALSE
      )
    }
    return(sequential_palette_catalog[[requested]])
  }

  if (match_row$kind == "status") {
    if (!requested %in% names(status_palette_catalog)) {
      stop(
        "resolve_palette(): unknown status palette '", requested, "' for context '", context, "'. ",
        "Available palettes: ", paste(names(status_palette_catalog), collapse = ", "),
        call. = FALSE
      )
    }
    return(status_palette_catalog[[requested]])
  }

  # discrete context: a status palette gives a fixed color per value
  if (requested %in% names(status_palette_catalog)) {
    return(status_palette_catalog[[requested]])
  }

  if (requested == "viridis_discrete") {
    if (is.null(n)) {
      stop("resolve_palette(): context '", context, "' needs `n` (number of colors) for a ",
           "viridis_discrete palette.", call. = FALSE)
    }
    return(scales::viridis_pal(option = "cividis")(n))
  }

  if (!requested %in% names(discrete_palette_catalog)) {
    stop(
      "resolve_palette(): unknown discrete palette '", requested, "' for context '", context, "'. ",
      "Available palettes: ",
      paste(c(names(discrete_palette_catalog), "viridis_discrete", names(status_palette_catalog)), collapse = ", "),
      call. = FALSE
    )
  }

  palette <- discrete_palette_catalog[[requested]]

  if (!is.null(n)) {
    if (n > length(palette)) {
      message("resolve_palette(): '", requested, "' has ", length(palette), " colors but ", n,
              " are needed for '", context, "' - adding ", n - length(palette), " generic colors.")
      palette <- c(palette, scales::hue_pal()(n - length(palette)))
    }
    palette <- palette[seq_len(n)]
  }

  palette
}
#' Colors for the levels of a grouping column
#'
#' The first two columns of params$group_by_vars use the "group_var1" and
#' "group_var2" palette contexts (e.g. Patient -> colorbrewer1,
#' Timepoint -> tol_muted); any other column uses "group_umap". Colors are
#' assigned in the order of params$display_order[[group_col]] when it is set,
#' so a level keeps the same color in every plot, even when some levels are
#' absent (e.g. BL in the Milo analysis). Named palettes are matched by level.
#'
#' @param group_col Metadata column name.
#' @param group_levels Levels to color, in plot order.
#'
#' @return Named character vector, one color per level.
resolve_group_colors <- function(group_col,
                                 group_levels) {

  group_vars <- params$group_by_vars %||% character(0)
  position <- match(group_col, group_vars)
  context <- if (!is.na(position) && position <= 2) paste0("group_var", position) else "group_umap"

  all_levels <- union(params$display_order[[group_col]], group_levels)
  pal <- resolve_palette(context)

  # named palette: colors are matched to the plotted levels by name
  if (!is.null(names(pal))) {
    missing_levels <- setdiff(group_levels, names(pal))
    if (length(missing_levels) > 0) {
      stop("resolve_group_colors(): the named palette for '", group_col,
           "' has no color for: ", paste(missing_levels, collapse = ", "),
           call. = FALSE)
    }
    return(pal[group_levels])
  }

  # positional palette: needs one color per level in display order
  if (length(pal) < length(all_levels)) {
    warning("resolve_group_colors(): the '", context, "' palette has fewer colors than the ",
            length(all_levels), " levels of '", group_col, "' - using 'group_umap' instead.",
            call. = FALSE)
    pal <- resolve_palette("group_umap")
    if (length(pal) < length(all_levels)) {
      stop("resolve_group_colors(): '", group_col, "' has more levels (", length(all_levels),
           ") than any configured palette provides.",
           call. = FALSE)
    }
  }

  setNames(unname(pal)[seq_along(all_levels)], all_levels)[group_levels]
}
#' Look up which palette a context uses, or search the registry
#'
#' @param pattern Optional regex (case-insensitive), matched against both
#'   `context` and `description`. Can be several patterns (rows matching
#'   any are returned). NULL (default) prints the full registry.
#'
#' @return Invisibly, the (possibly filtered) registry data.frame.
#'
#' @examples
#' \dontrun{
#' palette_help()
#' palette_help("cluster")
#' palette_help(c("heatmap", "umap"))
#' }
palette_help <- function(pattern = NULL) {
  df <- palette_registry
  if (!is.null(pattern)) {
    pattern <- paste(pattern, collapse = "|")
    match_context <- grepl(pattern, df$context, ignore.case = TRUE)
    match_desc    <- grepl(pattern, df$description, ignore.case = TRUE)
    df <- df[match_context | match_desc, ]
  }
  print(df, n = Inf)
  invisible(df)
}
#' Add a custom named discrete palette to the catalog
#'
#' @param name Character. New catalog name.
#' @param colors Character vector of 6-digit hex colors.
#'
#' @return Invisibly, `colors`.
#'
#' @examples
#' \dontrun{
#' add_palette("lab_house", c("#1B9E77", "#D95F02", "#7570B3"))
#' params$palettes$cluster <- "lab_house"
#' }
add_palette <- function(name,
                        colors) {
  if (!is.character(name) || length(name) != 1 || is.na(name) || !nzchar(name)) {
    stop("add_palette(): `name` must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.character(colors) || length(colors) == 0 ||
      !all(grepl("^#[0-9A-Fa-f]{6}$", colors))) {
    stop("add_palette(): `colors` must be a non-empty character vector of 6-digit hex codes, ",
         "e.g. '#1B9E77'.", call. = FALSE)
  }
  all_names <- c(names(continuous_palette_catalog), names(discrete_palette_catalog),
                 names(diverging_palette_catalog), names(sequential_palette_catalog),
                 names(status_palette_catalog))
  if (name %in% all_names) {
    stop("add_palette(): palette '", name, "' already exists. Choose a new name rather than ",
         "overwriting an existing palette.", call. = FALSE)
  }
  discrete_palette_catalog[[name]] <<- colors
  message("Added palette '", name, "' (", length(colors), " colors) to discrete_palette_catalog.")
  invisible(colors)
}
#' Add a fixed, named palette for dataset-specific values (patient IDs, etc)
#'
#' Companion to add_palette(): a palette where each value has a fixed color,
#' matched by name rather than assigned in order. Usable for "status"
#' contexts and for discrete contexts such as group_var1 (e.g. one fixed
#' color per patient).
#'
#' @param name Character. New catalog name, e.g. "patient_dataset_specific".
#' @param colors A named character vector of hex colors, 6 or 8 digits (8
#'   includes an alpha channel, e.g. "#1B9E77FF"). Names are the fixed
#'   values this palette applies to (patient IDs, timepoint labels).
#'
#' @return Invisibly, `colors`.
#'
#' @examples
#' \dontrun{
#' add_status_palette("patient_dataset_specific", c(
#'   "P001" = "#E41A1C", "P002" = "#377EB8", "P003" = "#4DAF4A"
#' ))
#' params$palettes$group_var1 <- "patient_dataset_specific"
#' }
add_status_palette <- function(name,
                               colors) {
  if (!is.character(name) || length(name) != 1 || is.na(name) || !nzchar(name)) {
    stop("add_status_palette(): `name` must be a single non-empty character string.", call. = FALSE)
  }
  if (!is.character(colors) || length(colors) == 0 ||
      is.null(names(colors)) || anyNA(names(colors)) ||
      any(!nzchar(names(colors))) || anyDuplicated(names(colors))) {
    stop("add_status_palette(): `colors` must be a named character vector with unique, non-empty ",
         "names (the fixed ID/label each color applies to).",
         call. = FALSE)
  }
  if (!all(grepl("^#[0-9A-Fa-f]{6}([0-9A-Fa-f]{2})?$", colors))) {
    stop("add_status_palette(): `colors` must contain 6- or 8-digit hex codes, e.g. '#1B9E77' or ",
         "'#1B9E77FF'.", call. = FALSE)
  }
  all_names <- c(names(continuous_palette_catalog), names(discrete_palette_catalog),
                 names(diverging_palette_catalog), names(sequential_palette_catalog),
                 names(status_palette_catalog))
  if (name %in% all_names) {
    stop("add_status_palette(): palette '", name, "' already exists. Choose a new name rather than ",
         "overwriting an existing palette.", call. = FALSE)
  }
  status_palette_catalog[[name]] <<- colors
  message("Added palette '", name, "' (", length(colors), " values) to status_palette_catalog.")
  invisible(colors)
}
#' Colors for cluster labels
#'
#' A label in `cluster_color_map` gets that color, so a population keeps its
#' color across rounds and workflows even when its cluster number changes.
#' The other labels take the "cluster" palette in order, without the colors
#' `cluster_color_map` uses (generic colors are added if the palette runs
#' out).
#'
#' @param cluster_labels Cluster labels in plot order (NA for an unlabeled
#'   cluster).
#' @param cluster_color_map Optional named vector, label -> color, e.g.
#'   params$cluster_color_map.
#'
#' @return Named character vector, label -> color.
#'
#' @examples
#' \dontrun{
#' resolve_cluster_colors(c("FU-like", "Naive", "CD127-"), c("Naive" = "#009E73"))
#' }
resolve_cluster_colors <- function(cluster_labels,
                                   cluster_color_map = NULL) {
  fixed <- !is.na(cluster_labels) & cluster_labels %in% names(cluster_color_map)
  n_other <- sum(!fixed)
  palette <- setdiff(resolve_palette("cluster"), unname(cluster_color_map))
  if (length(palette) < n_other) palette <- c(palette, scales::hue_pal()(n_other - length(palette)))
  
  colors <- character(length(cluster_labels))
  colors[fixed]  <- unname(cluster_color_map[cluster_labels[fixed]])
  colors[!fixed] <- palette[seq_len(n_other)]
  setNames(colors, cluster_labels)
}
