#' Create pseudobulk metadata
#'
#' Summarizes metadata at the replicate-level. This is amenable for covariate correction using pseudobulk DE.
#'
#' @param input a single-cell matrix to be converted, with features (genes) in rows
#'   and cells in columns. Alternatively, a \code{Seurat}, \code{monocole3}, or
#'   or \code{SingleCellExperiment} object can be directly input.
#' @param meta the accompanying meta data whereby the rownames match the column
#'   names of \code{input}.
#' @param replicate_col the vector in \code{meta} containing the replicate
#'   information. Defaults to \code{replicate}.
#' @param cell_type_col the vector in \code{meta} containing the cell type
#'   information. Defaults to \code{cell_type}.
#' @param label_col the vector in \code{meta} containing the experimental
#'   label. Defaults to \code{label}.
#' @param min_cells the minimum number of cells in a cell type to retain it.
#'   Defaults to \code{3}.
#' @param min_reps the minimum number of replicates in a cell type to retain it.
#'   Defaults to \code{2}.
#' @param covariates character vector of covariates. Defaults to \code{NULL}.
#'
#' @return a list of dataframes with metadata summariozed at the replicate level, for each cell type.
#' @importFrom magrittr %<>%
#' @importFrom tibble remove_rownames
#' @importFrom dplyr %>% mutate n_distinct filter count group_by select summarise pull across all_of arrange
#' @importFrom purrr map
#'
#' @examples
to_pseudobulk_meta <- function(input,
                               meta = meta,
                               replicate_col = "replicate",
                               cell_type_col = "cell_type",
                               label_col = "label",
                               min_cells = 3,
                               min_reps = 2,
                               covariates = NULL) {
  library(dplyr)
  library(tidyr)

  # check inputs
  inputs <- check_inputs(
    input,
    meta = meta,
    replicate_col = replicate_col,
    cell_type_col = cell_type_col,
    label_col = label_col,
    covariates = covariates
  )
  expr <- inputs$expr
  meta <- inputs$meta

  # convert to characters
  meta %<>% mutate(
    replicate = as.character(replicate),
    cell_type = as.character(cell_type),
    label = as.character(label)
  )

  # keep only cell types with enough cells
  keep <- meta %>%
    count(cell_type, label) %>%
    group_by(cell_type) %>%
    filter(all(n >= min_cells)) %>%
    pull(cell_type) %>%
    unique()

  # process data into replicate x covariate(s) x cell_type dataframes
  pseudobulks_meta <- keep %>%
    map(~ {
      cell_type <- .
      meta0 <- meta %>% filter(cell_type == !!cell_type)
      # catch cell types without replicates or conditions
      if (n_distinct(meta0$label) < 2) {
        return(NA)
      }
      replicate_counts <- distinct(meta0, label, replicate) %>%
        group_by(label) %>%
        summarise(replicates = n_distinct(replicate)) %>%
        pull(replicates)
      if (any(replicate_counts < min_reps)) {
        return(NA)
      }
      # Summarize metadata by replicate
      meta0 %<>% distinct(across(all_of(covariates)), label, replicate) %>%
        remove_rownames() %>%
        mutate(
          label = as.factor(label),
          group_sample = paste0(replicate, ":", label)
        ) %>%
        select(group_sample, replicate, label, everything()) %>%
        arrange(label)
      return(meta0)
    }) %>%
    setNames(keep)

  # drop NAs
  pseudobulks_meta %<>% extract(!is.na(.))
  return(pseudobulks_meta)
}
