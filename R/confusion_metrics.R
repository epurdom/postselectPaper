#' Annotate DE results with truth and compute performance metrics
#'
#' This function annotates a combined DE results data frame (sim_res_DE)
#' with columns indicating whether each gene-cluster pair is truly DE or non-DE,
#' based on the ground truth (overlap_matrix), and computes basic metrics:
#' true positives (TP), false positives (FP), false discovery proportion (FDP),
#' power, etc.
#'
#' @param sim_res_DE data.frame of DE results, must contain 'gene', 'cluster_id',
#'   and 'adj_p_val' columns.
#' @param overlap_matrix List with components:
#'   \itemize{
#'     \item cell_count_assign_de_gene: Matrix of overlap counts between
#'           assigned clusters and ground truth DE genes
#'     \item assign_size_clust: Named list/vector of cluster sizes
#'     \item gt_size_de_gene: Named list/vector of ground truth DE gene set sizes
#'   }
#' @param cut_off_true Numeric. Overlap threshold to consider a gene DE in a cluster.
#'   If the overlap percentage exceeds this threshold, the gene is marked as truly DE.
#' @param cut_off_false Numeric. Overlap threshold below which a gene is definitely
#'   not DE in a cluster. If the overlap percentage is below this threshold, the gene
#'   is marked as not DE.
#' @param sig_threshold Numeric. Threshold for adjusted p-value significance. Default is 0.1.
#' @param overlap_type Character string. Type of overlap denominator to use:
#'   \itemize{
#'     \item "union": Sum of cluster size and ground truth DE gene size
#'     \item "clust": Cluster size only
#'     \item "gt": Ground truth DE gene size only
#'   }
#'   Default is "union".
#' @param fdr_def Character. How to compute FDR/metrics: \code{"all_hypotheses"} (per
#'   gene-cluster pair) or \code{"gene_level_intersection_3"} (gene-level). Default is \code{"all_hypotheses"}.
#' @param use_indeterminate Logical. If \code{TRUE}, FDP denominator uses all significant
#'   findings (including indeterminate); if \code{FALSE}, only TP + FP. Default is \code{FALSE}.
#'
#' @return List with two components:
#'   \itemize{
#'     \item sim_res_DE: The input data.frame annotated with columns 'is_DE',
#'           'is_not_DE', 'perc_overlap', and 'is_found_sig'.
#'     \item metrics: A named list with performance metrics: TP, FP, FN, FDP, POWER;
#'           for \code{fdr_def = "all_hypotheses"} also includes N (count of non-DE).
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Example usage with simulated data
#' sim_res_DE <- data.frame(
#'   gene = c("GENE1", "GENE2", "GENE3"),
#'   cluster_id = c("C1", "C1", "C2"),
#'   adj_p_val = c(0.01, 0.15, 0.05)
#' )
#' 
#' overlap_matrix <- list(
#'   cell_count_assign_de_gene = matrix(c(50, 10, 5, 45), nrow = 2,
#'                                      dimnames = list(c("GENE1", "GENE2"), c("C1", "C2"))),
#'   assign_size_clust = list(C1 = 100, C2 = 100),
#'   gt_size_de_gene = list(GENE1 = 50, GENE2 = 50)
#' )
#' 
#' result <- confusion_from_overlap_matrix(
#'   sim_res_DE, overlap_matrix,
#'   cut_off_true = 0.1,
#'   cut_off_false = 0.05,
#'   sig_threshold = 0.1,
#'   overlap_type = "union"
#' )
#' print(result$metrics)
#' }
#' 
confusion_from_overlap_matrix <- function(sim_res_DE, overlap_matrix,
                                          cut_off_true,
                                          cut_off_false,
                                          sig_threshold = 0.1,
                                          overlap_type = c("union", "clust", "gt"),
                                          fdr_def = "all_hypotheses",
                                          use_indeterminate = FALSE) {
  overlap_type <- match.arg(overlap_type)
  
  # Initialize annotation columns
  sim_res_DE$is_DE <- FALSE
  sim_res_DE$is_not_DE <- TRUE
  sim_res_DE$perc_overlap <- NA_real_
  
  # Determine entries that are eligible for evaluation
  valid_idx <- sim_res_DE$gene %in% rownames(overlap_matrix$cell_count_assign_de_gene)
  
  if (any(valid_idx)) {
    # Calculate the overlap numerator for each valid gene-cluster pair
    overlap_counts <- mapply(
      function(gene, clust) overlap_matrix$cell_count_assign_de_gene[gene, clust], # size of intersection of DE genes and assigned clusters
      sim_res_DE$gene[valid_idx], sim_res_DE$cluster_id[valid_idx]
    )
    
    # Calculate the denominator for each gene-cluster pair based on overlap_type
    denominator_counts <- mapply(
      function(gene, clust) switch(
        overlap_type,
        "union" = overlap_matrix$assign_size_clust[[clust]] + 
                 overlap_matrix$gt_size_de_gene[[gene]] - 
                 overlap_matrix$cell_count_assign_de_gene[gene, clust], 
        "clust" = overlap_matrix$assign_size_clust[[clust]],
        "gt" = overlap_matrix$gt_size_de_gene[[gene]],
        stop("Invalid overlap type")
      ),
      sim_res_DE$gene[valid_idx], sim_res_DE$cluster_id[valid_idx]
    )
    
    # Compute percentage overlap
    perc_overlaps <- overlap_counts / denominator_counts
    sim_res_DE$perc_overlap[valid_idx] <- perc_overlaps
    
    # Mark DE if overlap exceeds cut_off_true, and not DE if below cut_off_false
    sim_res_DE$is_DE[valid_idx] <- perc_overlaps > cut_off_true
    sim_res_DE$is_not_DE[valid_idx] <- perc_overlaps < cut_off_false
  } else {
    # No valid gene-cluster pairs, perc_overlap will be NA
    sim_res_DE$perc_overlap <- NA_real_
  }
  
  # Mark significant findings by adjusted p-value
  sim_res_DE$is_found_sig <- sim_res_DE$adj_p_val < sig_threshold
  
  if (fdr_def == "all_hypotheses") {
    metrics <- compute_fdr_all_hypotheses(sim_res_DE, use_indeterminate)
  } else if (fdr_def == "gene_level_intersection_3") {
    metrics <- compute_fdr_gene_level_intersection_3(sim_res_DE)
  } else {
    stop("Invalid FDR definition")
  }
  return(list(sim_res_DE = sim_res_DE, metrics = metrics))
}

#' Compute FDR and power using all-hypotheses (gene-cluster) definition
#'
#' @param sim_res_DE Annotated DE results with \code{is_DE}, \code{is_not_DE}, \code{is_found_sig}.
#' @param use_indeterminate If \code{TRUE}, FDP denominator is all significant findings.
#' @return Named list with TP, FP, FN, N, FDP, POWER.
#' @noRd
compute_fdr_all_hypotheses <- function(sim_res_DE, use_indeterminate = FALSE) {
  TP <- sum(sim_res_DE$is_DE & sim_res_DE$is_found_sig, na.rm = TRUE)
  FP <- sum(sim_res_DE$is_not_DE & sim_res_DE$is_found_sig, na.rm = TRUE)
  FN <- sum(sim_res_DE$is_DE & !sim_res_DE$is_found_sig, na.rm = TRUE)
  N <- sum(sim_res_DE$is_not_DE, na.rm = TRUE)

  if (use_indeterminate) {
    denom_found <- sum(sim_res_DE$is_found_sig)
  } else {
    denom_found <- TP + FP
  }
  
  FDP <- if (denom_found > 0) FP / denom_found else 0
  POWER <- if ((TP + FN) > 0) TP / (TP + FN) else 1

  metrics <- list(TP = TP, FP = FP, FN = FN, N = N, FDP = FDP, POWER = POWER)
  return(metrics)
}

#' Compute FDR and power at gene level (significant if found in any cluster; true DE if in any cluster)
#'
#' @param sim_res_DE Annotated DE results with \code{is_DE}, \code{is_found_sig}.
#' @return Named list with TP, FP, FN, FDP, POWER.
#' @noRd
compute_fdr_gene_level_intersection_3 <- function(sim_res_DE) {  
  # find all genes that were found significant at the gene level
  significant_genes <- unique(sim_res_DE$gene[sim_res_DE$is_found_sig])

  # find all genes with at least 1 cluster that is truly DE
  truly_de_genes <- unique(sim_res_DE$gene[sim_res_DE$is_DE])

  TP <- sum(significant_genes %in% truly_de_genes)
  FP <- sum(!(significant_genes %in% truly_de_genes)) # not in truly de genes
  FN <- sum(!(truly_de_genes %in% significant_genes)) # not in significant genes
  FDP <- if ((TP + FP) > 0) FP / (TP + FP) else 0
  POWER <- if ((TP + FN) > 0) TP / (TP + FN) else 1

  metrics <- list(TP = TP, FP = FP, FN = FN, FDP = FDP, POWER = POWER)
  return(metrics)
}

#' Map genes to simulation groups (not_DE or cut_20 / cut_3 / cut_10 by modulo on gene index)
#'
#' Expects DE gene names like \code{simulated_gene-<number>}; group inferred from number %% 3.
#' @param sim_res_DE DE results data frame with \code{gene} column.
#' @param set_de_genes Character vector of ground-truth DE gene names.
#' @return Named character vector: gene -> \code{"not_DE"}, \code{"cut_20"}, \code{"cut_3"}, or \code{"cut_10"}.
#' @noRd
create_gene_group_mapping <- function(sim_res_DE, set_de_genes) {
  not_DE_gene_names <- sim_res_DE$gene[!(sim_res_DE$gene %in% set_de_genes)]
  not_DE_genes <- rep("not_DE", length(not_DE_gene_names))
  names(not_DE_genes) <- not_DE_gene_names
  de_gene_numbers <- as.numeric(sub("simulated_gene-", "", set_de_genes))
  de_set_modulo <- de_gene_numbers %% 3
  de_set_label <- c("cut_20","cut_3", "cut_10")[de_set_modulo + 1]
  names(de_set_label) <- set_de_genes
  gene_group_mapping <- c(not_DE_genes, de_set_label)
  return(gene_group_mapping)
}

#' Compute confusion metrics stratified by gene group (e.g. cut_20, cut_3, cut_10, not_DE)
#'
#' @param sim_res_DE DE results data frame.
#' @param set_de_genes Ground-truth DE gene names.
#' @param overlap_matrix Overlap list for \code{\link{confusion_from_overlap_matrix}}.
#' @param cut_off_true,cut_off_false Overlap thresholds for true/not DE.
#' @param sig_threshold Significance threshold for adjusted p-value.
#' @param overlap_definition Passed as \code{overlap_type}.
#' @param fdr_calc_type Passed as \code{fdr_def}.
#' @param use_R_tilde If \code{TRUE}, do not use indeterminate in FDP (\code{use_indeterminate = FALSE}).
#' @return Named vector of metrics per group (e.g. \code{TP_cut_20}, \code{FDP_not_DE}, ...).
#' @noRd
#' @importFrom stats setNames
get_row_metrics_by_gene_group <- function(sim_res_DE, set_de_genes, overlap_matrix, cut_off_true, cut_off_false, sig_threshold, overlap_definition, fdr_calc_type, use_R_tilde) {
  gene_group_mapping <- create_gene_group_mapping(sim_res_DE, set_de_genes)
  sim_res_DE$gene_group <- gene_group_mapping[sim_res_DE$gene]
  all_sub_group_metrics <- unlist(
    lapply(
      unique(sim_res_DE$gene_group),
      function(ex_group) {
        sim_res_DE_group <- sim_res_DE[sim_res_DE$gene_group == ex_group, ]
        sub_group_confusion <- confusion_from_overlap_matrix(
          sim_res_DE_group,
          overlap_matrix,
          cut_off_true = cut_off_true,
          cut_off_false = cut_off_false,
          sig_threshold = sig_threshold,
          overlap_type = overlap_definition,
          fdr_def = fdr_calc_type,
          use_indeterminate = !use_R_tilde
        )
        setNames(sub_group_confusion$metrics,
                 paste(names(sub_group_confusion$metrics), ex_group, sep = "_"))
      }
    ),
    recursive = FALSE
  )
  return(all_sub_group_metrics)
}

#' FDR/metrics by gene group for augSimAnalyse31020 hierarchical simulation
#'
#' @inheritParams get_row_metrics_by_gene_group
#' @param sims_description List with \code{set_de_genes}.
#' @return Named vector of metrics per gene group.
#' @noRd
compute_fdr_augSimAnalyse31020_hierarchical <- function(sim_res_DE, sims_description, overlap_matrix, cut_off_true, cut_off_false, sig_threshold, overlap_definition, fdr_calc_type, use_R_tilde) {
  set_de_genes <- sims_description$set_de_genes
  all_sub_group_metrics <- get_row_metrics_by_gene_group(sim_res_DE, set_de_genes, overlap_matrix, cut_off_true, cut_off_false, sig_threshold, overlap_definition, fdr_calc_type, use_R_tilde)
  return(all_sub_group_metrics)
}

#' Map genes to binary state: DE vs not_DE
#'
#' @param sim_res_DE DE results with \code{gene} column.
#' @param set_de_genes Ground-truth DE gene names.
#' @return Named character vector: gene -> \code{"not_DE"} or \code{"DE"}.
#' @noRd
create_gene_state_mapping <- function(sim_res_DE, set_de_genes) {
  not_DE_gene_names <- unique(sim_res_DE$gene[!(sim_res_DE$gene %in% set_de_genes)])
  not_DE_genes <- rep("not_DE", length(not_DE_gene_names))
  names(not_DE_genes) <- not_DE_gene_names
  de_set_label <- rep("DE", length(set_de_genes))
  names(de_set_label) <- set_de_genes

  gene_state_mapping <- c(not_DE_genes, de_set_label)
  return(gene_state_mapping)
}

#' Compute confusion metrics stratified by gene DE state (DE vs not_DE)
#'
#' @inheritParams get_row_metrics_by_gene_group
#' @return Named vector of metrics per state (e.g. \code{TP_DE}, \code{FDP_not_DE}, ...).
#' @noRd
get_row_metrics_by_gene_DE_state <- function(sim_res_DE, set_de_genes, overlap_matrix, cut_off_true, cut_off_false, sig_threshold, overlap_definition, fdr_calc_type, use_R_tilde) {
  gene_state_mapping <- create_gene_state_mapping(sim_res_DE, set_de_genes)
  sim_res_DE$gene_state <- gene_state_mapping[sim_res_DE$gene]
  all_sub_group_metrics <- unlist(
    lapply(
      unique(sim_res_DE$gene_state),
      function(ex_state) {
        sim_res_DE_state <- sim_res_DE[sim_res_DE$gene_state == ex_state, ]
        sub_group_confusion <- confusion_from_overlap_matrix(
          sim_res_DE_state,
          overlap_matrix,
          cut_off_true = cut_off_true,
          cut_off_false = cut_off_false,
          sig_threshold = sig_threshold,
          overlap_type = overlap_definition,
          fdr_def = fdr_calc_type,
          use_indeterminate = !use_R_tilde
        )
        setNames(sub_group_confusion$metrics,
                 paste(names(sub_group_confusion$metrics), ex_state, sep = "_"))
      }
    ),
    recursive = FALSE
  )
  return(all_sub_group_metrics)
}

#' FDR/metrics by gene DE state for augSimLeiden1 simulation
#'
#' @inheritParams compute_fdr_augSimAnalyse31020_hierarchical
#' @return Named vector of metrics per DE state.
#' @noRd
compute_fdr_augSimLeiden1 <- function(sim_res_DE, sims_description, overlap_matrix, cut_off_true, cut_off_false, sig_threshold, overlap_definition, fdr_calc_type, use_R_tilde) {
  set_de_genes <- sims_description$set_de_genes
  all_sub_group_metrics <- get_row_metrics_by_gene_DE_state(sim_res_DE, set_de_genes, overlap_matrix, cut_off_true, cut_off_false, sig_threshold, overlap_definition, fdr_calc_type, use_R_tilde)
  return(all_sub_group_metrics)
}

#' Dispatch FDR/metrics computation by simulation subcategory
#'
#' @inheritParams compute_fdr_augSimAnalyse31020_hierarchical
#' @param subcategory Character. \code{"augSimAnalyse31020_hierarchical"}, \code{"augSimLeiden1"}, or \code{"None"} (invalid).
#' @return Named vector of metrics (contents depend on subcategory).
#' @noRd
compute_fdr_by_subcategory <- function(sim_res_DE, sims_description, overlap_matrix, cut_off_true, cut_off_false, sig_threshold, overlap_definition, fdr_calc_type, use_R_tilde, subcategory = "None") {
  if (subcategory == "augSimAnalyse31020_hierarchical") {
    return(compute_fdr_augSimAnalyse31020_hierarchical(sim_res_DE, sims_description, overlap_matrix, cut_off_true, cut_off_false, sig_threshold, overlap_definition, fdr_calc_type, use_R_tilde))
  } else if (subcategory == "augSimLeiden1") {
    return(compute_fdr_augSimLeiden1(sim_res_DE, sims_description, overlap_matrix, cut_off_true, cut_off_false, sig_threshold, overlap_definition, fdr_calc_type, use_R_tilde))
  } else {
    stop("Invalid subcategory")
  }
}

