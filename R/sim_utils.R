#' Print quick FDP/TPP control results for an analysis run
#'
#' Runs cluster-by-cluster overlap control via \code{\link{run_clust_by_clust_w_overlap}}
#' and prints the selected control summary (local, global, or clusters-global).
#'
#' @param anls_res List with \code{res_DE} (per-cluster DE results) and \code{de_overlap_info}.
#' @param cut_off_true,cut_off_false Numeric. Overlap cut-offs for DE vs non-DE gene sets.
#' @param sig_threshold Numeric. Significance threshold for DE calls.
#' @param overlap_type Character. \code{"union"}, \code{"clust"}, or \code{"gt"}.
#' @param ctrl_type Character. Which control to print: \code{"ctrl_clusts_loc"}, \code{"ctrl_glob"}, or \code{"ctrl_clusts_glob"} (default \code{"ctrl_glob"}).
#'
#' @return Invisibly returns the control table; used for side effect (printing).
#'
#' @export
print_quick_glob_ctrl <- function(anls_res, cut_off_true, cut_off_false, sig_threshold, overlap_type, ctrl_type = "ctrl_glob") {
  ctrl_res_curr_anls <- run_clust_by_clust_w_overlap(anls_res$res_DE, cut_off_true, cut_off_false, sig_threshold, overlap_type, anls_res$de_overlap_info)
  print(paste0("Control type: ", ctrl_type))
  to_print <- ctrl_res_curr_anls[[ctrl_type]]
  print("Ctrl results:")
  print(to_print)
  invisible(to_print)
}

#' Create UMAP comparison plots by treatment (old vs new coordinates)
#'
#' Builds a 2-row grid: row 1 = original UMAP (all + per treatment), row 2 = new
#' UMAP. Saves the combined plot to \code{save_dir/save_name}.
#'
#' @param umap_kmeans_df Data frame with columns \code{UMAP1_old}, \code{UMAP2_old}, \code{UMAP1}, \code{UMAP2}, \code{treatment}.
#' @param category_var Unused; kept for interface compatibility.
#' @param category_values Character vector. Treatment values to plot in separate panels.
#' @param save_dir Character. Directory path for saving the plot.
#' @param save_name Character. Filename for the saved plot.
#'
#' @return Invisible; saves the plot via \code{ggplot2::ggsave}.
#'
#' @importFrom ggplot2 ggplot aes geom_point theme_classic labs ggsave
#' @importFrom patchwork wrap_plots
#' @export
create_treatment_umap_comparison <- function(umap_kmeans_df, category_var, category_values, save_dir, save_name) {
  plot_list <- list()

  plot_list[[1]] <- ggplot2::ggplot(umap_kmeans_df, ggplot2::aes(x = UMAP1_old, y = UMAP2_old, color = treatment)) +
    ggplot2::geom_point(size = 0.8, alpha = 0.7) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = "All Treatments (Original UMAP)")

  for (i in seq_along(category_values)) {
    treatment_value <- category_values[i]
    data_subset <- umap_kmeans_df[umap_kmeans_df$treatment == treatment_value, ]
    plot_list[[i + 1]] <- ggplot2::ggplot(data_subset, ggplot2::aes(x = UMAP1_old, y = UMAP2_old, color = treatment)) +
      ggplot2::geom_point(size = 0.8, alpha = 0.7) +
      ggplot2::theme_classic() +
      ggplot2::labs(title = paste("Treatment =", treatment_value, "(Original UMAP)"))
  }

  plot_list[[length(category_values) + 2]] <- ggplot2::ggplot(umap_kmeans_df, ggplot2::aes(x = UMAP1, y = UMAP2, color = treatment)) +
    ggplot2::geom_point(size = 0.8, alpha = 0.7) +
    ggplot2::theme_classic() +
    ggplot2::labs(title = "All Treatments (New UMAP)")

  for (i in seq_along(category_values)) {
    treatment_value <- category_values[i]
    data_subset <- umap_kmeans_df[umap_kmeans_df$treatment == treatment_value, ]
    plot_list[[i + length(category_values) + 2]] <- ggplot2::ggplot(data_subset, ggplot2::aes(x = UMAP1, y = UMAP2, color = treatment)) +
      ggplot2::geom_point(size = 0.8, alpha = 0.7) +
      ggplot2::theme_classic() +
      ggplot2::labs(title = paste("Treatment =", treatment_value, "(New UMAP)"))
  }

  n_plots_per_row <- length(category_values) + 1
  umap_plots <- patchwork::wrap_plots(plot_list, ncol = n_plots_per_row, nrow = 2)
  ggplot2::ggsave(umap_plots, filename = file.path(save_dir, save_name), width = 5 * n_plots_per_row, height = 10)
  invisible(NULL)
}

#' Per-cluster and overall DE analysis metrics (TP/FP/FN, FDR, TPR)
#'
#' Computes per-cluster and overall (Total, Union) counts of found/ground-truth
#' genes, true/false positives/negatives, and masked/unmasked FDR and TPR.
#'
#' @param de_found_sets List. Per-cluster character vectors of DE gene names found.
#' @param de_gt_sets List. Per-cluster character vectors of ground-truth DE genes.
#' @param clust_set List. Per-cluster character vectors of cell IDs (used for structure only).
#' @param sim_gene_list Character. Simulated (DE) gene names.
#' @param non_sim_gene_list Character. Non-simulated gene names.
#' @param fn Optional. File path to save the result matrix as RDS; if \code{NULL}, not saved.
#'
#' @return Numeric matrix. Rows = clusters plus \code{Total} and \code{Union}; columns include \code{curr_tp}, \code{curr_fp}, \code{curr_fn}, \code{fdr_masked}, \code{tpr_masked}, etc.
#'
#' @export
create_per_clust_res_de <- function(de_found_sets, de_gt_sets, clust_set, sim_gene_list, non_sim_gene_list, fn = NULL) {
  print("Calculating per-cluster TP/FP/FN results")
  per_clust_results <- lapply(names(clust_set), function(clust_name) {
    curr_found_genes <- de_found_sets[[clust_name]]
    curr_gt_genes <- de_gt_sets[[clust_name]]
    curr_found_sim_genes <- intersect(curr_found_genes, sim_gene_list)
    curr_found_non_sim_genes <- intersect(curr_found_genes, non_sim_gene_list)
    curr_tp <- intersect(curr_found_genes, curr_gt_genes)
    curr_fp <- setdiff(curr_found_genes, curr_gt_genes)
    curr_fn <- setdiff(curr_gt_genes, curr_found_genes)
    list(curr_found_genes = curr_found_genes, curr_gt_genes = curr_gt_genes, curr_found_sim_genes = curr_found_sim_genes,
         curr_found_non_sim_genes = curr_found_non_sim_genes, curr_tp = curr_tp, curr_fp = curr_fp, curr_fn = curr_fn)
  })
  names(per_clust_results) <- names(clust_set)

  # convert from sets to lengths of those sets
  per_clust_results_lengths <- lapply(per_clust_results, function(x) {
    list(curr_found_genes = length(x$curr_found_genes), curr_gt_genes = length(x$curr_gt_genes),
         curr_found_sim_genes = length(x$curr_found_sim_genes),
         curr_found_non_sim_genes = length(x$curr_found_non_sim_genes),
         curr_tp = length(x$curr_tp), curr_fp = length(x$curr_fp), curr_fn = length(x$curr_fn))
  })
  # convert to a numeric matrix
  per_clust_results_lengths_matrix <- do.call(rbind, lapply(per_clust_results_lengths, unlist))

  # add a total column row to the matrix
  per_clust_results_lengths_matrix <- rbind(per_clust_results_lengths_matrix, 
                                           Total = colSums(per_clust_results_lengths_matrix))

  print("Calculating overall results")
  # Initialize empty lists for overall results
  overall_results <- list()
  for (result_name in names(per_clust_results[[1]])) {
    overall_results[[result_name]] <- character(0)
  }
  # Combine results across all clusters using union
  for (clust_name in names(per_clust_results)) {
    curr_results <- per_clust_results[[clust_name]]
    for (result_name in names(overall_results)) {
      overall_results[[result_name]] <- c(overall_results[[result_name]], curr_results[[result_name]])
    }
  }

  # union results
  overall_results_union <- list()
  for (result_name in names(overall_results)) {
    overall_results_union[[result_name]] <- length(unique(overall_results[[result_name]]))
  }
  overall_results_union$curr_found_sim_genes <- length(intersect(unique(overall_results$curr_found_genes), sim_gene_list))
  overall_results_union$curr_found_non_sim_genes <- length(setdiff(unique(overall_results$curr_found_genes), sim_gene_list))
  per_clust_results_lengths_matrix <- rbind(per_clust_results_lengths_matrix, 
                                           Union = unlist(overall_results_union))

  # add new columns to per_clust_results_lengths_matrix for fdr metrics
  per_clust_results_lengths_matrix <- cbind(per_clust_results_lengths_matrix,
                                          fdr_masked = per_clust_results_lengths_matrix[,"curr_fp"] / 
                                            (per_clust_results_lengths_matrix[,"curr_tp"] + per_clust_results_lengths_matrix[,"curr_fp"]),
                                          fdr_unmasked = per_clust_results_lengths_matrix[,"curr_found_non_sim_genes"] / 
                                            per_clust_results_lengths_matrix[,"curr_found_genes"],
                                          tpr_masked = per_clust_results_lengths_matrix[,"curr_tp"] / 
                                            per_clust_results_lengths_matrix[,"curr_gt_genes"], 
                                          tpr_unmasked = per_clust_results_lengths_matrix[,"curr_found_sim_genes"] / 
                                            per_clust_results_lengths_matrix[,"curr_found_genes"])

  if (!is.null(fn)) {
    saveRDS(per_clust_results_lengths_matrix, fn)
  }
  return(per_clust_results_lengths_matrix)
}

#' Build FDP/TPP summary table for LaTeX or display
#'
#' Formats per-cluster DE counts into strings (e.g. \code{"x/y"} for TP/FP, \code{"a% / b%"} for FDP/TPP).
#'
#' @param de_table_res Matrix or data frame from \code{\link{create_per_clust_res_de}} with columns \code{curr_tp}, \code{curr_fp}, \code{curr_fn}, \code{curr_found_sim_genes}, \code{curr_found_non_sim_genes}.
#'
#' @return Data frame with same row names and columns \code{TP}, \code{FP}, \code{TPP}, \code{FDP} (formatted strings).
#'
#' @export
create_fdp_tpp_table <- function(de_table_res) {
  ratio_lists <- list()
  ratio_lists[["TP"]] <- list("curr_tp", "curr_found_sim_genes")
  ratio_lists[["FP"]] <- list("curr_fp", "curr_found_non_sim_genes")

  percentage_list <- list()
  percentage_list[["TPP"]] <- list("curr_tp", "curr_fn", "curr_found_sim_genes", "curr_fn")
  percentage_list[["FDP"]] <- list("curr_fp", "curr_tp", "curr_found_non_sim_genes", "curr_found_sim_genes")
  # creat strings for latex table

  de_table_res_str <- data.frame(matrix(nrow=nrow(de_table_res), ncol=length(ratio_lists) + length(percentage_list)))
  colnames(de_table_res_str) <- c(names(ratio_lists), names(percentage_list))
  rownames(de_table_res_str) <- rownames(de_table_res)
  for (clust_name in rownames(de_table_res)) {
    for (ratio_name in names(ratio_lists)) {
      de_table_res_str[clust_name, ratio_name] <- sprintf("%.0f/%.0f", 
                                              de_table_res[clust_name, ratio_lists[[ratio_name]][[1]]],
                                              de_table_res[clust_name, ratio_lists[[ratio_name]][[2]]])
    }
    for (percentage_name in names(percentage_list)) { 
      fdp <- 100 * de_table_res[clust_name, percentage_list[[percentage_name]][[1]]] / (de_table_res[clust_name, percentage_list[[percentage_name]][[1]]] + de_table_res[clust_name, percentage_list[[percentage_name]][[2]]])
      tpp <- 100 * de_table_res[clust_name, percentage_list[[percentage_name]][[3]]] / (de_table_res[clust_name, percentage_list[[percentage_name]][[3]]] + de_table_res[clust_name, percentage_list[[percentage_name]][[4]]])
      if (is.nan(fdp)) {
        fdp <- 0
      }
      if (is.nan(tpp)) {
        tpp <- 0
      }
      de_table_res_str[clust_name, percentage_name] <- sprintf("%.1f%% / %.1f%%", fdp, tpp)
    }
  }
  return(de_table_res_str)
}




#' Run pseudobulk differential expression for a clustering
#'
#' Copies cluster/sample/group from \code{colData} to \code{cluster_id}, \code{sample_id}, \code{group_id}, aggregates with \code{muscat::aggregateData}, then runs \code{muscat::pbDS}.
#'
#' @param sce \code{SingleCellExperiment} with \code{counts} and \code{colData} containing cluster, sample, and group columns.
#' @param leiden_clusters Unused; cluster assignment is taken from \code{sce[[cluster_id_col]]}.
#' @param cluster_id_col Character. \code{colData} column for cluster.
#' @param sample_id_col Character. \code{colData} column for sample.
#' @param group_id_col Character. \code{colData} column for condition/group.
#' @param method \code{"DESeq2"} or \code{"edgeR"} (default \code{"DESeq2"}).
#'
#' @return List of data frames. Per-cluster DE tables from \code{muscat::pbDS} (i.e. \code{res_msct$table}); each element has \code{gene}, \code{p_adj.glb}, etc.
#'
#' @importFrom muscat aggregateData pbDS
#' @importFrom S4Vectors metadata
#' @export
run_de_analysis <- function(sce, leiden_clusters, cluster_id_col, sample_id_col, group_id_col, method = "DESeq2") {
  # Assign cluster, sample and group IDs
  print("Assigning cluster, sample and group IDs")
  sce$cluster_id <- sce[[cluster_id_col]]
  sce$sample_id <- sce[[sample_id_col]]
  sce$group_id <- sce[[group_id_col]]

  # Make sure group_id is a factor with proper levels
  sce$group_id <- factor(sce$group_id)
  sce$cluster_id <- factor(sce$cluster_id)
  sce$cell_type <- factor(sce$cell_type)

  # aggregate into pseudobulk
  print("Aggregating into pseudobulk")
  pseudobulk <- muscat::aggregateData(
    sce,
    assay = "counts",
    fun = "sum",
    by = c("cluster_id", "sample_id")
  )

  # Create experiment_info data frame
  print("Creating experiment_info data frame")
  sample_info <- unique(data.frame(
    sample_id = sce$sample_id,
    group_id = sce$group_id,
    stringsAsFactors = TRUE  # Ensure these are factors
  ))
  S4Vectors::metadata(pseudobulk)$experiment_info <- sample_info

  # Perform differential expression analysis
  print("Performing differential expression analysis")
  res_msct <- muscat::pbDS(
    pseudobulk, 
    method = method,
    filter = "both",
    verbose = TRUE
  )
  print("Returning DE results")
  return(res_msct$table)
}

#' Build logical matrix of DE calls per cluster per gene
#'
#' @param de_by_clust List. Per-cluster DE result tables; each element has \code{p_adj.glb} and \code{gene}.
#' @param gene_list Character. Gene names (rows of output matrix).
#' @param sig_threshold Numeric. Genes with \code{p_adj.glb < sig_threshold} are considered DE.
#' @param clust_names Character. Cluster names (columns of output).
#'
#' @return Logical matrix. Rows = genes (from \code{gene_list}), columns = clusters; \code{TRUE} = DE in that cluster.
#'
#' @export
get_de_matrix <- function(de_by_clust, gene_list, sig_threshold, clust_names) {
  de_matrix <- list()
  for (clust_name in clust_names) {
    if (!(clust_name %in% names(de_by_clust))) {
      print(paste("Cluster", clust_name, "not found in de_by_clust"))
      de_matrix[[clust_name]] <- rep(FALSE, length(gene_list))
    } else {
      de_clust <- de_by_clust[[clust_name]]
      de_genes <- de_clust[de_clust$p_adj.glb < sig_threshold,c('gene')]
      de_matrix[[clust_name]] <- gene_list %in% de_genes
    }
  }
  # convert to dataframe
  de_matrix <- do.call(rbind, de_matrix)
  return(de_matrix)
}

#' Empirical DE gene sets per cluster
#'
#' @param de_by_clust List. Per-cluster DE result tables with \code{p_adj.glb} and \code{gene}.
#' @param gene_list Character. Full gene list (used when a cluster is missing).
#' @param sig_threshold Numeric. Significance cutoff (genes with \code{p_adj.glb < sig_threshold}).
#' @param clust_names Character. Cluster names for which to extract DE sets.
#'
#' @return List. One element per cluster; each element is a character vector of DE gene names.
#'
#' @export
get_de_sets <- function(de_by_clust, gene_list, sig_threshold, clust_names) {
  de_sets <- list()
  for (clust_name in clust_names) {
    clust_name <- as.character(clust_name)
    if (!(clust_name %in% names(de_by_clust))) {
      print(paste("Cluster", clust_name, "not found in de_by_clust"))
      de_sets[[clust_name]] <- list()
    } else {
      de_clust <- de_by_clust[[clust_name]]
      de_genes <- de_clust[de_clust$p_adj.glb < sig_threshold,c('gene')]
      de_sets[[clust_name]] <- de_genes
    }
  }
  # convert to dataframe
  return(de_sets)
}

#' Ground-truth DE indicator matrices from cluster-neighborhood overlap
#'
#' Computes cluster-by-gene overlap proportions and returns binary matrices for
#' DE (proportion >= \code{cut_off_prop_true}) and non-DE (proportion < \code{cut_off_prop_false}).
#'
#' @param clust_mat Matrix. Cluster membership (rows = clusters, columns = cells).
#' @param nbhd_matrix Matrix. DE neighborhood; dimensions must align so \code{clust_mat %*% nbhd_matrix} is clusters x genes.
#' @param cut_off_prop_true,cut_off_prop_false Numeric. Proportion thresholds for DE and non-DE.
#' @param overlap_type Character. Denominator for proportion: \code{"clust"}, \code{"nbhd"}, or \code{"union"} (default \code{"clust"}).
#'
#' @return List with \code{de_gt_matrix}, \code{non_de_gt_matrix}, \code{prop_final_matrix}.
#'
#' @importFrom Matrix rowSums colSums
#' @export
get_de_gt_matrix <- function(clust_mat, nbhd_matrix, cut_off_prop_true, cut_off_prop_false, overlap_type = "clust") {
  inter_matrix <- clust_mat %*% nbhd_matrix
  clust_size_vec <- Matrix::rowSums(clust_mat)
  nbhd_size_vec <- Matrix::colSums(nbhd_matrix)
  clust_size_matrix <- matrix(rep(clust_size_vec, length(nbhd_size_vec)), nrow = length(clust_size_vec))
  nbhd_size_matrix <- t(matrix(rep(nbhd_size_vec, length(clust_size_vec)), nrow = length(nbhd_size_vec)))
  union_matrix <- clust_size_matrix + nbhd_size_matrix - inter_matrix
  if (overlap_type == "clust") {
    prop_final_matrix <- inter_matrix / clust_size_matrix
  } else if (overlap_type == "nbhd") {
    prop_final_matrix <- inter_matrix / nbhd_size_matrix
  } else if (overlap_type == "union") {
    prop_final_matrix <- inter_matrix / union_matrix
  }
  de_gt_matrix <- matrix(0, nrow = nrow(prop_final_matrix), ncol = ncol(prop_final_matrix))
  de_gt_matrix[prop_final_matrix >= cut_off_prop_true] <- 1
  non_de_gt_matrix <- matrix(0, nrow = nrow(prop_final_matrix), ncol = ncol(prop_final_matrix))
  non_de_gt_matrix[prop_final_matrix < cut_off_prop_false] <- 1
  return(list(de_gt_matrix = de_gt_matrix, non_de_gt_matrix = non_de_gt_matrix, prop_final_matrix = prop_final_matrix))
}

#' Ground-truth DE and non-DE gene sets per cluster from overlap proportions
#'
#' Same overlap logic as \code{\link{get_de_gt_matrix}} but returns per-cluster
#' character vectors of gene names instead of matrices.
#'
#' @param clust_mat Matrix. Cluster membership (rows = clusters, columns = cells).
#' @param nbhd_matrix Matrix. DE neighborhood (dimensions compatible with \code{clust_mat %*% nbhd_matrix}).
#' @param cut_off_prop_true,cut_off_prop_false Numeric. Proportion cut-offs for DE and non-DE.
#' @param gene_list Character. Gene names (one per column of the cluster-by-gene proportion matrix).
#' @param overlap_type Character. \code{"clust"}, \code{"nbhd"}, or \code{"union"} (default \code{"clust"}).
#'
#' @return List with \code{de_genes_list}, \code{non_de_genes_list}, \code{prop_final_matrix}.
#'
#' @importFrom Matrix rowSums colSums
#' @export
get_de_gt_set <- function(clust_mat, nbhd_matrix, cut_off_prop_true, cut_off_prop_false, gene_list, overlap_type = "clust") {
  inter_matrix <- clust_mat %*% nbhd_matrix
  clust_size_vec <- Matrix::rowSums(clust_mat)
  nbhd_size_vec <- Matrix::colSums(nbhd_matrix)
  clust_size_matrix <- base::sweep(matrix(1, nrow=length(clust_size_vec), ncol=length(nbhd_size_vec)), 1, clust_size_vec, "*")
  nbhd_size_matrix <- base::sweep(matrix(1, nrow=length(clust_size_vec), ncol=length(nbhd_size_vec)), 2, nbhd_size_vec, "*")
  union_matrix <- clust_size_matrix + nbhd_size_matrix - inter_matrix
  if (overlap_type == "clust") {
    prop_final_matrix <- inter_matrix / clust_size_matrix
  } else if (overlap_type == "nbhd") {
    prop_final_matrix <- inter_matrix / nbhd_size_matrix
  } else if (overlap_type == "union") {
    prop_final_matrix <- inter_matrix / union_matrix
  }
  de_genes_list <- list()
  for (i in seq_len(nrow(prop_final_matrix))) {
    de_genes_list[[i]] <- gene_list[prop_final_matrix[i, ] >= cut_off_prop_true]
  }
  non_de_genes_list <- list()
  for (i in seq_len(nrow(prop_final_matrix))) {
    non_de_genes_list[[i]] <- gene_list[prop_final_matrix[i, ] < cut_off_prop_false]
  }
  return(list(de_genes_list = de_genes_list, non_de_genes_list = non_de_genes_list, prop_final_matrix = prop_final_matrix))
}

#' Per-cluster DE summary: proportion of significant genes that are simulated
#'
#' @param de_by_clust List. Per-cluster DE result tables with \code{p_adj.glb} and \code{gene}.
#' @param sig_threshold Numeric. FDR threshold for significance.
#'
#' @return Data frame. Rows = clusters; columns \code{clust_name}, \code{prop_sim_fdr}.
#'
#' @export
get_de_features <- function(de_by_clust, sig_threshold) {
  de_features <- list()
  for (clust_name in names(de_by_clust)) {
    de_clust <- de_by_clust[[clust_name]]
    de_genes_fdr <- de_clust[de_clust$p_adj.glb < sig_threshold, ]
    prop_sim_fdr <- 1-sum(grepl("simulated", de_genes_fdr$gene)) / nrow(de_genes_fdr)
    if (is.numeric(clust_name)) {
      print(paste("Cluster", clust_name, "is numeric"))
      clust_name <- paste0("Cluster ", clust_name)
    }
    de_features[[clust_name]] <- list(clust_name = clust_name, prop_sim_fdr = prop_sim_fdr)
  }
  # convert to dataframe
  de_features <- do.call(rbind, de_features)
  return(de_features)
}

#' Full DE pipeline from SCE and clustering
#'
#' Runs \code{\link{run_de_analysis}}, builds found and ground-truth DE sets via
#' \code{\link{get_de_sets}} and \code{\link{get_de_gt_set}}, then
#' \code{\link{create_per_clust_res_de}} and \code{\link{create_fdp_tpp_table}}.
#' \code{curr_sce} must have \code{rowData} with \code{name}, \code{is_simulated}, and (for simulated genes) \code{is_de_cell}; \code{colData} with \code{sample}, \code{fake_condition}.
#'
#' @param curr_sce \code{SingleCellExperiment} with \code{counts}, \code{rowData}, \code{colData} as above.
#' @param leiden_assignment_list Factor or character. Cluster assignment per cell.
#' @param sig_threshold Numeric. Significance threshold for DE.
#' @param overlap_type Character. Overlap definition for ground truth: \code{"clust"}, \code{"nbhd"}, or \code{"union"}.
#' @param cut_off_prop_true,cut_off_prop_false Numeric. Proportion cut-offs for DE vs non-DE sets.
#' @param method \code{"DESeq2"} or \code{"edgeR"} (default \code{"DESeq2"}).
#'
#' @return List with \code{full_analysis}, \code{de_table_res}, \code{de_table_res_str}, \code{leiden_de_found_sets}, \code{leiden_de_gt_set}, \code{leiden_non_de_gt_set}, \code{leiden_prop_final_matrix}.
#'
#' @importFrom SummarizedExperiment rowData colData
#' @export
run_de_pipeline_from_sce_grouping <- function(curr_sce, leiden_assignment_list, sig_threshold, overlap_type, cut_off_prop_true, cut_off_prop_false, method = "DESeq2") {
  print("Prepping for running the de analysis")
  curr_gene_list <- SummarizedExperiment::rowData(curr_sce)$name
  is_simulated_gene <- SummarizedExperiment::rowData(curr_sce)$is_simulated
  curr_sim_genes <- curr_gene_list[is_simulated_gene]
  curr_non_sim_gene_list <- curr_gene_list[!is_simulated_gene]
  leiden_assign_name <- paste0("leiden_assignment")
  SummarizedExperiment::colData(curr_sce)[[leiden_assign_name]] <- leiden_assignment_list

  print("Running the de analysis")
  leiden_de_analysis <- run_de_analysis(curr_sce, leiden_assignment_list, leiden_assign_name, "sample", "fake_condition", method = method)
  clust_names <- unique(leiden_assignment_list)
  leiden_de_found_sets <- get_de_sets(leiden_de_analysis, curr_gene_list, sig_threshold, clust_names = clust_names)
  
  # prep for getting ground truth de sets
  print("Prepping for getting ground truth de sets")
  leiden_assignment_matrix <- t(model.matrix(~leiden_assignment_list-1))
  rownames(leiden_assignment_matrix) <- gsub("leiden_assignment_list", "", rownames(leiden_assignment_matrix))
  leiden_cluster_sets <- lapply(rownames(leiden_assignment_matrix), function(clust_name) {
    colnames(leiden_assignment_matrix)[leiden_assignment_matrix[clust_name,] == 1]
  })
  names(leiden_cluster_sets) <- rownames(leiden_assignment_matrix)
  
  print("Getting nbhd matrix")
  if (length(curr_sim_genes) > 0) {
    nbhd_sublist <- SummarizedExperiment::rowData(curr_sce)$is_de_cell[curr_gene_list %in% curr_sim_genes]
    nbhd_matrix <- t(do.call(rbind, lapply(nbhd_sublist, function(x) as.numeric(x)))) # transpose to match the format of the other matrices

    # get ground truth de sets
    print("Getting ground truth de sets")
    leiden_de_gt_res <- get_de_gt_set(leiden_assignment_matrix, nbhd_matrix, cut_off_prop_true, cut_off_prop_false, curr_sim_genes, overlap_type = overlap_type)
    leiden_de_gt_set= leiden_de_gt_res[[1]] 
    names(leiden_de_gt_set) <- rownames(leiden_assignment_matrix)
    leiden_non_de_gt_set= leiden_de_gt_res[[2]] 
    names(leiden_non_de_gt_set) <- rownames(leiden_assignment_matrix)
    leiden_prop_final_matrix= leiden_de_gt_res[[3]]
    names(leiden_prop_final_matrix) <- rownames(leiden_assignment_matrix)
  } else {
    # make empty lists for all the de gt sets
    leiden_de_gt_set <- lapply(rownames(leiden_assignment_matrix), function(clust_name) list())
    names(leiden_de_gt_set) <- rownames(leiden_assignment_matrix)
    leiden_non_de_gt_set <- lapply(rownames(leiden_assignment_matrix), function(clust_name) list())
    names(leiden_non_de_gt_set) <- rownames(leiden_assignment_matrix)
    leiden_prop_final_matrix <- lapply(rownames(leiden_assignment_matrix), function(clust_name) list())
    names(leiden_prop_final_matrix) <- rownames(leiden_assignment_matrix)
  }

  ## now create per cluster de results
  print("Creating per cluster de results")
  fn <- NULL #paste0("2025-scRNA-postsel-DE/lemur_sim/results/celltype_bias/vary_leid_res/pclus_", res, "_", overlap_type, ".rds")
  de_table_res <- create_per_clust_res_de(leiden_de_found_sets, leiden_de_gt_set, leiden_cluster_sets, curr_sim_genes, curr_non_sim_gene_list, fn = fn)

  ## now create fdp tpp table of results for latex
  print("Creating fdp tpp table of results for latex")
  de_table_res_str <- create_fdp_tpp_table(de_table_res)  
  return(list(full_analysis = leiden_de_analysis, de_table_res = de_table_res, de_table_res_str = de_table_res_str, leiden_de_found_sets = leiden_de_found_sets, leiden_de_gt_set = leiden_de_gt_set, leiden_non_de_gt_set = leiden_non_de_gt_set, leiden_prop_final_matrix = leiden_prop_final_matrix))
}

#' Compute and save UMAP coordinates from SCE (PCA or Harmony)
#'
#' Subsets to non-simulated genes, builds a Seurat object, runs PCA or Harmony
#' (and PCA before Harmony if needed), then UMAP. Saves UMAP coordinates to
#' \code{dir_analysis_dump}.
#'
#' @param sce \code{SingleCellExperiment} with \code{rowData(sce)$is_simulated}.
#' @param dir_analysis_dump Character. Directory to save UMAP RDS.
#' @param reduction_type Character. \code{"pca"} or \code{"harmony"} (default \code{"pca"}).
#' @param num_pcs Integer. Number of PCs / UMAP dimensions (default \code{50}).
#' @param file_prefixes Character. Prefix for saved filename (default \code{""}). If \code{NULL}, uses \code{file_prefixes} from the calling environment.
#'
#' @return Invisible; saves \code{<file_prefixes>umap_coords_data_<reduction_type>.rds} in \code{dir_analysis_dump}.
#'
#' @importFrom Seurat CreateSeuratObject AddMetaData NormalizeData FindVariableFeatures ScaleData RunPCA RunUMAP Embeddings
#' @importFrom harmony RunHarmony
#' @importFrom SingleCellExperiment counts
#' @importFrom SummarizedExperiment rowData colData
#' @export
get_data_umap_from_sce <- function(sce, dir_analysis_dump, reduction_type = "pca", num_pcs = 50, file_prefixes = "") {
  if (is.null(file_prefixes)) file_prefixes <- get("file_prefixes", envir = parent.frame())
  data_sce <- sce[SummarizedExperiment::rowData(sce)$is_simulated == FALSE, ]
  print("Preprocessing seurat object")
  seurat_obj_data <- Seurat::CreateSeuratObject(counts = SingleCellExperiment::counts(data_sce))
  seurat_obj_data <- Seurat::AddMetaData(seurat_obj_data, metadata = data.frame(SummarizedExperiment::colData(data_sce)))
  seurat_obj_data <- Seurat::NormalizeData(seurat_obj_data)
  seurat_obj_data <- Seurat::FindVariableFeatures(seurat_obj_data)
  seurat_obj_data <- Seurat::ScaleData(seurat_obj_data)

  print("Running seurat pca and clustering")
  if (reduction_type == "pca") {
    seurat_obj_data <- Seurat::RunPCA(seurat_obj_data, npcs = num_pcs)
  } else if (reduction_type == "harmony") {
    seurat_obj_data <- Seurat::RunPCA(seurat_obj_data, npcs = num_pcs)
    seurat_obj_data <- harmony::RunHarmony(seurat_obj_data, group.by.vars = c("sample", "batch"))
  }

  print("Running UMAP for Original Data")
  seurat_obj_data <- Seurat::RunUMAP(
    seurat_obj_data,
    dims = 1:num_pcs,
    reduction = reduction_type,
    n.neighbors = 30,
    min.dist = 0.3,
    metric = "cosine",
    n.epochs = NULL,
    learning.rate = 1,
    spread = 1,
    repulsion.strength = 1
  )
  umap_coords_data_pca <- Seurat::Embeddings(seurat_obj_data, reduction = "umap")
  umap_pca_file_name <- paste0(file_prefixes, "umap_coords_data_", reduction_type, ".rds")
  print(paste0("Saving umap coords to file: ", file.path(dir_analysis_dump, umap_pca_file_name)))
  saveRDS(umap_coords_data_pca, file.path(dir_analysis_dump, umap_pca_file_name))
  invisible(NULL)
}


