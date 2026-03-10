#' Extract and merge per-simulation DE results, overlap, and PVE metrics
#'
#' Loads DE p-values by cluster, overlap matrix, and PVE metrics for a single
#' simulation ID, then merges them into one data frame. Applies
#' \code{adjust_all_pvals} for global testing and \code{two_stage_adjustment}
#' for min_holm, fisher, and cauchy screening (from \code{\link{combine_de_pvals_by_cluster}} / hypothesis adjustment).
#'
#' @param id_check_single Character. Single simulation run identifier (e.g. \code{"1"}).
#' @param phen_type_removal Character. Phenotype/removal label used in file names.
#' @param analysis_results_dir Character. Directory for analysis outputs (overlap, pa_de, PVE RDS files).
#' @param de_outputs_dir Character. Directory containing DE p-values by cluster RDS files.
#' @param overlap_matrix_prefix Character. Filename prefix for overlap matrix: \code{<prefix><phen_type_removal><id>.rds}.
#' @param de_pvals_by_cluster_prefix Character. Filename prefix for DE p-values by cluster RDS.
#' @param pa_de_save_prfx Character. Filename prefix for \code{pa_de} RDS (contains \code{num_de_genes}).
#' @param PVE_metrics_save_prfx Character. Filename prefix for PVE metrics RDS (\code{tPVE}, \code{sPVE} by gene).
#'
#' @return A data.frame with gene-cluster rows, DE results, overlap counts, PVE, and 2-stage
#'   columns for \code{min_holm}, \code{fisher}, and \code{cauchy}. Returns \code{NULL} if the
#'   overlap matrix file is missing.
#'
#' @seealso \code{\link{combine_de_pvals_by_cluster}}, \code{\link{confusion_from_overlap_matrix}}
#' @importFrom reshape2 melt
#' @export
per_sim_data_extraction <- function(id_check_single, phen_type_removal, analysis_results_dir, de_outputs_dir, overlap_matrix_prefix, de_pvals_by_cluster_prefix, pa_de_save_prfx, PVE_metrics_save_prfx) {
  pa_de_int_fn <- paste0(analysis_results_dir, pa_de_save_prfx, id_check_single, ".rds")
  pa_de_int <- readRDS(pa_de_int_fn)
  num_de_genes <- pa_de_int$num_de_genes

  overlap_matrix_fn <- paste0(
    analysis_results_dir,
    overlap_matrix_prefix,
    phen_type_removal,
    id_check_single,
    ".rds"
  )
  if (!file.exists(overlap_matrix_fn)) {
    message("File does not exist: ", overlap_matrix_fn)
    return(NULL)
  }
  overlap_matrix <- readRDS(overlap_matrix_fn)
  overlap_count_matrix <- overlap_matrix$cell_count_assign_de_gene
  overlap_count_long <- reshape2::melt(overlap_count_matrix, varnames = c("gene", "cluster_id"), value.name = "overlap_count")

  sim_res_DE_fn <- paste0(de_outputs_dir, de_pvals_by_cluster_prefix, phen_type_removal, id_check_single, ".rds")
  sim_res_DE <- readRDS(sim_res_DE_fn)
  combined_sim_res_DE_tested <- adjust_all_pvals(sim_res_DE)
  combined_sim_res_DE_tested$tested <- TRUE

  combined_sim_res_DE_all <- merge(
    combined_sim_res_DE_tested, overlap_count_long,
    by = c("gene", "cluster_id"), all.y = TRUE, all.x = TRUE
  )
  combined_sim_res_DE_all$gt_nbhd_count <- overlap_matrix$gt_size_de_gene[combined_sim_res_DE_all$gene]
  combined_sim_res_DE_all$found_cluster_count <- overlap_matrix$assign_size_clust[combined_sim_res_DE_all$cluster_id]

  combined_sim_res_DE_all$tested[is.na(combined_sim_res_DE_all$tested)] <- FALSE
  combined_sim_res_DE_all$gt_nbhd_count[is.na(combined_sim_res_DE_all$gt_nbhd_count)] <- 0
  combined_sim_res_DE_all$overlap_count[is.na(combined_sim_res_DE_all$overlap_count)] <- 0

  combined_sim_res_DE_all$id_check <- id_check_single
  combined_sim_res_DE_all$phen_type_removal <- phen_type_removal
  combined_sim_res_DE_all$num_de_genes <- num_de_genes
  
  # add tPVE and sPVE values for each gene
  PVE_metrics_fn <- paste0(analysis_results_dir, PVE_metrics_save_prfx, phen_type_removal, id_check_single, ".rds")
  PVE_metrics <- readRDS(PVE_metrics_fn)
  combined_sim_res_DE_all$tPVE <- PVE_metrics$tPVE[as.character(combined_sim_res_DE_all$gene)]
  combined_sim_res_DE_all$sPVE <- PVE_metrics$sPVE[as.character(combined_sim_res_DE_all$gene)]


  poss_screen_methods <- c("min_holm", "fisher", "cauchy")
  for (screen_method in poss_screen_methods) {
    combined_sim_res_DE_2stage <- two_stage_adjustment(sim_res_DE, screen_method = screen_method)
    tested_col <- paste0("tested_2stage_", screen_method)
    combined_sim_res_DE_2stage[[tested_col]] <- TRUE
    vals_to_merge_in <- combined_sim_res_DE_2stage[, c("gene", "cluster_id", "adj_p_val_2stage", "screen_pval", "screen_adj_pval", tested_col)]
    combined_sim_res_DE_all <- merge(
      combined_sim_res_DE_all, vals_to_merge_in,
      by = c("gene", "cluster_id"), all.x = TRUE
    )
    screen_pval_col_new <- paste0("screen_pval_", screen_method)
    screen_adj_pval_col_new <- paste0("screen_adj_pval_", screen_method)
    adj_pval_2stage_col_new <- paste0("adj_p_val_2stage_", screen_method)

    combined_sim_res_DE_all[[tested_col]][is.na(combined_sim_res_DE_all[[tested_col]])] <- FALSE

    names(combined_sim_res_DE_all)[names(combined_sim_res_DE_all) == "screen_pval"] <- screen_pval_col_new
    names(combined_sim_res_DE_all)[names(combined_sim_res_DE_all) == "screen_adj_pval"] <- screen_adj_pval_col_new
    names(combined_sim_res_DE_all)[names(combined_sim_res_DE_all) == "adj_p_val_2stage"] <- adj_pval_2stage_col_new
  }

  return(combined_sim_res_DE_all)
}




#' Compute all-hypothesis (cluster-level) confusion metrics and FDP / POWER
#'
#' Computes true/false positives and negatives from a master power table by
#' defining true DE via overlap percentage and significance via a p-value column.
#' Aggregates TP, FN, FP, POWER, FDP (and optional PVE averages) per \code{id_check}.
#'
#' @param master_pwr_table_used Data frame with at least \code{id_check}, \code{overlap_count},
#'   \code{gt_nbhd_count}, \code{found_cluster_count}, and a p-value column (see \code{p_val_col_name}).
#'   Optionally \code{sPVE} and \code{tPVE} for average PVE of found DEs.
#' @param sig_threshold_fixed Numeric. Significance threshold (alpha); findings with p-value below this are "found DE".
#' @param p_val_col_name Character. Name of the column containing the DE-test p-value (or adjusted p-value). Default \code{"adj_p_val"}.
#' @param overlap_def Character. Overlap denominator: \code{"union"} (default), \code{"clust"}, or \code{"gt"}.
#' @param cut_off_true Numeric. Gene-cluster pairs with \code{perc_overlap > cut_off_true} are considered truly DE. Default 0.
#' @param cut_off_false Numeric. Gene-cluster pairs with \code{perc_overlap <= cut_off_false} are considered not DE. Default 0.
#'
#' @return A data frame with one row per \code{id_check} and columns: \code{TP}, \code{FN}, \code{FP},
#'   \code{FN_all}, \code{POWER}, \code{POWER_all}, \code{FDP}, and optionally \code{sPVE_avg_found_DE},
#'   \code{tPVE_avg_found_DE}. \code{FN} excludes not-tested; \code{FN_all} includes them. \code{POWER}
#'   and \code{FDP} are set to \code{NA} or \code{0} when denominators are zero.
#'
#' @seealso \code{\link{confusion_from_overlap_matrix}}
#' @importFrom dplyr group_by summarize
#' @export
get_conf_metrics_tested_all_hypotheses <- function(master_pwr_table_used, sig_threshold_fixed, p_val_col_name = "adj_p_val", overlap_def = "union", cut_off_true = 0.00, cut_off_false = 0.00) {
  if (overlap_def == "union") {
    denom_overlap_count <- master_pwr_table_used$gt_nbhd_count + master_pwr_table_used$found_cluster_count - master_pwr_table_used$overlap_count
    master_pwr_table_used$perc_overlap <- master_pwr_table_used$overlap_count / denom_overlap_count
  } else if (overlap_def == "clust") {
    denom_overlap_count <- master_pwr_table_used$found_cluster_count
    master_pwr_table_used$perc_overlap <- master_pwr_table_used$overlap_count / denom_overlap_count
  } else if (overlap_def == "gt") {
    denom_overlap_count <- master_pwr_table_used$gt_nbhd_count
    master_pwr_table_used$perc_overlap <- master_pwr_table_used$overlap_count / denom_overlap_count
  } else {
    stop("Invalid overlap definition")
  }
  master_pwr_table_used$perc_overlap[is.na(master_pwr_table_used$perc_overlap)] <- 0
  master_pwr_table_used$is_DE <- master_pwr_table_used$perc_overlap > cut_off_true
  master_pwr_table_used$is_not_DE <- master_pwr_table_used$perc_overlap <= cut_off_false

  master_pwr_table_used$not_tested <- is.na(master_pwr_table_used[[p_val_col_name]])
  master_pwr_table_used[[p_val_col_name]][is.na(master_pwr_table_used[[p_val_col_name]])] <- 1
  master_pwr_table_used$found_DE <- master_pwr_table_used[[p_val_col_name]] < sig_threshold_fixed

  master_pwr_table_used$TP <- master_pwr_table_used$is_DE & master_pwr_table_used$found_DE
  master_pwr_table_used$FN <- master_pwr_table_used$is_DE & !master_pwr_table_used$found_DE & !master_pwr_table_used$not_tested
  master_pwr_table_used$FN_all <- master_pwr_table_used$is_DE & !master_pwr_table_used$found_DE

  master_pwr_table_used$FP <- master_pwr_table_used$is_not_DE & master_pwr_table_used$found_DE

  # if master_pwr_table_used does not have sPVE and tPVE columns, add them
  if (!("sPVE" %in% colnames(master_pwr_table_used))) {
    master_pwr_table_used$sPVE <- NA_real_
  }
  if (!("tPVE" %in% colnames(master_pwr_table_used))) {
    master_pwr_table_used$tPVE <- NA_real_
  }

  conf_metrics_tested <-   master_pwr_table_used %>%
    dplyr::group_by(id_check) %>%
    dplyr::summarize(
      TP = sum(TP, na.rm = TRUE),
      FN = sum(FN, na.rm = TRUE),
      FP = sum(FP, na.rm = TRUE),
      FN_all = sum(FN_all, na.rm = TRUE),
      sPVE_avg_found_DE = if (sum(found_DE, na.rm = TRUE) > 0) mean(sPVE[found_DE], na.rm = TRUE) else NA_real_,
      tPVE_avg_found_DE = if (sum(found_DE, na.rm = TRUE) > 0) mean(tPVE[found_DE], na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )
  denom_power <- conf_metrics_tested$TP + conf_metrics_tested$FN
  conf_metrics_tested$POWER <- ifelse(denom_power > 0, conf_metrics_tested$TP / denom_power, NA_real_)
  denom_power_all <- conf_metrics_tested$TP + conf_metrics_tested$FN_all
  conf_metrics_tested$POWER_all <- ifelse(denom_power_all > 0, conf_metrics_tested$TP / denom_power_all, NA_real_)
  denom_fdp <- conf_metrics_tested$FP + conf_metrics_tested$TP
  conf_metrics_tested$FDP <- ifelse(denom_fdp > 0, conf_metrics_tested$FP / denom_fdp, 0)
  return(conf_metrics_tested)
}




#' Compute gene-level overlap-based confusion metrics (ov_TP, ov_FP, ov_FN) and ov_FDP / ov_POWER
#'
#' First aggregates gene-cluster rows to gene level (per \code{id_check} and \code{gene}) using
#' screening p-values and overlap-based true DE, then defines gene-level TP/FP/FN and
#' sample-level \code{ov_TP}, \code{ov_FN}, \code{ov_FP}, \code{ov_POWER}, \code{ov_FDP} per \code{id_check}.
#'
#' @param master_pwr_table_used Data frame with at least \code{id_check}, \code{gene},
#'   \code{overlap_count}, \code{gt_nbhd_count}, \code{found_cluster_count}, \code{num_de_genes},
#'   and columns for screening and DE p-values (see \code{screen_p_val_col_name}, \code{p_val_col_name}).
#' @param sig_threshold_fixed Numeric. Significance threshold (alpha) for both screen and DE decisions.
#' @param screen_p_val_col_name Character. Column name for screening p-value (e.g. combined/gene-level). Default \code{"adj_p_val"}.
#' @param p_val_col_name Character. Column name for DE-test p-value (cluster-level). Default \code{"adj_p_val"}.
#' @param overlap_def Character. Overlap denominator: \code{"union"} (default), \code{"clust"}, or \code{"gt"}.
#' @param cut_off_true Numeric. Gene-cluster pairs with \code{perc_overlap > cut_off_true} are truly DE. Default 0.
#' @param cut_off_false Numeric. Gene-cluster pairs with \code{perc_overlap <= cut_off_false} are not DE. Default 0.
#'
#' @return A data frame with one row per \code{id_check} and columns \code{ov_TP}, \code{ov_FN}, \code{ov_FP},
#'   \code{ov_POWER}, \code{ov_FDP}, \code{num_de_genes}. \code{ov_POWER} and \code{ov_FDP} are set to
#'   \code{NA} or \code{0} when denominators are zero.
#'
#' @seealso \code{\link{get_conf_metrics_tested_all_hypotheses}}, \code{\link{confusion_from_overlap_matrix}}
#' @importFrom dplyr group_by summarize first
#' @export
get_conf_metrics_tested_by_gene <- function(master_pwr_table_used, sig_threshold_fixed, screen_p_val_col_name = "adj_p_val", p_val_col_name = "adj_p_val", overlap_def = "union", cut_off_true = 0.00, cut_off_false = 0.00) {
  if (overlap_def == "union") {
    denom_overlap_count <- master_pwr_table_used$gt_nbhd_count + master_pwr_table_used$found_cluster_count - master_pwr_table_used$overlap_count
    master_pwr_table_used$perc_overlap <- master_pwr_table_used$overlap_count / denom_overlap_count
  } else if (overlap_def == "clust") {
    denom_overlap_count <- master_pwr_table_used$found_cluster_count
    master_pwr_table_used$perc_overlap <- master_pwr_table_used$overlap_count / denom_overlap_count
  } else if (overlap_def == "gt") {
    denom_overlap_count <- master_pwr_table_used$gt_nbhd_count
    master_pwr_table_used$perc_overlap <- master_pwr_table_used$overlap_count / denom_overlap_count
  } else {
    stop("Invalid overlap definition")
  }

  master_pwr_table_used$screen_p_val <- master_pwr_table_used[[screen_p_val_col_name]]

  # for genes where no overlap was provided, the percent overlap is 0
  master_pwr_table_used$perc_overlap[is.na(master_pwr_table_used$perc_overlap)] <- 0
  
  master_pwr_table_used$not_tested <- is.na(master_pwr_table_used[[p_val_col_name]])
  master_pwr_table_used[[p_val_col_name]][is.na(master_pwr_table_used[[p_val_col_name]])] <- 1

  # genes are DE (tested) when overlap > cut_off_true and were tested
  master_pwr_table_used$is_DE <- master_pwr_table_used$perc_overlap > cut_off_true & !master_pwr_table_used$not_tested
  # genes are DE (all) when overlap > cut_off_true, regardless of tested
  master_pwr_table_used$is_DE_all <- master_pwr_table_used$perc_overlap > cut_off_true

  # genes are not DE when the overlap is less than the false cutoff
  master_pwr_table_used$is_not_DE <- master_pwr_table_used$perc_overlap <= cut_off_false

  # genes are found to be DE when the p-value is less than the significance threshold
  master_pwr_table_used$found_DE <- master_pwr_table_used[[p_val_col_name]] <= sig_threshold_fixed




  message("Calculating gene-level confusion metrics - over all clusters for a given simulation and gene")
  gene_level_agg <- master_pwr_table_used %>%
    dplyr::group_by(id_check, gene) %>%
    dplyr::summarize(
        TP_gene_clust_pair = sum(is_DE & found_DE, na.rm = TRUE), # true positives are genes that are DE and found to be DE
        FP_gene_clust_pair = sum((!is_DE) & found_DE, na.rm = TRUE), # false positives are genes that are not DE and found to be DE
        FN_gene_clust_pair = sum(is_DE & (!found_DE), na.rm = TRUE), # false negatives are genes that are DE and not found to be DE
        FN_gene_clust_pair_all = sum(is_DE_all & (!found_DE), na.rm = TRUE), # false negatives are genes that are DE and not found to be DE
        min_screen_p_val = min(screen_p_val, na.rm = TRUE), # aggregate the screen p-values by taking min per gene
        gene_level_null = sum(!is_not_DE) == 0, # if all clusters for a given gene are not DE, then the gene is a null gene (this is null for screen pvalue)
        num_de_genes = dplyr::first(num_de_genes),
        .groups = "drop"
    )

  # defining screen level decisions and TP/FN/FP
  gene_level_agg$screen_rejected_gene <- gene_level_agg$min_screen_p_val <= sig_threshold_fixed
  gene_level_agg$TP_gene_screen <- (!gene_level_agg$gene_level_null) & gene_level_agg$screen_rejected_gene
  gene_level_agg$FP_gene_screen <- gene_level_agg$gene_level_null & gene_level_agg$screen_rejected_gene
  gene_level_agg$FN_gene_screen <- !gene_level_agg$gene_level_null & !gene_level_agg$screen_rejected_gene
  
  # defining augmented V_g, S_g, T_g, R_g
  gene_level_agg$V_g <- gene_level_agg$FP_gene_clust_pair + as.numeric(gene_level_agg$FP_gene_screen)
  gene_level_agg$S_g <- gene_level_agg$TP_gene_clust_pair + as.numeric(gene_level_agg$TP_gene_screen)
  gene_level_agg$T_g <- gene_level_agg$FN_gene_clust_pair + as.numeric(gene_level_agg$FN_gene_screen)
  gene_level_agg$R_g <- gene_level_agg$V_g + gene_level_agg$S_g
  
  gene_level_agg$TP_gene_level_ov <- (gene_level_agg$V_g == 0) & (gene_level_agg$R_g > 0)
  gene_level_agg$FP_gene_level_ov <- (gene_level_agg$V_g > 0) 
  gene_level_agg$FN_gene_level_ov <- (gene_level_agg$T_g > 0) & (gene_level_agg$R_g == 0)

  message("Calculating sample level confusion metrics")
  conf_metrics_tested <- gene_level_agg %>%
    dplyr::group_by(id_check) %>%
    dplyr::summarize(
      ov_TP = sum(TP_gene_level_ov, na.rm = TRUE),
      ov_FN = sum(FN_gene_level_ov, na.rm = TRUE),
      ov_FP = sum(FP_gene_level_ov, na.rm = TRUE),
      num_de_genes = dplyr::first(num_de_genes),
      .groups = "drop"
    )
  denom_ov_power <- conf_metrics_tested$ov_TP + conf_metrics_tested$ov_FN
  conf_metrics_tested$ov_POWER <- ifelse(denom_ov_power > 0, conf_metrics_tested$ov_TP / denom_ov_power, NA_real_)
  denom_ov_fdp <- conf_metrics_tested$ov_FP + conf_metrics_tested$ov_TP
  conf_metrics_tested$ov_FDP <- ifelse(denom_ov_fdp > 0, conf_metrics_tested$ov_FP / denom_ov_fdp, 0)
  return(conf_metrics_tested)
}