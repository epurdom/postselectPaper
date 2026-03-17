#' Count clusters that are single-condition
#'
#' Loads \code{pa_de}, cluster assignment, and \code{fake_condition} from disk and
#' counts how many clusters have at most one condition with at least
#' \code{cut_off_single_cond} cells.
#'
#' @param curr_id_check Character. Run ID used in file names.
#' @param analysis_type Character. Key into \code{anls_patterns} (e.g. \code{harmony_}).
#' @param cut_off_single_cond Numeric. Minimum cells per condition to count cluster as "both conditions" (default \code{100}).
#' @param anls_info_dir Character. Directory for analysis RDS files.
#' @param pa_de_save_prfx Character. Prefix for \code{pa_de} RDS filenames.
#' @param clustering_prefix Character. Prefix for cluster assignment RDS filenames.
#' @param anls_patterns Named list. Maps \code{analysis_type} to filename infix.
#'
#' @return Integer. Number of clusters with at most one condition above the cutoff.
#'
#' @export
get_num_single_cond_clusts <- function(curr_id_check, analysis_type, cut_off_single_cond,
    anls_info_dir, pa_de_save_prfx, clustering_prefix, anls_patterns) {
  pa_de_int_fn <- paste0(anls_info_dir, pa_de_save_prfx, curr_id_check, ".rds")
  pa_de_int <- readRDS(pa_de_int_fn)

  cluster_assignment_fn <- paste0(anls_info_dir, clustering_prefix, anls_patterns[[analysis_type]], curr_id_check, ".rds")
  cluster_assignment <- readRDS(cluster_assignment_fn)

  fake_condition_fn <- pa_de_int$fake_condition_fn
  fake_condition <- readRDS(fake_condition_fn)

  cont_table <- table(fake_condition, cluster_assignment)
  is_single_cond <- apply(cont_table, 2, function(x) sum( x > cut_off_single_cond) <2)
  num_single_cond_clusts <- sum(is_single_cond)
  return(num_single_cond_clusts)
}

#' Adjusted Rand Index between current and reference clustering
#'
#' Loads cluster assignments from disk. Reference can be ground truth
#' (\code{ref_analysis_type == "gt"}) or another analysis type.
#'
#' @param curr_id_check Character. Run ID.
#' @param analysis_type Character. Key for current clustering in \code{anls_patterns}.
#' @param ref_analysis_type Character. \code{"gt"} for ground truth, or key for reference clustering in \code{anls_patterns}.
#' @param anls_info_dir Character. Directory for analysis RDS files.
#' @param pa_de_save_prfx Character. Prefix for \code{pa_de} RDS filenames.
#' @param clustering_prefix Character. Prefix for cluster assignment RDS filenames.
#' @param anls_patterns Named list. Maps analysis type to filename infix.
#'
#' @return Numeric. ARI from \code{mclust::adjustedRandIndex}.
#'
#' @importFrom mclust adjustedRandIndex
#' @export
get_ref_ari <- function(curr_id_check, analysis_type, ref_analysis_type,
    anls_info_dir, pa_de_save_prfx, clustering_prefix, anls_patterns) {
  pa_de_int_fn <- paste0(anls_info_dir, pa_de_save_prfx, curr_id_check, ".rds")
  pa_de_int <- readRDS(pa_de_int_fn)

  cluster_assignment_fn <- paste0(anls_info_dir, clustering_prefix, anls_patterns[[analysis_type]], curr_id_check, ".rds")
  cluster_assignment <- readRDS(cluster_assignment_fn)

  if (ref_analysis_type == "gt") {
    if (!file.exists(pa_de_int$reduced_par_list_fn)) {
      stop("Reduced par list file does not exist")
    }
    reduced_par_list_fn <- pa_de_int$reduced_par_list_fn
    reduced_par <- readRDS(reduced_par_list_fn)
    cluster_assignment_ref <- reduced_par$cd$cluster_id
  } else {
    cluster_assignment_ref_fn <- paste0(anls_info_dir, clustering_prefix, anls_patterns[[ref_analysis_type]], curr_id_check, ".rds")
    cluster_assignment_ref <- readRDS(cluster_assignment_ref_fn)
  }

  ari_clusts <- mclust::adjustedRandIndex(cluster_assignment, cluster_assignment_ref)
  print(paste0("ARI for clusters: ", ari_clusts))
  return(ari_clusts)
}

#' Chi-squared statistic for condition x cluster contingency table
#'
#' Loads cluster assignment and \code{fake_condition} from disk and returns the
#' chi-squared test statistic for the condition-by-cluster table.
#'
#' @param curr_id_check Character. Run ID.
#' @param analysis_type Character. Key for \code{anls_patterns}.
#' @param anls_info_dir Character. Directory for analysis RDS files.
#' @param pa_de_save_prfx Character. Prefix for \code{pa_de} RDS filenames.
#' @param clustering_prefix Character. Prefix for cluster assignment RDS filenames.
#' @param anls_patterns Named list. Maps \code{analysis_type} to filename infix.
#'
#' @return Numeric. \code{chisq.test(cont_table)$statistic}.
#'
#' @importFrom stats chisq.test
#' @export
get_chisq_stat <- function(curr_id_check, analysis_type, anls_info_dir, pa_de_save_prfx, clustering_prefix, anls_patterns) {
  pa_de_int_fn <- paste0(anls_info_dir, pa_de_save_prfx, curr_id_check, ".rds")
  pa_de_int <- readRDS(pa_de_int_fn)

  cluster_assignment_fn <- paste0(anls_info_dir, clustering_prefix, anls_patterns[[analysis_type]], curr_id_check, ".rds")
  cluster_assignment <- readRDS(cluster_assignment_fn)

  fake_condition_fn <- pa_de_int$fake_condition_fn
  fake_condition <- readRDS(fake_condition_fn)

  cont_table <- table(fake_condition, cluster_assignment)
  chisq_stat <- chisq.test(cont_table)$statistic
  return(chisq_stat)
}

#' DE and non-false gene sets for one cluster from overlap info
#'
#' Uses \code{cell_count_assign_de_gene}, \code{gt_size_de_gene}, and
#' \code{assign_size_clust} from overlap info to compute overlap proportions and
#' define DE and non-false gene sets by cut-offs.
#'
#' @param clust_name Character. Name of the cluster.
#' @param cut_off_true Numeric. Proportion threshold above which a gene is in the DE set.
#' @param cut_off_false Numeric. Proportion threshold above which a gene is in the non-false set.
#' @param overlap_type Character. Denominator for proportion: \code{"union"}, \code{"clust"}, or \code{"gt"}.
#' @param de_gene_overlap_info List from \code{\link{calc_overlap_mSim}} or \code{\link{calc_overlap_augData}}.
#'
#' @return A list with \code{de_genes_set} and \code{non_false_genes_set} (character vectors of gene names).
#'
#' @export
get_de_genes_sets_w_overlap <- function(clust_name, cut_off_true, cut_off_false, overlap_type, de_gene_overlap_info) {
# get the size of the intersection between the cluster and the ground truth clusters
  # overlap_size <- cont_table[clust_name,]
  nominator <- de_gene_overlap_info$cell_count_assign_de_gene[,clust_name]
  gs_idx <- rownames(de_gene_overlap_info$cell_count_assign_de_gene)
  gt_size_de_gene <- de_gene_overlap_info$gt_size_de_gene
  assign_size_clust <- de_gene_overlap_info$assign_size_clust[clust_name]
  assign_size_clust <- as.numeric(assign_size_clust)
  if (length(assign_size_clust) == 0L || is.na(assign_size_clust)) assign_size_clust <- 0

  # depending on the overlap type, get the denominator
  # gt_clust_sums <- colSums(cont_table)
  # data_clust_size <- sum(cont_table[clust_name,])
  denominator <- switch(overlap_type,
    "union" = assign_size_clust+ gt_size_de_gene,
    "clust" = assign_size_clust,
    "gt" = gt_size_de_gene,
    stop("Invalid overlap type")
  )

  # calculate the percentage overlap
  perc_overlap <- nominator / denominator

  # calculate the sets of DE and non-false genes using the cut-offs
  de_genes_set <- gs_idx[perc_overlap > cut_off_true]
  non_false_genes_set <- gs_idx[perc_overlap > cut_off_false]
  return(list(de_genes_set = de_genes_set, non_false_genes_set = non_false_genes_set))
}

#' Mean PVE (cluster) among top fraction of genes by variance
#'
#' Loads \code{var_fracs} from disk, ranks genes by cluster PVE, and returns the
#' mean cluster PVE over the top \code{pct_top_var_genes} fraction.
#'
#' @param curr_id_check Character. Run ID.
#' @param analysis_type Character. Key for \code{anls_patterns}.
#' @param pct_top_var_genes Numeric. Fraction of genes (by rank) to average (default \code{0.5}).
#' @param anls_info_dir Character. Directory for analysis RDS files.
#' @param var_fracs_save_prfx Character. Prefix for variance-fractions RDS filenames.
#' @param anls_patterns Named list. Maps \code{analysis_type} to filename infix.
#'
#' @return Numeric. Mean cluster PVE for the top genes.
#'
#' @export
get_avg_top_pve_de_genes <- function(curr_id_check, analysis_type, pct_top_var_genes,
    anls_info_dir, var_fracs_save_prfx, anls_patterns) {
  var_fracs_fn <- paste0(anls_info_dir, var_fracs_save_prfx, anls_patterns[[analysis_type]], curr_id_check, ".rds")
  var_fracs <- readRDS(var_fracs_fn)
  var_fracs_sorted <- var_fracs[order(var_fracs$cluster, decreasing = TRUE),]
  num_top_var_genes <- round(pct_top_var_genes*nrow(var_fracs_sorted))
  avg_pve_de_genes <- mean(var_fracs_sorted$cluster[1:num_top_var_genes])
  return(avg_pve_de_genes)
}

#' Redundancy index and number of DE genes found
#'
#' Loads CCA analysis from disk and returns the redundancy index plus the count
#' of unique genes significant at \code{sig_threshold}.
#'
#' @param curr_id_check Character. Run ID.
#' @param analysis_type Character. Key for \code{anls_patterns}.
#' @param res_DE List. Per-cluster DE result tables (must have \code{p_adj.glb} and \code{gene}).
#' @param num_genes Unused; kept for interface compatibility.
#' @param sig_threshold Numeric. FDR threshold for counting DE genes (default \code{0.10}).
#' @param anls_info_dir Character. Directory for analysis RDS files.
#' @param cca_anls_save_prfx Character. Prefix for CCA analysis RDS filenames.
#' @param anls_patterns Named list. Maps \code{analysis_type} to filename infix.
#'
#' @return List with \code{redundancy_idx} (from \code{cca_anls$Y.redun}) and \code{num_de_found}.
#'
#' @export
get_redundancy_idx <- function(curr_id_check, analysis_type, res_DE, num_genes, sig_threshold,
    anls_info_dir, cca_anls_save_prfx, anls_patterns) {
  cca_anls_fn <- paste0(anls_info_dir, cca_anls_save_prfx, anls_patterns[[analysis_type]], curr_id_check, ".rds")
  cca_anls <- readRDS(cca_anls_fn)

  # get number of DE genes found for comparison
  for (clust_name_iter in names(res_DE)){
    res_DE[[clust_name_iter]]$clust_name <- clust_name_iter
  }
  united_res <- data.frame(do.call(rbind, res_DE))
  num_de_found <- length(unique(united_res[united_res$p_adj.glb < sig_threshold,"gene"]))
  
  return(list(redundancy_idx = cca_anls$Y.redun, num_de_found = num_de_found))
}


#' Cluster-by-cluster FDP/TPP control using overlap-defined DE and non-false sets
#'
#' For each cluster, gets DE and non-false gene sets via
#' \code{\link{get_de_genes_sets_w_overlap}}, merges DE results, BH-adjusts
#' p-values, and computes local and global FDP/TPP via \code{\link{calc_fdp}}.
#'
#' @param res_DE List. Per-cluster DE result tables; each element has \code{gene}, \code{p_val}, etc.
#' @param cut_off_true,cut_off_false Numeric. Overlap proportion cut-offs for DE and non-false sets.
#' @param sig_threshold Numeric. Significance threshold for DE calls.
#' @param overlap_type Character. \code{"union"}, \code{"clust"}, or \code{"gt"} (denominator for overlap proportion).
#' @param de_gene_overlap_info List from \code{\link{calc_overlap_mSim}} or \code{\link{calc_overlap_augData}}.
#' @param dd_genes_sets Optional list of precomputed DE/non-false sets per cluster; if \code{NULL}, computed via \code{\link{get_de_genes_sets_w_overlap}}.
#' @param once_per_gene Logical. If \code{TRUE}, count each gene at most once for global TPP (default \code{FALSE}).
#'
#' @return List with \code{ctrl_clusts_loc}, \code{ctrl_glob}, \code{ctrl_clusts_glob} (data frames with \code{tp}, \code{fp}, \code{fdp}, \code{tpp}, etc.).
#'
#' @importFrom dplyr group_by filter ungroup
#' @importFrom rlang .data
#' @importFrom stats p.adjust
#' @export
run_clust_by_clust_w_overlap <- function(res_DE, cut_off_true, cut_off_false, sig_threshold, 
  overlap_type, de_gene_overlap_info, dd_genes_sets = NULL, once_per_gene = FALSE){
  # initialize the data frame and list to be returned
  sim_res_DE <- data.frame(gene = character(), clust_name = character(), is_de = numeric(),
                            not_de = numeric(), clust_adj_pvals = numeric(), p_val = numeric())
  ctrl_clusts_df_sim <- data.frame(clust_name = character(),
                                  tp = numeric(), fp = numeric(), fdp = numeric())
  clust_name <- names(res_DE)[1]
  total_num_de_genes <- 0
  for (clust_name in names(res_DE)){
    # find the gene sets for true and false positives if not given in the function call
    if (is.null(dd_genes_sets)){
      dd_genes_set <- get_de_genes_sets_w_overlap(clust_name, cut_off_true, cut_off_false,
                                                overlap_type, de_gene_overlap_info)
    } else {
      dd_genes_set <- dd_genes_sets[[clust_name]]
    }
    de_genes_set <- dd_genes_set$de_genes_set
    non_false_genes_set <- dd_genes_set$non_false_genes_set
    curr_num_de_genes <- length(de_genes_set)
    total_num_de_genes <- total_num_de_genes + curr_num_de_genes

    # include the DE, non_False information run
    clust_res_DE <- res_DE[[clust_name]]
    clust_res_DE$clust_adj_pvals <- stats::p.adjust(clust_res_DE$p_val, method = "BH")
    clust_res_DE <- clust_res_DE[!is.na(clust_res_DE$clust_adj_pvals),]
    clust_res_DE$is_de <- as.numeric(clust_res_DE$gene %in% de_genes_set)
    clust_res_DE$not_de <- as.numeric(!(clust_res_DE$gene %in% non_false_genes_set))
    clust_res_DE$clust_name <- clust_name
    clust_res_DE$is_sig_clust <- as.numeric(clust_res_DE$clust_adj_pvals < sig_threshold)# yes, anoying double negative
    sim_res_DE <- rbind(sim_res_DE, clust_res_DE)

    # store the FDP for the cluster on the local family
    ctrl_clusts_loc <- data.frame(calc_fdp(clust_res_DE[clust_res_DE$is_sig_clust == 1,]), 
                                  clust_name = clust_name, num_de_genes = curr_num_de_genes, once_per_gene = once_per_gene)
    ctrl_clusts_loc$tpp <- ctrl_clusts_loc$tp / ctrl_clusts_loc$num_de_genes
    ctrl_clusts_df_sim <- rbind(ctrl_clusts_df_sim, ctrl_clusts_loc)
  }
  sim_res_DE$glob_adj_pvals <- stats::p.adjust(sim_res_DE$p_val, method = "BH")
  sim_res_DE$is_sig_glob <- as.numeric(sim_res_DE$glob_adj_pvals < sig_threshold)
  if (once_per_gene) {
    sim_res_DE <- dplyr::ungroup(dplyr::filter(dplyr::group_by(sim_res_DE, .data$gene), .data$glob_adj_pvals == min(.data$glob_adj_pvals, na.rm = TRUE)))
    num_de_genes_tot <- length(de_gene_overlap_info$gt_size_de_gene)
  } else {
    num_de_genes_tot <- total_num_de_genes
  }
  ctrl_glob <- data.frame(calc_fdp(sim_res_DE[sim_res_DE$is_sig_glob == 1, ]))
  ctrl_glob$once_per_gene <- once_per_gene
  ctrl_glob$tpp <- ctrl_glob$tp / num_de_genes_tot
  ctrl_clusts_glob <- data.frame(calc_fdp(sim_res_DE[sim_res_DE$is_sig_clust == 1, ], once_per_gene = once_per_gene))
  ctrl_clusts_glob$tpp <- ctrl_clusts_glob$tp / num_de_genes_tot
  ctrl_clusts_glob$once_per_gene <- once_per_gene
  return(list(ctrl_clusts_loc = ctrl_clusts_df_sim, ctrl_glob = ctrl_glob, ctrl_clusts_glob = ctrl_clusts_glob))
}

#' False discovery proportion (FDP) and counts for significant DE calls
#'
#' Assumes \code{sub_res_DE} has columns \code{is_de}, \code{not_de}. A DE gene can be counted in multiple clusters.
#'
#' @param sub_res_DE Data frame. Rows = significant genes; columns \code{is_de}, \code{not_de} (numeric 0/1).
#' @param once_per_gene Unused; kept for interface compatibility.
#' @param allow_duplicates Unused; kept for interface compatibility.
#'
#' @return Data frame with one row: \code{tp}, \code{fp}, \code{fdp} (\code{fp/(tp+fp)}).
#'
#' @export
calc_fdp <- function(sub_res_DE, once_per_gene = FALSE, allow_duplicates = TRUE) {
  # with improper clustering, a DE gene can be a multiple false positives
  fp <- sum(sub_res_DE$not_de)
  tp <- sum(sub_res_DE$is_de)

  if (tp+fp > 0){
    fdp <- fp / (tp+fp)
  } else {
    fdp <- 0
  }
  return(data.frame(tp = tp, fp = fp, fdp = fdp))
}

# sim_res_DE, "is_sig_glob", , allow_duplicates = allow_duplicates
# sub_res <- sim_res_DE
# sig_col <- "is_sig_glob"
# allow_duplicates <- TRUE
# de_genes_all <- de_genes_all
#' FDP and TPP given empirical significant/non-significant split
#'
#' Splits \code{sub_res} by \code{sig_col}, then computes tp, fp, fn, fdp, and tpp.
#' If \code{allow_duplicates} is \code{FALSE}, counts unique genes and uses
#' \code{de_genes_all} for the false-negative denominator.
#'
#' @param sub_res Data frame with \code{is_de}, \code{not_de}, and the column named \code{sig_col}.
#' @param sig_col Character. Column name for significance (1 = significant).
#' @param de_genes_all Character. Full set of DE genes (for fn when \code{allow_duplicates = FALSE}).
#' @param allow_duplicates Logical. If \code{TRUE}, sum counts; if \code{FALSE}, count unique genes (default \code{TRUE}).
#'
#' @return Data frame with one row: \code{tp}, \code{fp}, \code{fdp}, \code{tpp}.
#'
#' @export
calc_fdp_given_emp <- function(sub_res, sig_col, de_genes_all = NULL, allow_duplicates = TRUE) {
  # with improper clustering, a DE gene can be a multiple false positives
  sub_res_sig <- sub_res[sub_res[[sig_col]]==1,]
  sub_res_not_sig <- sub_res[sub_res[[sig_col]]==0,]
  fp <- 0
  tp <- 0
  fn <- 0
  if (allow_duplicates){
    fp <- sum(sub_res_sig$not_de)
    tp <- sum(sub_res_sig$is_de)
    fn <- sum(sub_res_not_sig$is_de)
  } else {
    fp_set <- sub_res_sig[sub_res_sig$not_de==1,]
    fp <- length(unique(fp_set$gene))
    tp_set <- sub_res_sig[sub_res_sig$is_de==1,]
    tp <- length(unique(tp_set$gene))
    fn <- length(setdiff(de_genes_all, tp_set$gene))
  }

  if (tp+fp > 0){
    fdp <- fp / (tp+fp)
  } else {
    fdp <- 0
  }

  if (tp+fn > 0){
    tpp <- tp / (tp+fn)
  } else {
    tpp <- 1
  }
  return(data.frame(tp = tp, fp = fp, fdp = fdp, tpp = tpp))
}

#sub_res <- res_DE
# dd_genes_sets <- NULL
# once_per_gene <- FALSE
#' Cluster-by-cluster FDP/TPP using overlap sets and empirical FDP/TPP
#'
#' Same as \code{\link{run_clust_by_clust_w_overlap}} but uses
#' \code{\link{calc_fdp_given_emp}} for FDP/TPP computation (supports
#' \code{allow_duplicates = FALSE} for unique-gene counting).
#'
#' @param res_DE List. Per-cluster DE result tables.
#' @param cut_off_true,cut_off_false Numeric. Overlap proportion cut-offs.
#' @param sig_threshold Numeric. Significance threshold.
#' @param overlap_type Character. \code{"union"}, \code{"clust"}, or \code{"gt"}.
#' @param de_gene_overlap_info List from \code{\link{calc_overlap_mSim}} or \code{\link{calc_overlap_augData}}.
#' @param dd_genes_sets Optional list. Precomputed DE/non-false sets per cluster.
#' @param once_per_gene Logical. Count each gene at most once for global (default \code{FALSE}).
#' @param allow_duplicates Logical. Passed to \code{\link{calc_fdp_given_emp}} (default \code{TRUE}).
#'
#' @return List with \code{ctrl_clusts_loc}, \code{ctrl_glob}, \code{ctrl_clusts_glob}.
#'
#' @importFrom dplyr group_by filter ungroup
#' @importFrom rlang .data
#' @importFrom stats p.adjust
#' @export
run_clust_by_clust_w_overlap_given_emp <- function(res_DE, cut_off_true, cut_off_false, sig_threshold,
  overlap_type, de_gene_overlap_info, dd_genes_sets = NULL, once_per_gene = FALSE, allow_duplicates = TRUE) {
  # initialize the data frame and list to be returned
  sim_res_DE <- data.frame(gene = character(), clust_name = character(), is_de = numeric(),
                            not_de = numeric(), clust_adj_pvals = numeric(), p_val = numeric())
  ctrl_clusts_df_sim <- data.frame(clust_name = character(),
                                  tp = numeric(), fp = numeric(), fdp = numeric())
  clust_name <- names(res_DE)[1]
  total_num_de_genes <- 0
  de_genes_all <- rownames(de_gene_overlap_info$cell_count_assign_de_gene)
  for (clust_name in names(res_DE)){
    # find the gene sets for true and false positives if not given in the function call
    if (is.null(dd_genes_sets)){
      dd_genes_set <- get_de_genes_sets_w_overlap(clust_name, cut_off_true, cut_off_false,
                                                overlap_type, de_gene_overlap_info)
    } else {
      dd_genes_set <- dd_genes_sets[[clust_name]]
    }
    de_genes_set <- dd_genes_set$de_genes_set
    non_false_genes_set <- dd_genes_set$non_false_genes_set
    curr_num_de_genes <- length(de_genes_set)
    total_num_de_genes <- total_num_de_genes + curr_num_de_genes

    # include the DE, non_False information run
    clust_res_DE <- res_DE[[clust_name]]
    clust_res_DE$clust_adj_pvals <- stats::p.adjust(clust_res_DE$p_val, method = "BH")
    clust_res_DE <- clust_res_DE[!is.na(clust_res_DE$clust_adj_pvals),]
    clust_res_DE$is_de <- as.numeric(clust_res_DE$gene %in% de_genes_set)
    clust_res_DE$not_de <- as.numeric(!(clust_res_DE$gene %in% non_false_genes_set))
    clust_res_DE$clust_name <- clust_name
    clust_res_DE$is_sig_clust <- as.numeric(clust_res_DE$clust_adj_pvals < sig_threshold)# yes, anoying double negative
    sim_res_DE <- rbind(sim_res_DE, clust_res_DE)

    # store the FDP for the cluster on the local family
    ctrl_clusts_loc <- data.frame(calc_fdp_given_emp(clust_res_DE, "is_sig_clust", de_genes_all = de_genes_all, allow_duplicates = allow_duplicates), 
                                  clust_name = clust_name, num_de_genes = curr_num_de_genes, once_per_gene = once_per_gene)
    ctrl_clusts_df_sim <- rbind(ctrl_clusts_df_sim, ctrl_clusts_loc)
  }
  sim_res_DE$glob_adj_pvals <- stats::p.adjust(sim_res_DE$p_val, method = "BH")
  sim_res_DE$is_sig_glob <- as.numeric(sim_res_DE$glob_adj_pvals < sig_threshold)
  if (once_per_gene) {
    sim_res_DE <- dplyr::ungroup(dplyr::filter(dplyr::group_by(sim_res_DE, .data$gene), .data$glob_adj_pvals == min(.data$glob_adj_pvals, na.rm = TRUE)))
    num_de_genes_tot <- length(de_gene_overlap_info$gt_size_de_gene)
  } else {
    num_de_genes_tot <- total_num_de_genes
  }
  ctrl_glob <- data.frame(calc_fdp_given_emp(sim_res_DE, "is_sig_glob", de_genes_all = de_genes_all, allow_duplicates = allow_duplicates))
  ctrl_glob$once_per_gene <- once_per_gene
  ctrl_clusts_glob <- data.frame(calc_fdp_given_emp(sim_res_DE, "is_sig_clust", de_genes_all = de_genes_all, allow_duplicates = allow_duplicates))
  ctrl_clusts_glob$once_per_gene <- once_per_gene
  return(list(ctrl_clusts_loc = ctrl_clusts_df_sim, ctrl_glob = ctrl_glob, ctrl_clusts_glob = ctrl_clusts_glob))
}

