#' Combine differential expression p-values by cluster
#'
#' This is the main wrapper function that combines DE p-values from multiple
#' clusters using one of several methods: adjusting all p-values together,
#' taking the minimum p-value per cluster, or using Fisher's combined probability
#' test with minimum selection.
#'
#' @param de_pvals_by_cluster A named list where each element corresponds to a
#'   cluster and contains a data.frame with columns 'gene', 'cluster_id', and 'p_val'.
#' @param method Character string specifying the combination method. One of:
#'   \itemize{
#'     \item "adjust_all": Adjust all p-values across all clusters together.
#'     \item "min_per_cluster": Take minimum p-value per gene across clusters.
#'     \item "fisher_combined_and_min": Use Fisher's method to combine p-values
#'           per gene and report the cluster with minimum p-value.
#'     \item "2-stage-fisher": Two-stage procedure with Fisher-combined screening.
#'     \item "2-stage-min-holm": Two-stage procedure with minimum Holm-adjusted p-value screening.
#'     \item "2-stage-cauchy": Two-stage procedure with Cauchy combination screening.
#'   }
#'
#' @return A data.frame with columns:
#'   \itemize{
#'     \item gene: Gene identifier
#'     \item cluster_id: Cluster identifier (varies by method)
#'     \item adj_p_val: Adjusted p-value (Benjamini-Hochberg correction)
#'     \item Additional columns depending on method (e.g., fisher_combined_p_val, min_p_val)
#'   }
#'
#' @export
#' @importFrom dplyr group_by summarize
#' @importFrom stats p.adjust pchisq
#'
#' @examples
#' \dontrun{
#' # Example with mock data
#' de_pvals_by_cluster <- list(
#'   cluster1 = data.frame(gene = c("A", "B"), cluster_id = "cluster1", p_val = c(0.01, 0.05)),
#'   cluster2 = data.frame(gene = c("A", "B"), cluster_id = "cluster2", p_val = c(0.02, 0.03))
#' )
#' result <- combine_de_pvals_by_cluster(de_pvals_by_cluster, method = "adjust_all")
#' }
combine_de_pvals_by_cluster <- function(de_pvals_by_cluster, method = "adjust_all") {
  method <- match.arg(method, choices = c("adjust_all", "min_per_cluster", "fisher_combined_and_min", "2-stage-fisher", "2-stage-min-holm", "2-stage-cauchy"))
  
  if (method == "adjust_all") {
    return(adjust_all_pvals(de_pvals_by_cluster))
  } else if (method == "min_per_cluster") {
    return(min_per_cluster_pvals(de_pvals_by_cluster))
  } else if (method == "fisher_combined_and_min") {
    return(fisher_combined_pvals_and_min(de_pvals_by_cluster))
  } else if (method == "2-stage-fisher") {
    print("Adjusting p-values using 2-stage procedure with fisher screening")
    res_table_to_return <- two_stage_adjustment(de_pvals_by_cluster, screen_method = "fisher")
    res_table_to_return$adj_p_val <- res_table_to_return$adj_p_val_2stage
    # cast res_table_to_return$gene to character
    res_table_to_return$gene <- as.character(res_table_to_return$gene)
    res_table_to_return$cluster_id <- as.character(res_table_to_return$cluster_id)
    return(res_table_to_return)
  } else if (method == "2-stage-min-holm") {
    print("Adjusting p-values using 2-stage procedure with min holm adjusted screening")
    res_table_to_return <- two_stage_adjustment(de_pvals_by_cluster, screen_method = "min_holm")
    res_table_to_return$adj_p_val <- res_table_to_return$adj_p_val_2stage
    # cast res_table_to_return$gene to character
    res_table_to_return$gene <- as.character(res_table_to_return$gene)
    res_table_to_return$cluster_id <- as.character(res_table_to_return$cluster_id)
    return(res_table_to_return)
  } else if (method == "2-stage-cauchy") {
    print("Adjusting p-values using 2-stage procedure with Cauchy screening")
    res_table_to_return <- two_stage_adjustment(de_pvals_by_cluster, screen_method = "cauchy")
        res_table_to_return$adj_p_val <- res_table_to_return$adj_p_val_2stage
    # cast res_table_to_return$gene to character
    res_table_to_return$gene <- as.character(res_table_to_return$gene)
    res_table_to_return$cluster_id <- as.character(res_table_to_return$cluster_id)
    return(res_table_to_return)
  }
}

#' Combine p-values and compute 2-stage adjusted p-values
#'
#' For each gene, combine p-values across clusters (Fisher, min Holm, or Cauchy),
#' then compute adjusted p-values from a 2-stage procedure using by-gene screening
#' (positive genes should have at least one cluster with a significant p-value).
#'
#' @param de_pvals_by_cluster A named list where each element corresponds to a
#'   cluster and contains a data.frame with columns 'gene', 'cluster_id', and 'p_val'.
#' @param screen_method Character. Screening combination: \code{"min_holm"} (default),
#'   \code{"fisher"}, or \code{"cauchy"}.
#' @param n_cores Integer or NULL. Number of cores for parallel index computation; NULL = serial.
#'
#' @return A data.frame with columns:
#'   \itemize{
#'     \item gene: Gene identifier
#'     \item cluster_id: Cluster with minimum p-values for this gene
#'     \item fisher_combined_p_val: Combined p-values from Fisher's method
#'     \item adj_p_val_2stage: Adjusted p-values from 2-stage procedure
#'     \item fisher_adj_pval: BH-adjusted p-values from Fisher's method
#'   }
#' @importFrom matrixStats rowSums2
#' @importFrom stats pchisq p.adjust pcauchy
#' @importFrom foreach %dopar% foreach
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @importFrom stats na.omit
#' @keywords internal
two_stage_adjustment <- function(de_pvals_by_cluster, screen_method = "min_holm", n_cores = NULL) {
  # Combine all p-values from clusters in a gene by cluster matrix
  # Get set of all genes
  all_genes <- unique(unlist(lapply(de_pvals_by_cluster, function(x) x$gene)))
  pval_matrix <- matrix(NA, nrow = length(all_genes), ncol = length(de_pvals_by_cluster))
  colnames(pval_matrix) <- names(de_pvals_by_cluster)
  rownames(pval_matrix) <- all_genes
  for (i in seq_along(de_pvals_by_cluster)) {
    clust_name <- names(de_pvals_by_cluster)[i]
    curr_clust_res_DE <- de_pvals_by_cluster[[clust_name]]
    pval_matrix[, i] <- curr_clust_res_DE$p_val[match(rownames(pval_matrix), curr_clust_res_DE$gene)]
  }
  if (screen_method == "min_holm") {
    # Apply Holm adjustment per cluster and take minimum adjusted p-value per gene across clusters
    combined_pvals <- apply(pval_matrix, 1, function(x) min(p.adjust(x, method = "holm"), na.rm = TRUE))
  } else if (screen_method == "fisher") {
       # Use fisher method to combine p-values per gene
      my_pvalues_fisher_method <- function(pvalues) {
        # TODO Add a check that all pvalues are "valid"
        pvalues[pvalues == 0] <- 1e-285
        lnp <- data.matrix(log(pvalues))
        chisq <- (-2) * matrixStats::rowSums2(lnp, na.rm = T)
        df <- 2 * rowSums(!is.na(pvalues))
        stats::pchisq(chisq, df, lower.tail = FALSE)
      }
      combined_pvals <- my_pvalues_fisher_method(pval_matrix)
  } else if (screen_method == "cauchy") {
      cauchyP <- function(p, w = 1/length(p)) {
        T <- tan((0.5 - p) * pi)
        Tsum <- sum(T * w)
        return(1 - pcauchy(Tsum))
      }
      combined_pvals <- apply(pval_matrix, 1, function(x) cauchyP(stats::na.omit(x)))
  } else {
    stop("Invalid screen_method specified. Must be one of 'min_holm', 'fisher', or 'cauchy'.")
  }
  screen_pvalues_adj <- p.adjust(combined_pvals, method = "BH")
  # Compute adjusted p-values following 2-stage procedure
  idx_pscreen <- function(a_matrix, ps_adj_sorted, n_cores = n_cores) {
    pval_flat <- as.vector(a_matrix)
    ps_R_sorted <- c(0, ps_adj_sorted) * seq(0, length(ps_adj_sorted))
    if (is.null(n_cores)) {
      idx_mat <- matrix(NA, nrow = nrow(a_matrix), ncol = ncol(a_matrix))
      for (idx in which(!is.na(pval_flat))) {
        p <- pval_flat[idx]
        i <- findInterval(p, ps_R_sorted, left.open = TRUE)
        idx_mat[idx] <- i
        # idx_mat[idx] <- which(p <= ps_R_sorted)[1] - 1
      }
      return(idx_mat)
    } else {
      cl <- parallel::makeCluster(n_cores)
      doParallel::registerDoParallel(cl)
      idx_list <- foreach::foreach(idx = which(!is.na(pval_flat)), .packages = c()) %dopar% {
        p <- pval_flat[idx]
        i <- findInterval(p, ps_R_sorted, left.open = TRUE)
        return(c(idx, i))
      }
      doParallel::stopImplicitCluster()
      idx_mat <- matrix(NA, nrow = nrow(a_matrix), ncol = ncol(a_matrix))
      for (res in idx_list) {
        idx_mat[res[1]] <- res[2]
      }
      return(idx_mat)
    }
  }
  nu_fun <- function(a_matrix, ps_adj) {
    ps_adj_sorted <- sort(ps_adj)
    idx_mat_star <- idx_pscreen(a_matrix, ps_adj_sorted, n_cores = n_cores)
    ps_mat_istar <- matrix(c(0, ps_adj_sorted, 1)[as.vector(idx_mat_star + 1)], nrow = nrow(a_matrix), ncol = ncol(a_matrix)) # +1 for R indexing
    quotient_mat <- a_matrix / (idx_mat_star - 1)
    nu_mat <- pmin(ps_mat_istar, quotient_mat)
    return(nu_mat)
  }
  p_adj_2_stage <- function(p_bar_gk, ps_adj) {
    G <- nrow(p_bar_gk)
    a_matrix <- p_bar_gk * G
    nu_mat <- nu_fun(a_matrix, ps_adj)
    pval_adj_final <- pmax(nu_mat, ps_adj)
    return(pval_adj_final)
  }
  pval_matrix_holm_adj <- t(apply(pval_matrix, 1, function(x) p.adjust(x, method = "holm")))
  mat_pcg_adj <- p_adj_2_stage(pval_matrix_holm_adj, screen_pvalues_adj)
  # print(head(mat_pcg_adj))
  rownames(mat_pcg_adj) <- rownames(pval_matrix)
  colnames(mat_pcg_adj) <- colnames(pval_matrix)
  res_table <- stats::na.omit(as.data.frame.table(mat_pcg_adj))
  names(res_table) <- c("gene", "cluster_id", "adj_p_val_2stage")
  res_table$screen_pval <- combined_pvals[res_table$gene]
  res_table$screen_adj_pval <- screen_pvalues_adj[res_table$gene]
  return(res_table)
}
