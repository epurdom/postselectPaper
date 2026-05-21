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
#' # Example with mock data
#' de_pvals_by_cluster <- list(
#'   cluster1 = data.frame(gene = c("A", "B"), cluster_id = "cluster1", p_val = c(0.01, 0.05)),
#'   cluster2 = data.frame(gene = c("A", "B"), cluster_id = "cluster2", p_val = c(0.02, 0.03))
#' )
#' result <- combine_de_pvals_by_cluster(de_pvals_by_cluster, method = "adjust_all")
combine_de_pvals_by_cluster <- function(de_pvals_by_cluster, method = "adjust_all") {
  method <- match.arg(method, choices = c("adjust_all", "min_per_cluster", "fisher_combined_and_min", "2-stage-fisher", "2-stage-min-holm", "2-stage-cauchy", "2-stage-simes"))
  
  if (method == "adjust_all") {
    return(adjust_all_pvals(de_pvals_by_cluster))
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
  } else if (method == "2-stage-simes") {
    print("Adjusting p-values using 2-stage procedure with Simes screening")
    res_table_to_return <- two_stage_adjustment(de_pvals_by_cluster, screen_method = "simes")
    res_table_to_return$adj_p_val <- res_table_to_return$adj_p_val_2stage
    # cast res_table_to_return$gene to character
    res_table_to_return$gene <- as.character(res_table_to_return$gene)
    res_table_to_return$cluster_id <- as.character(res_table_to_return$cluster_id)
    return(res_table_to_return)
  } 
  else {
    stop("Invalid method specified. Must be one of 'adjust_all', '2-stage-fisher', '2-stage-min-holm', or '2-stage-cauchy'.")
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
#' @param path2pb Character or NULL. Path to pseudobulk data for globaltest screening (if screen_method = "globaltest"). If NULL or file does not exist, globaltest screening will be skipped and NA returned for all adjusted p-values.
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
two_stage_adjustment <- function(de_pvals_by_cluster, screen_method = "min_holm", n_cores = NULL, path2pb = NULL) {
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
  } else if (screen_method == "stouffer") {
    # Use Stouffer method to combine p-values per gene
    my_pvalues_stouffer_method <- function(pvalues) {
      # TODO Add a check that all pvalues are "valid"
      pvalues[pvalues == 0] <- 1e-285
      z_scores <- stats::qnorm(pvalues)
      z_sum <- rowSums(z_scores, na.rm = TRUE)
      n <- rowSums(!is.na(pvalues))
      z_mean <- z_sum / sqrt(n)
      stats::pnorm(z_mean, lower.tail = FALSE)
    }
    combined_pvals <- my_pvalues_stouffer_method(pval_matrix)
  } else if (screen_method == "cauchy") {
      cauchyP <- function(p, w = 1/length(p)) {
        T <- tan((0.5 - p) * pi)
        Tsum <- sum(T * w)
        return(1 - pcauchy(Tsum))
      }
      combined_pvals <- apply(pval_matrix, 1, function(x) cauchyP(stats::na.omit(x)))
  } else if (screen_method == "simes") {
      simesP <- function(p) {
        m <- length(p)
        p_sorted <- sort(p)
        simes_pvals <- p_sorted * m / seq_along(p_sorted)
        return(min(simes_pvals))
      }
      combined_pvals <- apply(pval_matrix, 1, function(x) simesP(stats::na.omit(x)))
  } else if (screen_method == "globaltest") {
    if (is.null(path2pb) | !file.exists(path2pb)) {
      warning("Path to pseudobulk data not provided or file does not exist. Cannot perform globaltest screening. Returning NA for all adjusted p-values.")
      combined_pvals <- rep(NA, nrow(pval_matrix))
    } else {
      pseudobulk <- readRDS(path2pb)
      combined_pvals <- compute_gene_globaltest_pvalues(pseudobulk)
      combined_pvals <- combined_pvals[rownames(pval_matrix)]
    }
  }
  else {
    stop("Invalid screen_method specified. Must be one of 'min_holm', 'fisher', 'cauchy', 'stouffer', 'globaltest' or 'simes'.")
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

#' Compute globaltest p-values for each gene across cell types
#'
#' For each gene, perform a globaltest to assess association between gene expression
#' and phenotype across all cell types. Before that, it replicates the default gene and sample filtering
#' used in the muscat package, normalizes counts and log-transform.
#'
#' @param pb A SummarizedExperiment object containing pseudobulk data with assays for each cell type and a 'group_id' column in colData for phenotype.
#' @return A named numeric vector of globaltest p-values for each gene.
#' @importFrom edgeR DGEList calcNormFactors cpm filterByExpr
#' @importFrom SummarizedExperiment assay assayNames colData
#' @importFrom scater isOutlier
#' @importFrom globaltest gt p.value
compute_gene_globaltest_pvalues <- function(pb) {
  cell_types <- SummarizedExperiment::assayNames(pb)
  gene_names <- rownames(pb)
  group_id <- setNames(SummarizedExperiment::colData(pb)$group_id, SummarizedExperiment::colData(pb)$sample)
  phenotype <- group_id[colnames(pb)]
  n_samples <- ncol(pb)
  n_cells <- do.call(cbind, pb@int_colData$n_cells)
  # Filter out samples/clusters, genes, normalize counts and log-transform
  normalized_assays <- setNames(lapply(cell_types, function(ct) {
    y <- pb
    formula <- ~ group_id
    cd <- as.data.frame(SummarizedExperiment::colData(pb))
    design <- model.matrix(formula, cd)
    colnames(design) <- levels(SummarizedExperiment::colData(pb)$group_id)
    rmv <- n_cells[ct, ] < 10 # Filter out samples with fewer than 10 cells in this cluster
    y <- pb[ , !rmv]
    d <- design[colnames(y), , drop = FALSE]
    ls <- colSums(SummarizedExperiment::assay(y, ct))
    ol <- scater::isOutlier(ls, log = TRUE, type = "lower", nmads = 3)
    d <- d[colnames(y <- y[, !ol]), , drop = FALSE]
    y <- y[rowSums(SummarizedExperiment::assay(y, ct)) != 0, ]
    if (ncol(y) == 0 || nrow(y) == 0) {
      return(matrix(NA_real_, nrow = length(gene_names), ncol = n_samples,
        dimnames = list(gene_names, colnames(pb))))
    }
    if (max(SummarizedExperiment::assay(y, ct)) > 100) {
      keep <- edgeR::filterByExpr(SummarizedExperiment::assay(y, ct), d) # Filter out genes with low counts across samples using edgeR's filterByExpr function
      y <- y[keep, ]
    }

    counts_sub <- SummarizedExperiment::assay(y, ct)
    if (ncol(counts_sub) == 0 || nrow(counts_sub) == 0) {
      return(matrix(NA_real_, nrow = length(gene_names), ncol = n_samples,
        dimnames = list(gene_names, colnames(pb))))
    }

    dge_sub <- suppressMessages(edgeR::DGEList(counts_sub, remove.zeros = TRUE))
    dge_sub <- edgeR::calcNormFactors(dge_sub)
    cpm_sub <- edgeR::cpm(dge_sub, log = TRUE)

    mat_full <- matrix(NA_real_, nrow = length(gene_names), ncol = n_samples,
      dimnames = list(gene_names, colnames(pb)))
    if (ncol(cpm_sub) > 0) {
      mat_full[rownames(cpm_sub), colnames(cpm_sub)] <- cpm_sub
    }
    mat_full
  }), cell_types)

  pvals <- vapply(gene_names, function(g) {
    Y <- do.call(cbind, lapply(normalized_assays, function(mat) mat[g, ]))
    colnames(Y) <- cell_types
    rownames(Y) <- colnames(pb)
    # if (any(is.na(Y))) print(paste("NA values for gene", g))

    keep <- colSums(!is.na(Y)) > 0
    if (!any(keep)) return(NA_real_)
    Y <- Y[, keep, drop = FALSE]

    tryCatch({
      df_test <- data.frame(phenotype = phenotype, Y)
      gt_result <- globaltest::gt(phenotype ~ ., data = df_test)
      globaltest::p.value(gt_result)
    }, error = function(e) {
      message("Error in globaltest: ", e$message)
      NA_real_
    })
  }, numeric(1))

  names(pvals) <- gene_names
  pvals
}