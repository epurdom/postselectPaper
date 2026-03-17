

#' Prepare data for augmented simulation
#'
#' Prepares data for augmented simulation by processing a \code{SingleCellExperiment}
#' or using a pre-processed Seurat object. Performs normalization, dimensionality
#' reduction (PCA), and optional batch correction via Harmony. Optionally saves
#' the processed Seurat object to disk.
#'
#' @param seed Optional integer. Random seed for reproducibility. If \code{NULL}, seed is not set.
#' @param randomization Character. \code{"cells"} or \code{"samples"} for fake treatment/control assignment.
#' @param clustering Character. \code{"graph"}, \code{"kmeans"}, \code{"celltype"}, or \code{"leiden"} for gene-neighborhood clusters.
#' @param cut_at Numeric vector. Cluster sizes (e.g. \code{c(1, 3, 10)}). Use \code{NA} to get default \code{c(1, 3, 10)}. Ignored for \code{leiden}.
#' @param lfc_mean Numeric or vector. Mean log fold change for simulated DE genes. Use \code{NA} for default \code{1.5}.
#' @param n_hvgs Integer. Number of highly variable genes to keep (from real data). Required when \code{processed_sim} is \code{NULL}.
#' @param datapath Character. Path to input RDS (SingleCellExperiment). Required when \code{processed_sim} is \code{NULL}.
#' @param dest_file Optional character. Path to save the processed Seurat object. If \code{NULL}, object is not saved.
#' @param condition Optional character vector. Value(s) of \code{main_covariate} to filter cells; if \code{NULL}, all cells used.
#' @param use_harmony Logical. Whether to run Harmony batch correction.
#' @param leiden_resolution Optional numeric. Resolution for Leiden when \code{clustering == "leiden"}.
#' @param leiden_clusters Optional integer vector. Precomputed Leiden cluster labels when \code{clustering == "leiden"}.
#' @param main_covariate Character. Name of condition/treatment column in \code{colData}.
#' @param sample_covariate Character. Name of sample column (and for Harmony).
#' @param assay_continuous Character. Assay name used for variance (e.g. HVG selection).
#' @param cell_type_column Character. Name of cell type column (used when \code{clustering == "celltype"}).
#' @param batch_covariate Character. Name of batch column (for Harmony).
#' @param processed_sim Optional pre-processed Seurat object; if provided, \code{datapath} is not used.
#' @param num_pcs Integer. Number of principal components to compute (default \code{50}).
#' @param verbose Logical. If \code{TRUE}, print debug messages (default \code{TRUE}).
#'
#' @return A list with all input parameters (possibly updated) plus:
#'   \itemize{
#'     \item \code{pca_embeds}: Matrix of PCA (or Harmony) embeddings, cells in columns.
#'     \item \code{sample}: Sample ID per cell.
#'     \item \code{nc}: Number of cells.
#'     \item \code{sf}: Cell size factors (library size / median).
#'   }
#'
#' @importFrom Matrix colSums
#' @importFrom Seurat CreateSeuratObject AddMetaData NormalizeData FindVariableFeatures
#'   ScaleData RunPCA FindNeighbors FindClusters Idents Embeddings
#' @importFrom SingleCellExperiment counts
#' @importFrom SummarizedExperiment colData assay
#' @importFrom harmony RunHarmony
#' @importFrom MatrixGenerics rowVars
#' @importFrom stats kmeans rnorm rnbinom runif median
#' @importFrom methods is
#' @export
prep_aug_sim <- function(seed = NULL,
                        randomization,
                        clustering,
                        cut_at = NA_real_,
                        lfc_mean = NA_real_,
                        n_hvgs,
                        datapath = NULL,
                        dest_file = NULL,
                        condition = NULL,
                        use_harmony,
                        leiden_resolution = NULL,
                        leiden_clusters = NULL,
                        main_covariate,
                        sample_covariate,
                        assay_continuous,
                        cell_type_column,
                        batch_covariate,
                        processed_sim = NULL,
                        num_pcs = 50,
                        verbose = TRUE) {
  # Function to print debug info only when verbose is TRUE
  debug_print <- function(...) {
    if (verbose) {
      message(paste0("[DEBUG] ", ...))
    }
  }

  # defining random seed
  if (!is.null(seed)) {
    debug_print(paste("Setting random seed to:", seed))
    set.seed(seed)
  } else {
    debug_print("No seed provided, not setting random seed")
  }

  # Set randomization and clustering methods
  randomization_levels <- match.arg(randomization, c("cells", "samples"))
  debug_print(paste("Randomization level set to:", randomization_levels))
  clustering_method <- match.arg(clustering, c("graph", "kmeans", "celltype", "leiden"))
  debug_print(paste("Clustering method set to:", clustering_method))

  # Set default cluster sizes if not provided
  if (all(is.na(cut_at))) {
    cut_at <- c(1, 3, 10)
    debug_print(paste("Using default cut_at values:", paste(cut_at, collapse = ", ")))
  } else {
    debug_print(paste("Using provided cut_at values:", paste(cut_at, collapse = ", ")))
  }

  # Set default log fold change mean if not provided
  if (all(is.na(lfc_mean))) {
    lfc_mean <- 1.5
    debug_print(paste("Using default lfc_mean:", lfc_mean))
  } else {
    debug_print(paste("Using provided lfc_mean:", lfc_mean))
  }

  # ---------------------------------------

  if (is.null(processed_sim)) {
    debug_print("No pre-processed simulation data provided, loading from file")
    sce <- readRDS(datapath)
    debug_print(paste("Loading the file:", datapath))
    debug_print(paste("Dataset loaded with dimensions:", nrow(sce), "genes x", ncol(sce), "cells"))

    if (!is(sce, "SingleCellExperiment")) {
      stop("Loaded object is not a SingleCellExperiment")
    }

    # Filter to specified condition if provided
    if (!is.null(condition)) {
      debug_print(paste("Filtering to default condition:", condition))
      condition <- paste0(condition, collapse = " ")
      sce <- sce[, colData(sce)[[main_covariate]] == condition, drop = FALSE]
      debug_print(paste("After filtering:", nrow(sce), "genes x", ncol(sce), "cells"))
    }

    # Select highly variable genes if requested
    if (n_hvgs < nrow(sce)) {
      debug_print(paste("Selecting top", n_hvgs, "highly variable genes"))
      x <- assay(sce, assay_continuous)
      hvg <- order(-rowVars(x))
      n_hvgs <- min(nrow(sce), n_hvgs)
      sce <- sce[hvg[seq_len(n_hvgs)], ]
      debug_print(paste("After HVG selection:", nrow(sce), "genes x", ncol(sce), "cells"))
    }

    # Perform dimensionality reduction for clustering
    debug_print("Creating Seurat object")
    processed_sim <- CreateSeuratObject(counts = counts(sce))
    processed_sim <- AddMetaData(processed_sim, metadata = as.data.frame(colData(sce)))
    debug_print("Normalizing data")
    processed_sim <- NormalizeData(processed_sim)
    debug_print("Finding variable features")
    processed_sim <- FindVariableFeatures(processed_sim, nfeatures = n_hvgs)
    debug_print("Scaling data")
    processed_sim <- ScaleData(processed_sim)
    debug_print("Running PCA")
    processed_sim <- RunPCA(processed_sim, npcs = min(nrow(sce), ncol(sce), num_pcs), verbose = TRUE)

    if (clustering_method == "leiden") {
      if (is.null(leiden_clusters)) {
        debug_print("Running Leiden clustering")
        if (length(unique(cut_at)) > 1) {
          stop("Leiden clustering only supports one cut_at value at leiden resolution")
        }
        processed_sim <- FindNeighbors(processed_sim, dims = 1:num_pcs, reduction = "pca")
        processed_sim <- FindClusters(processed_sim, resolution = leiden_resolution)
        leiden_clusters <- Idents(processed_sim)
      } else {
        debug_print("Using Leiden clustering from INPUT")
      }
      num_leiden_clusters <- length(unique(leiden_clusters))
      debug_print(paste("Created", num_leiden_clusters, "clusters"))
      cut_at <- c(num_leiden_clusters)
      debug_print("Leiden clustering complete")
    }
    if (use_harmony) {
      debug_print("Running Harmony")
      processed_sim <- harmony::RunHarmony(processed_sim, group.by.vars = c(sample_covariate, batch_covariate))
    }
    if (!is.null(dest_file)) {
      debug_print(paste("Saving processed Seurat object to:", dest_file))
      saveRDS(processed_sim, dest_file)
    }
  } else {
    debug_print("Using pre-processed Seurat object")
    if (clustering_method == "leiden") {
      if (is.null(leiden_resolution)) {
        debug_print("Using Leiden clustering from saved Seurat object - Assuming the leiden resolution is already set")
        leiden_clusters <- Idents(processed_sim)
      } else {
        debug_print(paste("Using saved Leiden clustering with resolution:", leiden_resolution))
      }
      cut_at <- c(length(unique(leiden_clusters)))
      debug_print(paste("Assigned", cut_at[[1]], "Leiden Clusters"))
    }
  }

  if (use_harmony) {
    pca <- list(embedding = t(Embeddings(processed_sim, reduction = "harmony")))
  } else {
    pca <- list(embedding = t(Embeddings(processed_sim, reduction = "pca")))
  }

  mat_to_sf <- processed_sim@assays$RNA@layers$counts
  sf <- Matrix::colSums(mat_to_sf)
  sf <- sf / median(sf)

  list(
    seed = seed,
    randomization = randomization,
    clustering = clustering_method,
    cut_at = cut_at,
    lfc_mean = lfc_mean,
    n_hvgs = n_hvgs,
    datapath = datapath,
    dest_file = dest_file,
    condition = condition,
    use_harmony = use_harmony,
    leiden_resolution = leiden_resolution,
    leiden_clusters = leiden_clusters,
    main_covariate = main_covariate,
    sample_covariate = sample_covariate,
    assay_continuous = assay_continuous,
    cell_type_column = cell_type_column,
    batch_covariate = batch_covariate,
    pca_embeds = pca$embedding,
    sample = processed_sim@meta.data$sample,
    nc = ncol(processed_sim),
    sf = sf
  )
}

#' Select DE neighborhood (cluster assignment and DE cell selection)
#'
#' Internal helper used by \code{\link{create_augmented_data}}. Given a cluster
#' configuration, selects which cluster is the DE neighborhood (cells in that cluster will be differentially expressed).
#'
#' @param n_clusters Integer. Number of clusters in the clustering that will be used to generate the DE neighborhood
#' @param clustering_bank Named list of cluster assignment vectors for kmeans/leiden.
#' @param clust_to_sample_from Integer vector or \code{NULL}. Which cluster IDs to sample as DE cluster; \code{NULL} = all.
#' @param verbose Logical. If \code{TRUE}, print debug messages (default \code{TRUE}).
#'
#' @return A list with \code{n_clusters}, \code{sel_cluster}, \code{is_de_cell}.
#'
#' @keywords internal
select_de_neighborhood <- function(n_clusters,
                                     clustering_bank,
                                     clust_to_sample_from,
                                     verbose = TRUE) {
  if (verbose) {
    message(paste0("[DEBUG] Generating DE neighborhood with n_clusters =", n_clusters))
  }
  if (!(as.character(n_clusters) %in% names(clustering_bank))) {
    stop("Cluster not found in clustering bank")
  }
  cluster_assign <- clustering_bank[[as.character(n_clusters)]]

  if (!is.null(clust_to_sample_from)) {
    if (length(clust_to_sample_from) == 1) {
      sel_cluster <- clust_to_sample_from[[1]]
    } else {
      sel_cluster <- base::sample(clust_to_sample_from, size = 1)
    }
  } else {
    sel_cluster <- base::sample(unique(cluster_assign), size = 1)
  }
  if (verbose) message(paste("Selected cluster:", sel_cluster))
  is_de_cell <- cluster_assign == sel_cluster
  list(n_clusters = n_clusters, sel_cluster = sel_cluster, is_de_cell = is_de_cell)
}

#' Generate DE counts for a given DE neighborhood
#'
#' Internal helper used by \code{\link{create_augmented_data}}. Given a logical
#' vector \code{is_de_cell} (from \code{\link{select_de_neighborhood}}),
#' simulates counts for one gene that is differentially expressed between fake
#' treatment and control only in DE cells.
#'
#' @param is_de_cell Logical vector. \code{TRUE} for cells in the DE neighborhood.
#' @param base_expr Numeric. Base log2 expression level.
#' @param lfc_mean Numeric. Mean log fold change (treatment vs control) for this gene.
#' @param lfc_sd Numeric. SD of log fold change (default \code{0.5}).
#' @param sample_sd Numeric. SD of per-sample random effect (default \code{0.1}).
#' @param overdispersion Numeric. Negative binomial size parameter (default \code{0.1}).
#' @param randomization_by_cell Character vector. \code{"fake_ctrl"} or \code{"fake_trt"} per cell.
#' @param sample Character or factor. Sample ID per cell.
#' @param sf Numeric vector. Cell size factors.
#' @param nc Integer. Number of cells.
#' @param verbose Logical. If \code{TRUE}, print debug messages (default \code{TRUE}).
#'
#' @return A list with \code{lfc}, \code{base_expr}, \code{log_expression_level},
#'   \code{sample_mean}, \code{counts}.
#'
#' @keywords internal
generate_de_counts <- function(is_de_cell,
                              base_expr,
                              lfc_mean,
                              lfc_sd = 0.5,
                              sample_sd = 0.1,
                              overdispersion = 0.1,
                              randomization_by_cell,
                              sample,
                              sf,
                              nc,
                              verbose = TRUE) {
  eta_ctrl <- 0
  eta_trt <- stats::rnorm(1, mean = lfc_mean, sd = lfc_sd)
  trt_eff <- colSums(lemur:::one_hot_encoding(randomization_by_cell)[c("fake_ctrl", "fake_trt"), ] * c(eta_ctrl, eta_trt))

  sample_mean <- stats::rnorm(length(unique(sample)), mean = 0, sd = sample_sd)
  sample_eff <- colSums(lemur:::one_hot_encoding(sample) * sample_mean)

  mu <- 2^(trt_eff * is_de_cell + sample_eff + base_expr) * sf
  counts <- stats::rnbinom(n = nc, mu = mu, size = 1 / overdispersion)

  list(
    lfc = eta_trt - eta_ctrl,
    base_expr = base_expr,
    log_expression_level = unname(trt_eff * is_de_cell + sample_eff + base_expr),
    sample_mean = sample_mean,
    counts = counts
  )
}

#' Generate one simulated DE gene for a given cluster configuration
#'
#' Internal helper used by \code{\link{create_augmented_data}}. Generates counts
#' for a single gene that is differentially expressed between fake treatment and
#' control only within one selected cell cluster. This is a convenience wrapper
#' that calls \code{\link{select_de_neighborhood}} then \code{\link{generate_de_counts}}.
#'
#' @param n_clusters Integer. Number of clusters (or cluster index) for this gene.
#' @param base_expr Numeric. Base log2 expression level.
#' @param lfc_mean Numeric. Mean log fold change (treatment vs control) for this gene.
#' @param lfc_sd Numeric. SD of log fold change (default \code{0.5}).
#' @param sample_sd Numeric. SD of per-sample random effect (default \code{0.1}).
#' @param overdispersion Numeric. Negative binomial size parameter (default \code{0.1}).
#' @param clustering_bank Named list of cluster assignment vectors for kmeans/leiden.
#' @param clust_to_sample_from Integer vector or \code{NULL}. Which cluster IDs to sample as DE cluster; \code{NULL} = all.
#' @param randomization_by_cell Character vector. \code{"fake_ctrl"} or \code{"fake_trt"} per cell.
#' @param sample Character or factor. Sample ID per cell.
#' @param sf Numeric vector. Cell size factors.
#' @param nc Integer. Number of cells.
#' @param verbose Logical. If \code{TRUE}, print debug messages (default \code{TRUE}).
#'
#' @return A list with \code{n_clusters}, \code{sel_cluster}, \code{is_de_cell}, \code{lfc},
#'   \code{base_expr}, \code{log_expression_level}, \code{sample_mean}, \code{counts}.
#'
#' @keywords internal
generate_gene_for_cluster <- function(n_clusters,
                                     base_expr,
                                     lfc_mean,
                                     lfc_sd = 0.5,
                                     sample_sd = 0.1,
                                     overdispersion = 0.1,
                                     clustering_bank,
                                     clust_to_sample_from,
                                     randomization_by_cell,
                                     sample,
                                     sf,
                                     nc,
                                     verbose = TRUE) {
  nbhd <- select_de_neighborhood(
    n_clusters = n_clusters,
    clustering_bank = clustering_bank,
    clust_to_sample_from = clust_to_sample_from,
    verbose = verbose
  )
  cnt <- generate_de_counts(
    is_de_cell = nbhd$is_de_cell,
    base_expr = base_expr,
    lfc_mean = lfc_mean,
    lfc_sd = lfc_sd,
    sample_sd = sample_sd,
    overdispersion = overdispersion,
    randomization_by_cell = randomization_by_cell,
    sample = sample,
    sf = sf,
    nc = nc,
    verbose = verbose
  )
  c(nbhd, cnt)
}

#' Create augmented data with simulated differentially expressed genes
#'
#' Builds fake treatment/control groups and generates synthetic genes that are
#' differentially expressed between conditions only within selected cell
#' clusters (e.g. k-means, graph, Leiden, or cell type). Uses the embedding and
#' parameters produced by \code{\link{prep_aug_sim}}.
#'
#' @param clustering Character. \code{"graph"}, \code{"kmeans"}, \code{"leiden"}, or \code{"celltype"}.
#' @param pca_embeds Matrix. PCA (or Harmony) embeddings, cells in columns (from \code{prep_aug_sim}).
#' @param sample Character or factor. Sample ID per cell (from \code{prep_aug_sim}).
#' @param nc Integer. Number of cells (from \code{prep_aug_sim}).
#' @param sf Numeric vector. Cell size factors (from \code{prep_aug_sim}).
#' @param randomization Character. \code{"cells"} or \code{"samples"} for fake treatment/control assignment.
#' @param cut_at Numeric vector. Cluster sizes used for each simulated gene (recycled to \code{n_de_genes}).
#' @param leiden_clusters Optional integer vector. Cluster labels when \code{clustering == "leiden"}.
#' @param clust_to_sample_from Optional integer vector. Which cluster IDs to sample as DE clusters; if \code{NULL}, all clusters.
#' @param n_de_genes Integer. Number of simulated DE genes to generate.
#' @param lfc_mean Numeric or vector. Mean log fold change (recycled; sign randomized per gene).
#' @param verbose Logical. If \code{TRUE}, print progress and per-gene debug info (default \code{TRUE}).
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{new_counts}: Matrix of simulated counts (\code{n_de_genes} x \code{nc}).
#'     \item \code{clustering_bank}: Named list of cluster vectors (for \code{kmeans}/\code{leiden}); used for downstream alignment.
#'     \item \code{de_args}: A tibble with one row per simulated gene: \code{name}, \code{is_simulated}, \code{cut_at}, \code{base_expr}, \code{lfc}, \code{sel_cluster}, \code{is_de_cell}.
#'     \item \code{fake_condition}: Character vector of \code{"fake_ctrl"} or \code{"fake_trt"} per cell.
#'   }
#'
#' @importFrom bluster makeKNNGraph
#' @importFrom BiocNeighbors AnnoyParam
#' @importFrom igraph cluster_walktrap membership cut_at
#' @importFrom tibble tibble
#' @importFrom purrr map
#' @export
create_augmented_data <- function(clustering,
                                  pca_embeds,
                                  sample,
                                  nc,
                                  sf,
                                  randomization,
                                  cut_at,
                                  leiden_clusters = NULL,
                                  clust_to_sample_from = NULL,
                                  n_de_genes,
                                  lfc_mean,
                                  verbose = TRUE) {
  debug_print <- function(...) {
    if (verbose) {
      message(paste0("[DEBUG] ", ...))
    }
  }

  clustering_method <- clustering
  sample_assignments <- sample
  num_cells_data <- nc
  randomization_levels <- randomization

  debug_print(paste("Creating fake treatment/control groups with randomization at", randomization, "level"))
  randomization_by_cell <- NULL
  if (randomization_levels == "cells") {
    randomization_by_cell <- sample(c("fake_ctrl", "fake_trt"), size = num_cells_data, replace = TRUE)
    debug_print(paste("Assigned", sum(randomization_by_cell == "fake_ctrl"), "cells to fake_ctrl and",
                    sum(randomization_by_cell == "fake_trt"), "cells to fake_trt"))
  } else if (randomization_levels == "samples") {
    samples <- unique(sample_assignments)
    trt_samples <- sample(samples, size = round(length(samples) / 2), replace = FALSE)
    randomization_by_cell <- ifelse(sample_assignments %in% trt_samples, "fake_trt", "fake_ctrl")
    debug_print(paste("Assigned", length(trt_samples), "samples to fake_trt and",
                    length(samples) - length(trt_samples), "samples to fake_ctrl"))
    debug_print(paste("This results in", sum(randomization_by_cell == "fake_ctrl"), "cells as fake_ctrl and",
                    sum(randomization_by_cell == "fake_trt"), "cells as fake_trt"))
  } else {
    stop("Invalid randomization level")
  }

  debug_print(paste("Creating clusters using method:", clustering_method))
  clustering_bank <- list()
  if (clustering_method == "kmeans") {
    debug_print(paste("Running kmeans clustering with k values:", paste(unique(cut_at), collapse = ", ")))
    clustering_bank <- lapply(unique(cut_at), function(k) {
      debug_print(paste("  Running kmeans with k =", k))
      stats::kmeans(t(pca_embeds), centers = k, iter.max = 100)$cluster
    })
    names(clustering_bank) <- as.character(unique(cut_at))
    debug_print("Kmeans clustering complete")
  } else if (clustering_method == "leiden") {
    debug_print("Using Leiden clustering")
    if (length(unique(cut_at)) > 1) {
      stop("Leiden clustering only supports one cut_at value at leiden resolution")
    }
    clustering_bank <- list()
    clustering_bank[[as.character(cut_at[[1]])]] <- leiden_clusters
    debug_print(paste("Assigned Leiden Clustering to", cut_at[[1]], "clusters"))
  }
  if (verbose) message(paste("List to sample from:", clust_to_sample_from))

  de_context <- list(
    clustering_method = clustering_method,
    clustering_bank = clustering_bank,
    clust_to_sample_from = clust_to_sample_from,
    randomization_by_cell = randomization_by_cell,
    sample = sample,
    sf = sf,
    nc = nc,
    verbose = verbose
  )

  debug_print(paste("Setting up parameters for", n_de_genes, "simulated DE genes"))
  de_args <- tibble(
    name = paste0("simulated_gene-", seq_len(n_de_genes)),
    is_simulated = TRUE,
    cut_at = rep_len(cut_at, length.out = n_de_genes),
    base_expr = stats::runif(n = n_de_genes, min = -7, max = 3),
    lfc = rep_len(lfc_mean, length.out = n_de_genes) * sample(c(-1, 1), size = n_de_genes, replace = TRUE),
    sel_cluster = rep(-1, n_de_genes),
    is_de_cell = purrr::map(seq_len(n_de_genes), function(.) rep(FALSE, nc))
  )
  debug_print("DE gene parameters initialized")

  debug_print("Generating simulated counts for DE genes")
  new_counts <- matrix(0, nrow = n_de_genes, ncol = nc)
  for (idx in seq_len(n_de_genes)) {
    args_gene <- list(
      n_clusters = de_args$cut_at[idx],
      base_expr = de_args$base_expr[idx],
      lfc_mean = de_args$lfc[idx],
      lfc_sd = 0,
      sample_sd = 0.1,
      overdispersion = 0.2
    )
    gv <- do.call(generate_gene_for_cluster, c(args_gene, de_context))
    new_counts[idx, ] <- gv$counts
    de_args$is_de_cell[[idx]] <- gv$is_de_cell
    de_args$sel_cluster[idx] <- gv$sel_cluster

    if (verbose) {
      up_cells_ctrl <- sum(gv$is_de_cell & randomization_by_cell == "fake_ctrl")
      up_cells_trt <- sum(gv$is_de_cell & randomization_by_cell == "fake_trt")
      debug_print(paste("  DE cells:", sum(gv$is_de_cell),
                       "(", up_cells_ctrl, "in ctrl,",
                       up_cells_trt, "in trt)"))
      debug_print(paste("  Gene mean count:", mean(gv$counts), "max:", max(gv$counts)))
    }

    if (idx %% 20 == 0) {
      debug_print(paste("Progress:", round(idx / n_de_genes * 100), "% complete"))
    }
  }
  debug_print("Finished generating all DE genes")
  list(de_args = de_args, new_counts = new_counts, clustering_bank = clustering_bank, fake_condition = randomization_by_cell)
}
