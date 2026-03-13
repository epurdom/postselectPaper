

#' Prepare data for augmented simulation
#'
#' Prepares data for augmented simulation by processing a \code{SingleCellExperiment}
#' or using a pre-processed Seurat object. Performs normalization, dimensionality
#' reduction (PCA), and optional batch correction via Harmony. Optionally saves
#' the processed Seurat object to disk.
#'
#' @param pa A list of simulation parameters. Required/optional elements:
#'   \itemize{
#'     \item \code{seed}: (optional) Random seed for reproducibility.
#'     \item \code{randomization}: \code{"cells"} or \code{"samples"} for fake treatment/control assignment.
#'     \item \code{clustering}: \code{"graph"}, \code{"kmeans"}, \code{"celltype"}, or \code{"leiden"} for gene-neighborhood clusters.
#'     \item \code{cut_at}: Numeric vector of cluster sizes (e.g. \code{c(1, 3, 10)}). Ignored for \code{leiden}.
#'     \item \code{lfc_mean}: Numeric (or vector) mean log fold change for simulated DE genes.
#'     \item \code{n_hvgs}: Number of highly variable genes to keep (from real data).
#'     \item \code{datapath}: Path to input file (RDS SingleCellExperiment) when \code{processed_sim} is \code{NULL}.
#'     \item \code{dest_file}: (optional) Path to save the processed Seurat object.
#'     \item \code{condition}: (optional) Value(s) of \code{config$main_covariate} to filter cells; if \code{NULL}, all cells used.
#'     \item \code{use_harmony}: Whether to run Harmony batch correction.
#'     \item \code{leiden_resolution}: (optional) Resolution for Leiden when \code{clustering == "leiden"}.
#'     \item \code{leiden_clusters}: (optional) Precomputed Leiden cluster labels when \code{clustering == "leiden"}.
#'   }
#' @param config A list of data configuration:
#'   \itemize{
#'     \item \code{main_covariate}: Name of condition/treatment column in \code{colData}.
#'     \item \code{sample_covariate}: Name of sample column.
#'     \item \code{assay_continuous}: Assay name used for variance (e.g. HVG selection).
#'     \item \code{cell_type_column}: Name of cell type column (used when \code{clustering == "celltype"}).
#'     \item \code{batch_covariate}: Name of batch column (for Harmony).
#'   }
#' @param processed_sim Optional pre-processed Seurat object; if provided, \code{datapath} is not used.
#' @param num_pcs Number of principal components to compute (default \code{50}).
#' @param verbose If \code{TRUE}, print debug messages (default \code{TRUE}).
#'
#' @return The same list \code{pa} with additional elements set:
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
#' @export
prep_aug_sim <- function(pa, config, processed_sim = NULL, num_pcs = 50, verbose = TRUE){
  # Function to print debug info only when verbose is TRUE
  debug_print <- function(...) {
    if (verbose) {
      message(paste0("[DEBUG] ", ...))
    }
  }

  # definining random seed
  if (!is.null(pa$seed)) {
    debug_print(paste("Setting random seed to:", pa$seed))
    set.seed(pa$seed)  # Set random seed for reproducibility
  } else {
    debug_print("No seed provided, not setting random seed")
  }

    # Set randomization and clustering methods
  randomization_levels <- match.arg(pa$randomization, c("cells", "samples"))
  debug_print(paste("Randomization level set to:", randomization_levels))
  clustering_method <- match.arg(pa$clustering, c("graph", "kmeans", "celltype", "leiden"))
  debug_print(paste("Clustering method set to:", clustering_method))

    # Set default cluster sizes if not provided
  if(all(is.na(pa$cut_at))){
    pa$cut_at <- c(1, 3, 10)  # Default cluster sizes if not provided
    debug_print(paste("Using default cut_at values:", paste(pa$cut_at, collapse=", ")))
  } else {
    debug_print(paste("Using provided cut_at values:", paste(pa$cut_at, collapse=", ")))
  }

  # Set default log fold change mean if not provided
  if(all(is.na(pa$lfc_mean))){
    pa$lfc_mean <- 1.5  # Default log fold change mean if not provided
    debug_print(paste("Using default lfc_mean:", pa$lfc_mean))
  } else {
    debug_print(paste("Using provided lfc_mean:", pa$lfc_mean))
  }
  

  # ---------------------------------------

  if(is.null(processed_sim)){
    debug_print("No pre-processed simulation data provided, loading from file")
    # sce <- get(load(pa$datapath))
    sce <- readRDS(pa$datapath)
    debug_print(paste("Loading the file:", pa$datapath))
    debug_print(paste("Dataset loaded with dimensions:", nrow(sce), "genes x", ncol(sce), "cells"))

    if (!is(sce, "SingleCellExperiment")) {
      stop("Loaded object is not a SingleCellExperiment")
    }

    # Filter to specified condition if provided
    if(!is.null(pa$condition)){
      debug_print(paste("Filtering to default condition:", pa$condition))
      pa$condition <- paste0(pa$condition, collapse = " ")
      sce <- sce[,colData(sce)[[config$main_covariate]] == pa$condition]  # Filter to specified condition
      debug_print(paste("After filtering:", nrow(sce), "genes x", ncol(sce), "cells"))
    }

    # Select highly variable genes if requested
    if(pa$n_hvgs < nrow(sce)){
      debug_print(paste("Selecting top", pa$n_hvgs, "highly variable genes"))
      x <- assay(sce, config$assay_continuous)
      hvg <- order(-rowVars(x))  # Order genes by variance (MatrixGenerics dispatches on assay type)
      n_hvgs <- min(nrow(sce), pa$n_hvgs)
      sce <- sce[hvg[seq_len(n_hvgs)],]  # Subset to top HVGs
      debug_print(paste("After HVG selection:", nrow(sce), "genes x", ncol(sce), "cells"))
    }

    # Perform dimensionality reduction for clustering
    # Seurat version
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

    if (pa$clustering == "leiden") {
      if (is.null(pa$leiden_clusters)){
        debug_print("Running Leiden clustering")
        # assert that there is only one cut
        if(length(unique(pa$cut_at)) > 1){
          stop("Leiden clustering only supports one cut_at value at leiden resolution")
        }
        processed_sim <- FindNeighbors(processed_sim, dims = 1:num_pcs, reduction = "pca")
        processed_sim <- FindClusters(processed_sim, resolution = pa$leiden_resolution)
        pa$leiden_clusters <- Idents(processed_sim)
      } else {
        debug_print("Using Leiden clustering from INPUT pa object")
      }
      num_leiden_clusters <- length(unique(pa$leiden_clusters))
      debug_print(paste("Created", num_leiden_clusters, "clusters"))
      # change the cut_at to the number of leiden clusters
      pa$cut_at <- c(num_leiden_clusters)
      debug_print("Leiden clustering complete")
    }
    # Run Harmony batch correction if requested
    if (pa$use_harmony) {
      debug_print("Running Harmony")
      processed_sim <- harmony::RunHarmony(processed_sim, group.by.vars = c(config$sample_covariate, config$batch_covariate))
    }
    if(!is.null(pa$dest_file)){
      debug_print(paste("Saving processed Seurat object to:", pa$dest_file))
      saveRDS(processed_sim, pa$dest_file)
    }
  } else {
    debug_print("Using pre-processed Seurat object")
    if (pa$clustering == "leiden") { # if the clustering method is leiden, we assume the leiden resolution is already set
      if (is.null(pa$leiden_resolution)){
        debug_print("Using Leiden clustering from saved Seurat object - Assuming the leiden resolution is already set")
        pa$leiden_clusters <- Idents(processed_sim)
      } else {
        debug_print(paste("Using saved Leiden clustering with resolution:", pa$leiden_resolution))
      }
        pa$cut_at <- c(length(unique(pa$leiden_clusters)))
        debug_print(paste("Assigned", pa$cut_at[[1]], "Leiden Clusters"))
    }
  }

  if (pa$use_harmony) {
    pca <- list(embedding = t(Embeddings(processed_sim, reduction = "harmony")))
  } else {
    pca <- list(embedding = t(Embeddings(processed_sim, reduction = "pca")))
  }


  # create size factors
  mat_to_sf <- processed_sim@assays$RNA@layers$counts
  sf <- Matrix::colSums(mat_to_sf)
  sf <- sf / median(sf)


  # store parameters
  pa$pca_embeds <- pca$embedding
  pa$sample <- processed_sim@meta.data$sample
  pa$nc <- ncol(processed_sim)
  pa$sf <- sf

  return(pa)
}

#' Create augmented data with simulated differentially expressed genes
#'
#' Builds fake treatment/control groups and generates synthetic genes that are
#' differentially expressed between conditions only within selected cell
#' clusters (e.g. k-means, graph, Leiden, or cell type). Uses the embedding and
#' parameters produced by \code{\link{prep_aug_sim}}.
#'
#' @param pa Parameter list from \code{\link{prep_aug_sim}}. Must include:
#'   \itemize{
#'     \item \code{pca_embeds}, \code{sample}, \code{nc}, \code{sf}: from \code{prep_aug_sim}.
#'     \item \code{randomization}: \code{"cells"} or \code{"samples"}.
#'     \item \code{clustering}: \code{"graph"}, \code{"kmeans"}, \code{"leiden"}, or \code{"celltype"}.
#'     \item \code{cut_at}: Cluster sizes used for each simulated gene (recycled to \code{n_de_genes}).
#'     \item \code{n_de_genes}: Number of simulated DE genes to generate.
#'     \item \code{lfc_mean}: Mean log fold change (recycled; sign randomized per gene).
#'     \item \code{leiden_clusters}: (optional) Cluster labels when \code{clustering == "leiden"}.
#'     \item \code{clust_to_sample_from}: (optional) Which cluster IDs to sample as DE clusters; if \code{NULL}, all clusters.
#'   }
#' @param config Config list; \code{cell_type_column} is used when \code{clustering == "celltype"}.
#' @param verbose If \code{TRUE}, print progress and per-gene debug info (default \code{TRUE}).
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{new_counts}: Matrix of simulated counts (\code{n_de_genes} x \code{nc}).
#'     \item \code{kmeans_clusterings}: Named list of cluster vectors (for \code{kmeans}/\code{leiden}); used for downstream alignment.
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
create_augmented_data <- function(pa, config, verbose = TRUE){
  # load the processed sim
  debug_print <- function(...) {
    if (verbose) {
      message(paste0("[DEBUG] ", ...))
    }
  }

  clustering_method <- pa$clustering
  pca_embeds <- pa$pca_embeds
  sample_assignments <- pa$sample
  num_cells_data <- pa$nc
  randomization_levels <- pa$randomization

  debug_print(paste("Creating fake treatment/control groups with randomization at", pa$randomization, "level"))
  randomization_by_cell <- NULL
  if(randomization_levels == "cells"){
    # Randomize at cell level
    randomization_by_cell <- sample(c("fake_ctrl", "fake_trt"), size = num_cells_data, replace = TRUE)
    debug_print(paste("Assigned", sum(randomization_by_cell == "fake_ctrl"), "cells to fake_ctrl and", 
                    sum(randomization_by_cell == "fake_trt"), "cells to fake_trt"))
  }else if (randomization_levels == "samples"){
    # Randomize at sample level
    samples <- unique(sample_assignments)
    # trt_samples <- sample(samples, size = round(length(samples)), replace = FALSE) ## original code seems wrong????
    trt_samples <- sample(samples, size = round(length(samples)/2), replace = FALSE) # balanced samples for ctrl and trt
    
    randomization_by_cell <- ifelse(sample_assignments %in% trt_samples, "fake_trt", "fake_ctrl")
    debug_print(paste("Assigned", length(trt_samples), "samples to fake_trt and", 
                    length(samples) - length(trt_samples), "samples to fake_ctrl"))
    debug_print(paste("This results in", sum(randomization_by_cell== "fake_ctrl"), "cells as fake_ctrl and", 
                    sum(randomization_by_cell == "fake_trt"), "cells as fake_trt"))
  } else {
    stop("Invalid randomization level")
  }

  debug_print(paste("Creating clusters using method:", clustering_method))
  if(clustering_method == "graph"){
    # Graph-based clustering using k-nearest neighbors and walktrap algorithm
    debug_print("Building KNN graph")
    graph <- bluster::makeKNNGraph(t(pca_embeds), k = 15, BNPARAM = BiocNeighbors::AnnoyParam())
    debug_print("Running walktrap clustering")
    clustering <- igraph::cluster_walktrap(graph)
    debug_print(paste("Created", length(unique(igraph::membership(clustering))), "clusters"))
  }else if(clustering_method == "kmeans"){
    # K-means clustering with centers specified by cut_at values
    debug_print(paste("Running kmeans clustering with k values:", paste(unique(pa$cut_at), collapse=", ")))
    kmeans_clusterings <- lapply(unique(pa$cut_at), function(k){
      debug_print(paste("  Running kmeans with k =", k))
      kmeans(t(pca_embeds), centers = k, iter.max = 100)$cluster
    })
    names(kmeans_clusterings) <- as.character(unique(pa$cut_at))
    debug_print("Kmeans clustering complete")
  }else if(clustering_method == "leiden"){
    debug_print("Using Leiden clustering")
    # assert that there is only one cut
    if(length(unique(pa$cut_at)) > 1){
      stop("Leiden clustering only supports one cut_at value at leiden resolution")
    }
    kmeans_clusterings <- list()
    kmeans_clusterings[[as.character(pa$cut_at[[1]])]] <- pa$leiden_clusters
    debug_print(paste("Assigned Leiden Clustering to", pa$cut_at[[1]], "clusters"))
  }else if(clustering_method == "celltype"){ 
    # Use existing cell type annotations
    debug_print("Using existing cell type annotations")
    celltypes <- colData(sce)[[config$cell_type_column]]
    celltypes[is.na(celltypes)] <- "MISSING_ANNOTATION"  # Handle missing cell types
    celltype_counts <- table(celltypes)
    debug_print(paste("Found", length(unique(celltypes)), "cell types"))
    if(verbose) {
      debug_print("Cell type distribution:")
      for(i in 1:length(celltype_counts)) {
        debug_print(paste("  ", names(celltype_counts)[i], ":", celltype_counts[i], "cells"))
      }
    }
  }
  print(paste("List to sample from:", pa$clust_to_sample_from))

  # Function to generate simulated genes with DE in specific clusters
  # This creates genes that are differentially expressed between conditions,
  # but only within specific cell clusters (targeted DE)
  generate_gene_for_cluster <- function(n_clusters = 10, base_expr = -2, lfc_mean = 1, lfc_sd = 0.5, sample_sd = 0.1, overdispersion = 0.1, ...){
    debug_print(paste("Generating DE gene with n_clusters =", n_clusters, "and lfc_mean =", lfc_mean))
    
    # Determine cell cluster assignments based on selected clustering method
    if(clustering_method == "kmeans"){
      print(paste("Using kmeans clustering for", n_clusters, "clusters"))
      if(as.character(n_clusters) %in% names(kmeans_clusterings)){
        print(paste("Using precomputed kmeans clustering for", n_clusters, "clusters"))
      }else{
        print(paste("Computing kmeans clustering for", n_clusters, "clusters"))
      }
    }
    cluster_assign <- if(clustering_method == "graph"){
      igraph::cut_at(clustering, no = n_clusters)
    }else if(clustering_method == "kmeans" || clustering_method == "leiden"){
      if(as.character(n_clusters) %in% names(kmeans_clusterings)){
        kmeans_clusterings[[as.character(n_clusters)]]
      }else{
        if ("clustering_method" == "leiden"){
          stop("Leiden clustering not being properly used for DE gene generation")
        }
        kmeans(t(pca_embeds), centers = n_clusters)$cluster
      }
    }else if(clustering_method == "celltype"){
      celltypes
    }
    
    # Randomly select one cluster to have DE
    # if clust_to_sample_from is specified, use that cluster
    if (!is.null(pa$clust_to_sample_from)) {
      print(paste("Sampling from clusters:", pa$clust_to_sample_from))
      if (length(pa$clust_to_sample_from) == 1) {
        sel_cluster <- pa$clust_to_sample_from[[1]]
      } else {
        sel_cluster <- sample(pa$clust_to_sample_from, size = 1)
      }
    } else {
      print(paste("Sampling from all clusters"))
      sel_cluster <- sample(unique(cluster_assign), size = 1)
    }
    print(paste("Selected cluster:", sel_cluster))
    is_de_cell <- cluster_assign == sel_cluster
    
    # Set up log fold change between conditions
    eta_ctrl <- 0
    eta_trt <- rnorm(1, mean = lfc_mean, sd = lfc_sd)
    
    # Calculate cell-specific effects:
    # 1. Treatment effect (only applied to cells in the selected cluster & fake_trt)
    trt_eff <- colSums(lemur:::one_hot_encoding(randomization_by_cell)[c("fake_ctrl", "fake_trt"),] * c(eta_ctrl, eta_trt))

    mouse_mean <- rnorm(length(unique(pa$sample)), mean = 0, sd = sample_sd) # create effects for each sample effects
    mouse_eff <- colSums(lemur:::one_hot_encoding(pa$sample) * mouse_mean) # add sample effects to each cell
    # 3. Cell size normalization
    # mat_to_sf <- counts(sce)
    # sf <- Matrix::colSums(mat_to_sf)
    # sf <- sf / median(sf)
    
    # Calculate expected counts using negative binomial model
    mu <- 2^(trt_eff * is_de_cell + mouse_eff + base_expr) * pa$sf
    counts <- rnbinom(n = pa$nc, mu = mu, size = 1/overdispersion)

    # Return simulation details for tracking
    list(n_clusters = n_clusters, sel_cluster = sel_cluster, is_de_cell = is_de_cell, 
         lfc = eta_trt - eta_ctrl, base_expr = base_expr,
         log_expression_level = unname(trt_eff * is_de_cell + mouse_eff + base_expr),
         mouse_mean = mouse_mean, counts = counts)
  }

  # Set up parameters for simulated DE genes
  debug_print(paste("Setting up parameters for", pa$n_de_genes, "simulated DE genes"))
  de_args <- tibble(name = paste0("simulated_gene-", seq_len(pa$n_de_genes)),
                    is_simulated =  TRUE,
                    cut_at = rep_len(pa$cut_at, length.out = pa$n_de_genes),
                    base_expr = runif(n = pa$n_de_genes, min = -7, max = 3),  # Random base expression
                    lfc = rep_len(pa$lfc_mean, length.out = pa$n_de_genes) * sample(c(-1, 1), size = pa$n_de_genes, replace = TRUE),  # Log fold changes
                    sel_cluster = rep(-1, pa$n_de_genes),
                    is_de_cell = purrr::map(seq_len(pa$n_de_genes), function(.) rep(FALSE, pa$nc)))  # Initialize DE cell indicators
  debug_print("DE gene parameters initialized")


  # Generate simulated counts for DE genes
  debug_print("Generating simulated counts for DE genes")
  new_counts <- matrix(0, nrow = pa$n_de_genes, ncol = pa$nc)
  for(idx in seq_len(pa$n_de_genes)){
    # Generate each simulated gene with cluster-specific DE
    gv <- generate_gene_for_cluster(n_clusters = de_args$cut_at[idx], 
                                    base_expr = de_args$base_expr[idx],
                                    lfc_mean = de_args$lfc[idx],
                                    lfc_sd = 0, sample_sd = 0.1, overdispersion = 0.2)  
    new_counts[idx,] <- gv$counts  # Store simulated counts
    de_args$is_de_cell[[idx]] <- gv$is_de_cell  # Store which cells should have DE
    de_args$sel_cluster[idx] <- gv$sel_cluster  # Store the selected cluster
    
    # Debug info about this gene
    if(verbose) {
      up_cells_ctrl <- sum(gv$is_de_cell & randomization_by_cell == "fake_ctrl")
      up_cells_trt <- sum(gv$is_de_cell & randomization_by_cell == "fake_trt")
      debug_print(paste("  DE cells:", sum(gv$is_de_cell), 
                       "(", up_cells_ctrl, "in ctrl,", 
                       up_cells_trt, "in trt)"))
      debug_print(paste("  Gene mean count:", mean(gv$counts), "max:", max(gv$counts)))
    }
    
    # Progress reporting for long runs
    if(idx %% 20 == 0) {
      debug_print(paste("Progress:", round(idx/pa$n_de_genes*100), "% complete"))
    }
  }
  debug_print("Finished generating all DE genes")
  return(list(de_args = de_args, new_counts = new_counts, kmeans_clusterings = kmeans_clusterings, fake_condition = randomization_by_cell))
}
