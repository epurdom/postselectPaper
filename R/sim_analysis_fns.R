#' Perform Leiden clustering on a Seurat object
#'
#' Runs graph-based clustering (FindNeighbors + FindClusters) on a specified
#' dimensionality reduction and returns factor labels \code{cluster1}, \code{cluster2}, etc.
#'
#' @param seurat_obj A Seurat object with the requested reduction already computed.
#' @param leiden_res Resolution for \code{FindClusters} (default \code{0.1}).
#' @param reduction_type Name of reduction to use (e.g. \code{"pca"}, \code{"harmony"}). Default \code{"pca"}.
#'
#' @return A factor vector of cluster assignments (\code{cluster1}, \code{cluster2}, ...).
#'
#' @details Stops with an error if \code{reduction_type} is not present in \code{seurat_obj@reductions}.
#'
#' @importFrom Seurat Embeddings FindNeighbors FindClusters Idents
#' @export
cluster_using_leiden <- function(seurat_obj, leiden_res = 0.1, reduction_type = "pca") {
  # Existence check for reduction
  if (!(reduction_type %in% names(seurat_obj@reductions))) {
    stop(paste0("Reduction '", reduction_type, "' not found in the Seurat object."))
  }
  # Determine number of components to use
  num_pcs <- ncol(Embeddings(seurat_obj, reduction = reduction_type))
  seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:num_pcs, reduction = reduction_type)
  seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = leiden_res)
  leiden_assignment_list <- Seurat::Idents(seurat_obj)
  leiden_assignment_list <- factor(paste0("cluster", as.numeric(leiden_assignment_list)))
  levels(leiden_assignment_list) <- paste0("cluster", seq_along(levels(leiden_assignment_list)))
  return(leiden_assignment_list)
}


#' Perform k-means clustering on a Seurat object
#'
#' Runs k-means on the embedding of a specified dimensionality reduction and
#' returns factor labels \code{cluster1}, \code{cluster2}, etc.
#'
#' @param seurat_obj A Seurat object with the requested reduction already computed.
#' @param kmeans_k Number of clusters (default \code{4}).
#' @param reduction_type Name of reduction to use (e.g. \code{"pca"}, \code{"harmony"}). Default \code{"pca"}.
#'
#' @return A factor vector of cluster assignments (\code{cluster1}, \code{cluster2}, ...).
#'
#' @details Stops with an error if \code{reduction_type} is not in \code{seurat_obj@reductions}. Uses \code{mclust::kmeans}.
#'
#' @importFrom Seurat Embeddings
#' @importFrom stats kmeans
#' @export
cluster_using_kmeans <- function(seurat_obj, kmeans_k = 4, reduction_type = "pca") {
  # Existence check for reduction
  if (!(reduction_type %in% names(seurat_obj@reductions))) {
    stop(paste0("Reduction '", reduction_type, "' not found in the Seurat object."))
  }
  embeds <- Seurat::Embeddings(seurat_obj, reduction = reduction_type)
  kmeans_assignment_list <- as.factor(stats::kmeans(embeds, centers = kmeans_k, iter.max = 100)$cluster)
  kmeans_assignment_list <- factor(paste0("cluster", as.numeric(kmeans_assignment_list)))
  levels(kmeans_assignment_list) <- paste0("cluster", seq_along(levels(kmeans_assignment_list)))
  return(kmeans_assignment_list)
}



#' Run PCA, Harmony, and clustering on a SingleCellExperiment
#'
#' Converts the SCE to Seurat, normalizes, runs PCA and Harmony (optional batch
#' correction), then performs Leiden or k-means clustering on both PCA and
#' Harmony embeddings.
#'
#' @param sce A \code{SingleCellExperiment} with \code{counts} and \code{colData}
#'   (e.g. \code{sample}, \code{batch} if \code{include_batch} is \code{TRUE}).
#' @param cluster_type \code{"leiden"} (default) or \code{"kmeans"}.
#' @param leiden_res Resolution for Leiden (default \code{0.1}).
#' @param kmeans_k Number of clusters for k-means (default \code{6}).
#' @param num_pcs Number of principal components (default \code{50}).
#' @param include_batch If \code{TRUE}, Harmony uses \code{sample} and \code{batch}; else only \code{sample}. Default \code{TRUE}.
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{pca_embeds}: PCA embedding matrix.
#'     \item \code{harmony_embeddings}: Harmony embedding matrix.
#'     \item \code{cluster_assignment_list}: Cluster factor from PCA.
#'     \item \code{cluster_assignment_list_harmony}: Cluster factor from Harmony.
#'     \item \code{lblnorm_counts}: Log-normalized counts matrix.
#'   }
#'
#' @importFrom Seurat CreateSeuratObject AddMetaData NormalizeData GetAssayData FindVariableFeatures ScaleData RunPCA Embeddings
#' @importFrom harmony RunHarmony
#' @importFrom SingleCellExperiment counts colData
#' @export
run_pca_harmony_leiden <- function(
  sce,
  cluster_type = "leiden",
  leiden_res = 0.1,
  kmeans_k = 6,
  num_pcs = 50,
  include_batch = TRUE
) {
  # Preprocess the SingleCellExperiment using Seurat
  print("Preprocessing seurat object")
  seurat_obj <- Seurat::CreateSeuratObject(counts = counts(sce))
  seurat_obj <- Seurat::AddMetaData(seurat_obj, metadata = data.frame(colData(sce)))
  
  # Normalize and extract log-library size normalized counts
  seurat_obj <- Seurat::NormalizeData(seurat_obj)
  lblnorm_counts <- Seurat::GetAssayData(seurat_obj, assay = "RNA", layer = "data")
  
  # Find variable genes and scale for PCA
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj)
  seurat_obj <- Seurat::ScaleData(seurat_obj)
  
  # PCA and extraction of embeddings
  print("Running seurat pca and clustering")
  seurat_obj <- Seurat::RunPCA(seurat_obj, npcs = num_pcs)
  pca_embeds <- Seurat::Embeddings(seurat_obj, reduction = "pca")
  
  # Run Harmony for batch correction (optionally including "batch" in the model)
  if (include_batch) {
    seurat_obj <- harmony::RunHarmony(seurat_obj, group.by.vars = c("sample", "batch"))
  } else {
    seurat_obj <- harmony::RunHarmony(seurat_obj, group.by.vars = c("sample"))
  }
  harmony_embeddings <- Seurat::Embeddings(seurat_obj, reduction = "harmony")
  
  # Clustering on PCA embeddings
  cluster_assignment_list <- NULL
  if (cluster_type == "leiden") {
    cluster_assignment_list <- cluster_using_leiden(seurat_obj, leiden_res = leiden_res, reduction_type = "pca")
  } else if (cluster_type == "kmeans") {
    cluster_assignment_list <- cluster_using_kmeans(seurat_obj, kmeans_k = kmeans_k, reduction_type = "pca")
  } else {
    stop("cluster_type must be either 'leiden' or 'kmeans'.")
  }
  
  # Clustering on Harmony embeddings
  cluster_assignment_list_harmony <- NULL
  if (cluster_type == "leiden") {
    cluster_assignment_list_harmony <- cluster_using_leiden(seurat_obj, leiden_res = leiden_res, reduction_type = "harmony")
  } else if (cluster_type == "kmeans") {
    cluster_assignment_list_harmony <- cluster_using_kmeans(seurat_obj, kmeans_k = kmeans_k, reduction_type = "harmony")
  }
  
  # Assemble return list
  to_return <- list(
    pca_embeds                  = pca_embeds,
    harmony_embeddings          = harmony_embeddings,
    cluster_assignment_list     = cluster_assignment_list,
    cluster_assignment_list_harmony = cluster_assignment_list_harmony,
    lblnorm_counts              = lblnorm_counts
  )
  return(to_return)
}








#' Overlap between clustering and ground-truth for mSim simulation
#'
#' For mSim-style simulations, computes how many cells from each ground-truth
#' DE-gene cluster fall into each predicted cluster.
#'
#' @param clustering_assignment Factor or character. Predicted cluster per cell.
#' @param reduced_par List with \code{cd} (e.g. \code{cluster_id}) and \code{gs_idx} (DE gene indices per cluster).
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{cell_count_assign_de_gene}: Matrix  of cell counts with dimensions DE genes x predicted clusters.
#'     \item \code{gt_size_de_gene}: Integer vector of true cluster size per DE gene.
#'     \item \code{assign_size_clust}: Integer vector of cell count per predicted cluster.
#'   }
#'
#' @export
calc_overlap_mSim <- function(clustering_assignment, reduced_par) {
  gt_clusters <- reduced_par$cd$cluster_id
  cont_table <-  table(clustering_assignment, gt_clusters)
  clust_assigns_names <- rownames(cont_table)
  gt_clusts_names <- colnames(cont_table)
  genes_de <- unique(unlist(reduced_par$gs_idx["de", ]))
  cell_count_assign_de_gene <- matrix(
    0, nrow = length(genes_de), ncol = length(clust_assigns_names),
    dimnames = list(genes_de, clust_assigns_names)
  )
  
  # For each predicted cluster, count how many DE gene cells it captures from ground-truth clusters
  for (clust_name in clust_assigns_names) {
    print(paste0("Clust name: ", clust_name))
    for (gt_clust_name in gt_clusts_names) {
      print(paste0("GT clust name: ", gt_clust_name))
      gt_de_genes <- unlist(reduced_par$gs_idx["de", gt_clust_name])
      print(paste0("GT de genes: ", sum((genes_de %in% gt_de_genes))))
      # Increment counts for DE genes present in this ground-truth cluster
      cell_count_assign_de_gene[genes_de %in% gt_de_genes, clust_name] <-
        cell_count_assign_de_gene[genes_de %in% gt_de_genes, clust_name] + cont_table[clust_name, gt_clust_name]
    }
  }

  # Calculate the total size (number of cells) of each DE gene's true ground-truth cluster
  gt_size_de_gene <- integer(length = length(genes_de))
  for (gt_clust_name in gt_clusts_names) {
    size_gt_clust <- sum(cont_table[, gt_clust_name])
    gt_size_de_gene[genes_de %in% unlist(reduced_par$gs_idx["de", gt_clust_name])] <- size_gt_clust
  }

  # Number of cells in each predicted cluster
  assign_size_clust <- integer(length = length(clust_assigns_names))
  for (clust_name in clust_assigns_names) {
    assign_size_clust[clust_name] <- sum(cont_table[clust_name, ])
  }

  de_gene_overlap_info <- list(
    cell_count_assign_de_gene = cell_count_assign_de_gene,
    gt_size_de_gene = gt_size_de_gene,
    assign_size_clust = assign_size_clust
  )
  return(de_gene_overlap_info)
}

#' Overlap between clustering and ground-truth for augData simulation
#'
#' For augmented-data simulations, counts how many cells in each predicted
#' cluster are DE for each simulated DE gene, using \code{de_args} from the
#' augmented SCE.
#'
#' @param clustering_assignment Factor or character. Predicted cluster per cell.
#' @param pa_de_int_df Data frame or list with \code{num_de_genes}.
#' @param de_args Data frame with columns \code{name}, \code{is_simulated}, \code{is_de_cell} (list column).
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{cell_count_assign_de_gene}: Matrix of DE cell counts with dimensions DE genes x clusters.
#'     \item \code{gt_size_de_gene}: Named integer of true DE cell count per gene.
#'     \item \code{assign_size_clust}: Table of cell counts per predicted cluster.
#'   }
#'
#' @export
calc_overlap_augData <- function(clustering_assignment, pa_de_int_df, de_args) {
  clust_assigns_names <- names(table(clustering_assignment))
  nde_genes_current <- pa_de_int_df$num_de_genes

  # Identify simulated DE genes
  de_genes <- de_args[de_args$is_simulated, ]$name
  if (nde_genes_current == 0) {
    de_genes <- character(0)
  } else {
    de_genes <- de_genes[1:nde_genes_current]
  }
  sim_neigh <- de_args[de_args$is_simulated, c("name", "is_de_cell")]
  sim_neigh <- sim_neigh[sim_neigh$name %in% de_genes, ]
  
  # Matrix: rows DE genes, columns clusters; each cell is count of DE cells for that gene/cluster
  cell_count_assign_de_gene <- matrix(
    0, nrow = length(de_genes), ncol = length(clust_assigns_names),
    dimnames = list(de_genes, clust_assigns_names)
  )
  gt_size_de_gene <- integer(length = length(de_genes))
  names(gt_size_de_gene) <- de_genes
  assign_size_clust <- table(clustering_assignment)

  for (gene in de_genes) {
    # The total number of true DE cells for this gene
    gt_size_de_gene[gene] <- sum(sim_neigh[sim_neigh$name == gene, ]$is_de_cell[[1]])
    # For each cluster, count cells that are in predicted cluster and are DE for this gene
    for (clust_name in clust_assigns_names) {
      cell_count_assign_de_gene[gene, clust_name] <-
        sum((clustering_assignment == clust_name) & (sim_neigh[sim_neigh$name == gene, ]$is_de_cell[[1]]))
    }
  }
  de_gene_overlap_info <- list(
    cell_count_assign_de_gene = cell_count_assign_de_gene,
    gt_size_de_gene = gt_size_de_gene,
    assign_size_clust = assign_size_clust
  )
  return(de_gene_overlap_info)
}

#' Overlap between clustering and ground-truth (mSim or augData)
#'
#' Dispatches to \code{\link{calc_overlap_mSim}} or \code{\link{calc_overlap_augData}}
#' based on \code{sim_type}. For mSim, reads \code{reduced_par} from
#' \code{pa_de_int_df$reduced_par_list_fn}; for augData, uses \code{rowData(used_sce)} as \code{de_args}.
#'
#' @param clustering_assignment Factor or character. Predicted cluster per cell.
#' @param used_sce \code{SingleCellExperiment} (used for augData to get \code{rowData}).
#' @param pa_de_int_df List or data frame with \code{num_de_genes}; for mSim also \code{reduced_par_list_fn}.
#' @param sim_type \code{"mSim"} or \code{"augData"}.
#'
#' @return The list returned by \code{calc_overlap_mSim} or \code{calc_overlap_augData}.
#'
#' @importFrom SingleCellExperiment rowData
#' @export
calc_overlap <- function(clustering_assignment, used_sce, pa_de_int_df, sim_type) {
  if (sim_type == "mSim") {
    reduced_par_list_fn <- pa_de_int_df$reduced_par_list_fn
    print(paste0("Reading reduced par list from: ", reduced_par_list_fn))
    reduced_par <- readRDS(reduced_par_list_fn)
    de_gene_overlap_info <- calc_overlap_mSim(clustering_assignment, reduced_par)
  } else if (sim_type == "augData") {
    de_args <- as.data.frame(SingleCellExperiment::rowData(used_sce))
    de_gene_overlap_info <- calc_overlap_augData(clustering_assignment, pa_de_int_df, de_args)
  } else {
    stop("sim_type must be \"mSim\" or \"augData\".")
  }
  de_gene_overlap_info
}

#' Pseudobulk differential expression (muscat) by cluster and sample
#'
#' Copies \code{colData} columns to \code{cluster_id}, \code{sample_id}, \code{group_id},
#' aggregates with \code{muscat::aggregateData}, then runs \code{muscat::pbDS}.
#'
#' @param sce \code{SingleCellExperiment} with \code{counts} and \code{colData} containing cluster, sample, and group columns.
#' @param cluster_id_col Character. \code{colData} column name for cluster.
#' @param sample_id_col Character. \code{colData} column name for sample.
#' @param group_id_col Character. \code{colData} column name for condition/group.
#' @param method \code{"DESeq2"} or \code{"edgeR"} (passed to \code{muscat::pbDS}). Default \code{"DESeq2"}.
#'
#' @return Data frame: DE table from \code{pbDS} (e.g. \code{res_msct$table[[1]]}).
#'
#' @importFrom muscat aggregateData pbDS
#' @importFrom S4Vectors metadata
#' @export
run_mSim_type_DE_analysis <- function(sce, cluster_id_col, sample_id_col, group_id_col, method = "DESeq2") {
  # Assign cluster, sample and group IDs
  sce$cluster_id <- sce[[cluster_id_col]]
  sce$sample_id <- sce[[sample_id_col]]
  sce$group_id <- sce[[group_id_col]]

  # Make sure group_id is a factor with proper levels
  sce$group_id <- factor(sce$group_id)
  sce$cluster_id <- factor(sce$cluster_id)
  sce$cell_type <- factor(sce$cell_type)

  # aggregate into pseudobulk
  pseudobulk <- muscat::aggregateData(  
    sce,
    assay = "counts",
    fun = "sum",
    by = c("cluster_id", "sample_id")
  ) 

  # Create experiment_info data frame
  sample_info <- unique(data.frame(
    sample_id = sce$sample_id,
    group_id = sce$group_id,
    stringsAsFactors = TRUE  # Ensure these are factors
  ))
  S4Vectors::metadata(pseudobulk)$experiment_info <- sample_info

  # Perform differential expression analysis
  res_msct <- muscat::pbDS(
    pseudobulk, 
    method = method, 
    verbose = TRUE 
  )
  return(res_msct$table[[1]])
}

#' Run full analysis for one clustering: overlap, DE, optional PVE
#'
#' Saves cluster assignment and overlap info to RDS, runs pseudobulk DE via
#' \code{\link{run_mSim_type_DE_analysis}}, adjusts p-values with
#' \code{\link{adjust_all_pvals}}, and optionally runs and saves PVE with
#' \code{\link{get_PVE_of_genes}}.
#'
#' @param clustering_assignment_curr_anls Factor. Cluster assignment per cell.
#' @param lblnorm_counts Matrix. Log-normalized counts (for PVE when \code{run_PVE} is \code{TRUE}).
#' @param used_sce \code{SingleCellExperiment}. Used for overlap and DE.
#' @param pa_de List with \code{num_de_genes} and (for mSim) \code{reduced_par_list_fn}; passed to \code{\link{calc_overlap}}.
#' @param sim_type \code{"mSim"} or \code{"augData"}.
#' @param curr_anls Character. Analysis label (e.g. \code{harmony_}, \code{noharm_}).
#' @param new_id_check Character. Run ID used in file names.
#' @param anls_info_dir Character. Directory for intermediate analysis RDS files.
#' @param des_dir Character. Directory for final DE result RDS files.
#' @param clustering_prefix Character. Prefix for cluster assignment RDS filenames.
#' @param de_overlap_info_save_prfx Character. Prefix for overlap info RDS filenames.
#' @param res_de_save_prfx Character. Prefix for DE result RDS filenames.
#' @param PVE_metrics_save_prfx Character. Prefix for PVE metrics RDS filenames.
#' @param anls_patterns Named list. Maps analysis label (e.g. \code{curr_anls}) to filename infix.
#' @param cut_off_true,cut_off_false Numeric. Reserved for downstream use; not used here.
#' @param sig_threshold Numeric. Threshold for adjusted p-value to define DE genes.
#' @param overlap_type Unused; kept for interface compatibility.
#' @param run_PVE If \code{TRUE}, run and save PVE metrics (default \code{FALSE}).
#'
#' @return A list with \code{res_DE} (DE table) and \code{de_overlap_info} (overlap list).
#'
#' @importFrom SingleCellExperiment rowData colData
#' @importFrom SummarizedExperiment assay
#' @export
run_analysis_for_clustering <- function(clustering_assignment_curr_anls, lblnorm_counts, used_sce, pa_de, sim_type,
                                      curr_anls, new_id_check,
                                      anls_info_dir, des_dir, clustering_prefix, de_overlap_info_save_prfx,
                                      res_de_save_prfx, PVE_metrics_save_prfx, anls_patterns,
                                      cut_off_true, cut_off_false, sig_threshold, overlap_type, run_PVE = FALSE) {

  print(paste0("Running analysis for clustering: ", curr_anls))
  saveRDS(clustering_assignment_curr_anls, file.path(anls_info_dir, paste0(clustering_prefix, anls_patterns[[curr_anls]], new_id_check, ".rds")))
  print("Calculating de overlap info")
  de_overlap_info_curr_anls <- calc_overlap(clustering_assignment_curr_anls, used_sce, pa_de, sim_type)
  saveRDS(de_overlap_info_curr_anls, file.path(anls_info_dir, paste0(de_overlap_info_save_prfx, anls_patterns[[curr_anls]], new_id_check, ".rds")))

  print("Running de analysis")
  used_sce$data_cluster <- clustering_assignment_curr_anls
  res_de_curr_anls <- run_mSim_type_DE_analysis(used_sce, "data_cluster", "sample", "fake_condition", method = "edgeR")
  saveRDS(res_de_curr_anls, file.path(des_dir, paste0(res_de_save_prfx, anls_patterns[[curr_anls]], new_id_check, ".rds")))
  combined_sim_res_DE_tested <- adjust_all_pvals(res_de_curr_anls)
  gene_clusts_DE <- combined_sim_res_DE_tested[combined_sim_res_DE_tested$adj_p_val < sig_threshold, ]
  genes_DE <- unique(combined_sim_res_DE_tested$gene)[1:100]

  if (run_PVE) {
    used_sce$cluster_id <- clustering_assignment_curr_anls
    start_time <- Sys.time()
    pve_results <- get_PVE_of_genes(used_sce, lblnorm_counts, use_pseudobulk = TRUE)
    end_time <- Sys.time()
    print(paste0("Time taken to run PVE analysis: ", end_time - start_time))

    saveRDS(pve_results, file.path(anls_info_dir, paste0(PVE_metrics_save_prfx, anls_patterns[[curr_anls]], new_id_check, ".rds")))
  }
  
   ###### first attempt at PCA_Variance_Explained analysis
  # cell_embeddings <- clust_results[[embedding_name]]
  # cluster_id_cells <- cluster_assignment_curr_anls
  # fake_condition_cells <- used_sce$fake_condition
  
  # # Compute R^2 for each principal component (PC) regressed on cluster and fake condition
  # get_r2_per_pc <- function(embeddings, predictor) {
  #   # predictor: factor to regress against
  #   apply(embeddings, 2, function(pc_vec) {
  #     mod_full <- summary(lm(pc_vec ~ predictor))
  #     mod_full$r.squared
  #   })
  # }
  # r2_by_cluster <- get_r2_per_pc(cell_embeddings, cluster_id_cells)
  # r2_by_fake_condition <- get_r2_per_pc(cell_embeddings, fake_condition_cells)

  # genes_DE_lbl_counts <- lblnorm_counts[genes_DE, ]

  # pca_loadings <- clust_results[[loadings_name]]



  
  # # Raw variance of each gene explained by each PC via multivariable lm; partial SS from drop1.
  # pcs <- cell_embeddings[, colnames(pca_loadings), drop = FALSE]
  # gene_pc_r2_matrix <- matrix(NA_real_, nrow = length(genes_DE), ncol = ncol(pca_loadings),
  #                             dimnames = list(genes_DE, colnames(pca_loadings)))
  # for (gene in genes_DE) {
  #   print(paste0("Running analysis for gene: ", gene))
  #   y <- as.numeric(genes_DE_lbl_counts[gene, ])
  #   ok <- stats::complete.cases(y, pcs)
  #   if (sum(ok) < 2L) next
  #   fit <- stats::lm(y ~ ., data = data.frame(y = y[ok], pcs[ok, , drop = FALSE]))
  #   drop_tbl <- stats::drop1(fit, scope = ~ ., test = "F")
  #   n_ok <- sum(ok)
  #   var_explained <- drop_tbl[-1L, "Sum of Sq"] / (n_ok - 1L)
  #   gene_pc_r2_matrix[gene, names(var_explained)] <- var_explained
  # }

  
  
  return(list(res_DE = res_de_curr_anls, de_overlap_info = de_overlap_info_curr_anls))
}

#' Proportion of variance explained (PVE) by cluster and condition
#'
#' Fits a variance-partition model \code{~ (1|cluster_id) + (1|sample_id) + (1|fake_condition)}
#' via \code{variancePartition::fitExtractVarPartModel}. Optionally aggregates to
#' pseudobulk first with \code{scuttle::aggregateAcrossCells} and \code{logNormCounts}.
#' \code{colData(sce)} must include \code{cluster_id}, \code{sample_id}, \code{fake_condition}.
#'
#' @param sce \code{SingleCellExperiment} with \code{colData} containing \code{cluster_id}, \code{sample_id}, \code{fake_condition}.
#' @param log_counts_cells Matrix of log-normalized counts. If \code{use_pseudobulk} is \code{TRUE}, can be a placeholder; log-counts are derived from pseudobulk.
#' @param use_pseudobulk If \code{TRUE}, aggregate by cluster and sample then fit (default); else fit on cell-level \code{log_counts_cells}.
#'
#' @return A list with \code{tPVE} (type: cluster) and \code{sPVE} (state: fake_condition), each a named numeric vector (one value per gene).
#'
#' @importFrom scuttle aggregateAcrossCells logNormCounts
#' @importFrom variancePartition fitExtractVarPartModel
#' @importFrom SingleCellExperiment colData
#' @importFrom SummarizedExperiment assay
#' @importFrom S4Vectors DataFrame
#' @importFrom stats setNames
#' @export
get_PVE_of_genes <- function(sce, log_counts_cells, use_pseudobulk = TRUE) {  
  if (use_pseudobulk) {    # Reduces memory usage and computation time for large datasets
    pseudobulk_agg <- scuttle::aggregateAcrossCells(
      sce, 
      ids = DataFrame(cluster_id = SingleCellExperiment::colData(sce)$cluster_id, sample_id = SingleCellExperiment::colData(sce)$sample_id),
      use.assay.type = "counts",
      statistics = "sum"
    )
    pseudobulk_ln <- scuttle::logNormCounts(pseudobulk_agg, assay.type = "counts",size.factors = NULL)
    log_counts_cells <- assay(pseudobulk_ln, "logcounts")
    cd <- data.frame(colData(pseudobulk_ln))
  } else {
    # Cell-level approach (original)
    cd <- data.frame(SingleCellExperiment::colData(sce))
  }
  # Fit variance partition model
  cd <- cd[c("cluster_id", "sample_id", "fake_condition")]
  mod <- ~(1|cluster_id)+(1|sample_id)+(1|fake_condition)
  res <- variancePartition::fitExtractVarPartModel(log_counts_cells, mod, cd, quiet = TRUE)
  
  to_return <- lapply(list(
    tPVE = res$cluster_id, 
    sPVE = res$fake_condition
  ), setNames, rownames(res))
  names(to_return) <- c("tPVE", "sPVE")
  return(to_return)
}

#' Type (cluster) vs state metric per gene via limma
#'
#' Fits \code{~ sample_id + cluster_id} to \code{assay(sce, "logcounts")} and
#' returns \code{1 - adj.P.Val} for the cluster coefficient (type effect).
#' Uses \code{pbapply} for parallelization.
#'
#' @param sce \code{SingleCellExperiment} with \code{logcounts} assay and \code{colData} containing \code{sample_id}, \code{cluster_id}.
#' @param log_counts_cells Unused; log-counts are taken from \code{sce}.
#' @param use_pseudobulk Unused; kept for interface compatibility.
#'
#' @return Named numeric vector: per-gene \code{1 - adj.P.Val} for the cluster effect.
#'
#' @importFrom limma lmFit eBayes topTable
#' @importFrom pbapply pbapply
#' @importFrom SingleCellExperiment colData
#' @importFrom SummarizedExperiment assay
#' @importFrom stats model.matrix
#' @export
get_tF_genes <- function(sce, log_counts_cells, use_pseudobulk = TRUE) {
    # Cell-level approach (original)
  log_counts_cells <- assay(sce, "logcounts")
  cd <- data.frame(colData(sce))
  f <- ~ sample_id + cluster_id
  mm <- stats::model.matrix(f, data=cd)
  rownames(mm) <- colnames(sce)
  # Progress bar version: use pbapply::pbapply for tqdm-like progress

  res <- pbapply::pbapply(log_counts_cells, 1, function(g) {
        fit <- limma::lmFit(g, mm)
        fit <- limma::eBayes(fit, trend=FALSE)
        cs <- colnames(fit$cov.coefficients)
        nan <- !colnames(mm) %in% cs
        if (any(nan)) {
            fit <- limma::lmFit(g, mm[, !nan])
            fit <- limma::eBayes(fit, trend=FALSE)
        }
        cs <- grep("cluster", cs)
        tt <- limma::topTable(fit, coef=cs, sort.by="none")
        as.numeric(1-tt["adj.P.Val"])
  }, cl = 10)
  return(res)
}

#' Per-cluster pseudobulk DE (edgeR) and minimum adjusted p-value per gene
#'
#' Aggregates the SCE by \code{sample_id} and \code{cluster_id} with
#' \code{scuttle::aggregateAcrossCells}, then runs edgeR \code{glmQLFTest} per
#' cluster. Returns \code{-log10(min p_adj)} per gene across clusters.
#'
#' @param sce \code{SingleCellExperiment} with \code{colData} containing \code{sample_id}, \code{cluster_id}, \code{fake_condition}.
#' @param use_pseudobulk Unused; aggregation is always performed. Kept for interface compatibility.
#'
#' @return Named numeric vector of per-gene \code{-log10(min adjusted p-value)} (or a data frame in edge cases).
#'
#' @importFrom scuttle aggregateAcrossCells
#' @importFrom edgeR DGEList calcNormFactors estimateDisp glmQLFit glmQLFTest topTags
#' @importFrom pbapply pblapply
#' @importFrom SingleCellExperiment colData
#' @importFrom SummarizedExperiment assay
#' @export
get_sPBD_genes <- function(sce, use_pseudobulk = TRUE) {
  ids <- SingleCellExperiment::colData(sce)[c("sample_id", "cluster_id")]
  y <- scuttle::aggregateAcrossCells(sce, ids)
  idx <- split(seq_len(ncol(y)), y$cluster_id)
  res <- pbapply::pblapply(names(idx), function(k) {
        z <- y[, idx[[k]]]
        gs <- unique(z$fake_condition)
        if (length(gs) == 2) {
            mm <- stats::model.matrix(~z$fake_condition)
            z <- edgeR::DGEList(assay(z))
            z <- edgeR::calcNormFactors(z)
            z <- edgeR::estimateDisp(z, mm)
            fit <- edgeR::glmQLFit(z, mm)
            lrt <- edgeR::glmQLFTest(fit, 2)
            tbl <- edgeR::topTags(lrt, n=Inf, sort.by="none")$table
            old <- c("logFC", "PValue", "FDR")
            new <- c("lfc", "p_val", "p_adj")
            names(tbl)[match(old, names(tbl))] <- new
            data.frame(
                row.names=NULL,
                gene=rownames(tbl), 
                cluster_id=k, tbl)
        }
    })
    res <- do.call(rbind, res)
    if (!is.null(res)) {
        ps <- split(res$p_adj, res$gene)
        ps <- vapply(ps, min, numeric(1))
        to_return <- -log(ps)[rownames(sce)]
        return(to_return)
    }
    setNames(res, rownames(sce))
  return(res)
}

#' Find top type-vs-state (TvsS) genes for clustering
#'
#' Computes PVE (type and state) via \code{\link{get_PVE_of_genes}}, ranks genes
#' by \code{rank(tPVE) - rank(sPVE)}, and returns the top
#' \code{num_genes_for_TvsS_clustering} genes. Saves PVE metrics to \code{anls_info_dir}.
#'
#' @param used_sce \code{SingleCellExperiment}. \code{cluster_id} is set from \code{clust_results$cluster_assignment_list} inside the function.
#' @param clust_results List with \code{cluster_assignment_list} and \code{lblnorm_counts}.
#' @param num_genes_for_TvsS_clustering Number of top type-over-state genes to return.
#' @param new_id_check Run ID used when saving PVE RDS.
#' @param anls_info_dir Character. Directory for intermediate analysis RDS files.
#' @param PVE_metrics_save_prfx Character. Prefix for PVE metrics RDS filenames.
#' @param use_pseudobulk_TvsS Passed to \code{\link{get_PVE_of_genes}} (default \code{TRUE}).
#'
#' @return Character vector of gene names (top type-over-state genes).
#'
#' @importFrom SingleCellExperiment colData
#' @export
find_TvsS_genes <- function(used_sce, clust_results, num_genes_for_TvsS_clustering, new_id_check,
    anls_info_dir, PVE_metrics_save_prfx, use_pseudobulk_TvsS = TRUE) {
  used_sce$cluster_id <- clust_results$cluster_assignment_list
  log_counts_cells <- clust_results$lblnorm_counts
  PVE_metrics <- get_PVE_of_genes(used_sce, log_counts_cells, use_pseudobulk = use_pseudobulk_TvsS)
  saveRDS(PVE_metrics, file.path(anls_info_dir, paste0(PVE_metrics_save_prfx, new_id_check, ".rds")))
  # sPBD_metrics <- get_sPBD_genes(used_sce, use_pseudobulk = use_pseudobulk_TvsS)
  type_metric <- PVE_metrics$tPVE
  state_metric <- PVE_metrics$sPVE
  rank_type <- rank(type_metric)
  rank_state <- rank(state_metric)
  rank_diff <- rank_type - rank_state
  top_n_genes <- names(rank_diff[order(rank_diff, decreasing = TRUE)][1:num_genes_for_TvsS_clustering])
  return(top_n_genes)
}



#' Adjust all p-values across clusters (Benjamini-Hochberg)
#'
#' Pools gene-cluster p-values from all clusters and applies a single
#' Benjamini-Hochberg correction. Used internally by
#' \code{\link{run_analysis_for_clustering}}.
#'
#' @param de_pvals_by_cluster A named list; each element is a data frame with
#'   columns \code{gene}, \code{cluster_id}, and \code{p_val}.
#'
#' @return A data frame with columns \code{gene}, \code{cluster_id}, \code{p_val},
#'   and \code{adj_p_val} (BH-adjusted).
#'
#' @importFrom stats p.adjust
#' @keywords internal
adjust_all_pvals <- function(de_pvals_by_cluster) {
  # Initialize an empty data.frame to store the combined results
  sim_res_DE <- data.frame(
    gene = character(),
    cluster_id = character(),
    p_val = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop through each cluster and combine results
  for (clust_name in names(de_pvals_by_cluster)) {
    curr_clust_res_DE <- de_pvals_by_cluster[[clust_name]]
    # Filter out NA p-values
    curr_clust_res_DE <- curr_clust_res_DE[!is.na(curr_clust_res_DE$p_val), ]
    # Keep only relevant columns
    curr_clust_res_DE <- curr_clust_res_DE[, c("gene", "cluster_id", "p_val")]
    # Append to aggregated results
    sim_res_DE <- rbind(sim_res_DE, curr_clust_res_DE)
  }
  
  # Adjust p-values across all genes using Benjamini-Hochberg method
  sim_res_DE$adj_p_val <- stats::p.adjust(sim_res_DE$p_val, method = "BH")
  
  return(sim_res_DE)
}