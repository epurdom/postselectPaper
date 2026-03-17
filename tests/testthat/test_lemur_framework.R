# Unit tests for prep_aug_sim and create_augmented_data (lemur_framework.R)
# Uses synthetic data: minimal Seurat with 100 cells, 100 genes; generates 30 DE genes
# and checks structure, cluster assignment, and that count differences match lfc sign.

test_that("prep_aug_sim and create_augmented_data run with synthetic data and respect simulation contract", {
  set.seed(42)
  # Suppress expected Seurat/irlba warnings for small synthetic data
  suppressWarnings({

  # Synthetic data: 100 genes x 100 cells
  n_genes <- 100L
  n_cells <- 100L
  counts <- matrix(
    as.integer(exp(rnorm(n_genes * n_cells, mean = 2, sd = 1.5))),
    nrow = n_genes,
    ncol = n_cells
  )
  # Avoid zeros for PCA
  counts[counts < 1L] <- 1L
  rownames(counts) <- paste0("gene", seq_len(n_genes))
  colnames(counts) <- paste0("cell", seq_len(n_cells))

  # Sample and batch for colData (required for prep_aug_sim / Harmony; we skip Harmony)
  n_samples <- 10L
  sample_id <- paste0("s", rep(seq_len(n_samples), length.out = n_cells))
  meta <- data.frame(
    sample = sample_id,
    group = rep("A", n_cells),
    batch = paste0("b", (seq_len(n_cells) - 1L) %% 3L + 1L),
    row.names = colnames(counts)
  )

  # Build minimal Seurat with PCA (required when passing processed_sim)
  seurat_obj <- Seurat::CreateSeuratObject(counts = counts, meta.data = meta)
  seurat_obj <- Seurat::NormalizeData(seurat_obj)
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, nfeatures = min(50L, n_genes))
  seurat_obj <- Seurat::ScaleData(seurat_obj)
  npcs <- min(50L, n_genes - 1L, n_cells - 1L)
  seurat_obj <- Seurat::RunPCA(seurat_obj, npcs = npcs, verbose = FALSE)

  # prep_aug_sim uses processed_sim@assays$RNA@layers$counts (Seurat v5); skip if not available
  has_layers_counts <- tryCatch({
    dim(seurat_obj@assays$RNA@layers$counts)
    TRUE
  }, error = function(e) FALSE)
  skip_if(!has_layers_counts, "Seurat assay has no layers$counts (requires Seurat v5); skipping test")

  # prep_aug_sim args (no config; explicit main_covariate, sample_covariate, etc.)
  args_prep <- list(
    seed = 42L,
    randomization = "samples",
    clustering = "kmeans",
    cut_at = c(3, 10),
    lfc_mean = 1.5,
    n_hvgs = n_genes,
    datapath = NULL,
    dest_file = NULL,
    condition = NULL,
    use_harmony = FALSE,
    leiden_resolution = NULL,
    leiden_clusters = NULL,
    main_covariate = "group",
    sample_covariate = "sample",
    assay_continuous = "logcounts",
    cell_type_column = "cell_type",
    batch_covariate = "batch",
    processed_sim = seurat_obj,
    num_pcs = npcs,
    verbose = FALSE
  )
  pa <- do.call(postselect::prep_aug_sim, args_prep)

  expect_equal(pa$nc, n_cells)
  expect_gte(nrow(pa$pca_embeds), 1L)
  expect_lte(nrow(pa$pca_embeds), npcs)
  expect_equal(ncol(pa$pca_embeds), n_cells)
  expect_equal(length(pa$sample), n_cells)
  expect_equal(length(pa$sf), n_cells)

  # create_augmented_data: 30 DE genes
  n_de_genes <- 30L
  pa$n_de_genes <- n_de_genes
  args_aug <- pa[names(pa) %in% formalArgs(postselect::create_augmented_data)]
  args_aug$verbose <- FALSE
  augmented_data <- do.call(postselect::create_augmented_data, args_aug)

  # --- Structure ---
  expect_equal(dim(augmented_data$new_counts), c(n_de_genes, n_cells))
  expect_equal(nrow(augmented_data$de_args), n_de_genes)
  expect_true(all(c("name", "cut_at", "lfc", "sel_cluster", "is_de_cell") %in% names(augmented_data$de_args)))
  expect_equal(length(augmented_data$fake_condition), n_cells)
  expect_true(all(augmented_data$fake_condition %in% c("fake_ctrl", "fake_trt")))
  expect_type(augmented_data$kmeans_clusterings, "list")
  expect_true(length(augmented_data$kmeans_clusterings) >= 1L)

  de_args <- augmented_data$de_args
  kmeans_clusterings <- augmented_data$kmeans_clusterings
  new_counts <- augmented_data$new_counts
  fake_condition <- augmented_data$fake_condition

  # --- Cluster assignment: sel_cluster in valid range and is_de_cell matches clustering ---
  for (i in seq_len(n_de_genes)) {
    k <- de_args$cut_at[i]
    expect_true(k %in% c(3, 10), info = paste("gene", i, "cut_at"))
    clust_vec <- kmeans_clusterings[[as.character(k)]]
    expect_length(clust_vec, n_cells)
    expect_true(de_args$sel_cluster[i] %in% clust_vec, info = paste("gene", i, "sel_cluster"))
    expected_de_cells <- as.vector(clust_vec == de_args$sel_cluster[i])
    expect_equal(as.vector(as.logical(de_args$is_de_cell[[i]])), expected_de_cells, info = paste("gene", i, "is_de_cell"))
  }

  # --- Distribution: mean count in DE cells (trt vs ctrl) has same sign as lfc ---
  n_correct_sign <- 0L
  for (i in seq_len(n_de_genes)) {
    is_de <- as.logical(de_args$is_de_cell[[i]])
    if (!any(is_de)) next
    de_trt <- is_de & (fake_condition == "fake_trt")
    de_ctrl <- is_de & (fake_condition == "fake_ctrl")
    if (!any(de_trt) || !any(de_ctrl)) next
    mean_trt <- mean(new_counts[i, de_trt])
    mean_ctrl <- mean(new_counts[i, de_ctrl])
    diff_positive <- mean_trt > mean_ctrl
    if (de_args$lfc[i] > 0 && diff_positive) n_correct_sign <- n_correct_sign + 1L
    if (de_args$lfc[i] < 0 && !diff_positive) n_correct_sign <- n_correct_sign + 1L
  }
  # With fixed seed we expect most genes to show the right direction (allow a few failures from randomness)
  expect_true(n_correct_sign > n_de_genes * 0.5,
    info = paste("most genes should have mean count difference matching lfc sign; got", n_correct_sign, "of", n_de_genes))

  }) # end suppressWarnings
})
