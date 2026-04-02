#' Create simulated data based on simulation type
#'
#' Dispatches to the appropriate simulation pipeline: \code{"augData"} uses
#' \code{\link{create_full_augmented_data_w_num_sim_genes}}; \code{"mSim"} is
#' not implemented and will raise an error.
#'
#' @param sim_type Character. One of \code{"mSim"} or \code{"augData"}.
#' @param sim_prefix Character. Prefix for simulation file naming.
#' @param new_id_check Character. Run ID or tag for the current simulation.
#' @param num_de_genes Integer. Number of differentially expressed genes to simulate.
#' @param base_seurat_fn Character. Full path to base Seurat RDS file.
#' @param sims_info_dir Character. Directory for simulation parameter and output RDS files.
#' @param param_file_prefix Character. Prefix for parameter list RDS (e.g. \code{pa_preseurat_}).
#' @param de_args_save_prfx Character. Prefix for saving \code{de_args} RDS.
#' @param condition_save_prfx Character. Prefix for saving \code{fake_condition} RDS.
#' @param sample_ids_save_prfx Character. Prefix for saving sample IDs RDS.
#' @param gt_clusterings_save_prfx Character. Prefix for saving ground-truth clusterings RDS.
#' @param verbose Logical. If \code{TRUE}, print progress (default \code{TRUE}).
#' @param save_down Logical. If \code{TRUE}, save intermediate RDS files (default \code{TRUE}).
#'
#' @return A list as returned by the underlying simulation function (e.g. \code{used_sce} and \code{pa_de} for augData).
#'
#' @export
create_simulated_data <- function(sim_type, sim_prefix, new_id_check, num_de_genes,
                                  base_seurat_fn, sims_info_dir, param_file_prefix,
                                  de_args_save_prfx, condition_save_prfx, sample_ids_save_prfx, gt_clusterings_save_prfx,
                                  verbose = TRUE, save_down = TRUE) {
  if (sim_type == "mSim") {
    print("Creating simulated data for mSim")
    stop("Not implemented yet")
  } else if (sim_type == "augData") {
    sim_data <- create_full_augmented_data_w_num_sim_genes(
      sim_prefix = sim_prefix, new_id_check = new_id_check, num_sim_genes = num_de_genes,
      base_seurat_fn = base_seurat_fn, sims_info_dir = sims_info_dir, param_file_prefix = param_file_prefix,
      de_args_save_prfx = de_args_save_prfx, condition_save_prfx = condition_save_prfx,
      sample_ids_save_prfx = sample_ids_save_prfx, gt_clusterings_save_prfx = gt_clusterings_save_prfx,
      verbose = verbose, save_down = save_down)
  } else {
    stop("Invalid sim_type")
  }
  return(sim_data)
}

#' Create augmented single-cell data with a fixed number of simulated DE genes
#'
#' Loads the base Seurat object and parameter list, runs \code{\link{prep_aug_sim}}
#' and \code{\link{create_augmented_data}} (from this package), then builds a
#' \code{SingleCellExperiment} combining real and simulated counts. Optionally
#' saves \code{de_args}, \code{fake_condition}, sample IDs, and ground-truth
#' clusterings under \code{sims_info_dir}.
#'
#' @param sim_prefix Character. Simulation prefix (e.g. \code{augSimLeiden3}).
#' @param new_id_check Character. Run ID (e.g. \code{sim_prefix_id_nsim_500}).
#' @param num_sim_genes Integer. Number of simulated DE genes.
#' @param base_seurat_fn Character. Full path to base Seurat RDS file.
#' @param sims_info_dir Character. Directory for simulation parameter and output RDS files.
#' @param param_file_prefix Character. Prefix for parameter list RDS (e.g. \code{pa_preseurat_}).
#' @param de_args_save_prfx Character. Prefix for saving \code{de_args} RDS.
#' @param condition_save_prfx Character. Prefix for saving \code{fake_condition} RDS.
#' @param sample_ids_save_prfx Character. Prefix for saving sample IDs RDS.
#' @param gt_clusterings_save_prfx Character. Prefix for saving ground-truth clusterings RDS.
#' @param verbose Logical. If \code{TRUE}, print progress (default \code{TRUE}).
#' @param save_down Logical. If \code{TRUE}, save intermediate RDS files (default \code{TRUE}).
#'
#' @return A list with:
#'   \itemize{
#'     \item \code{used_sce}: \code{SingleCellExperiment} with real + simulated counts and \code{fake_condition} in \code{colData}.
#'     \item \code{pa_de}: list with \code{new_id_check}, \code{set_de_genes}, \code{leiden_clusters}, file paths, etc.
#'   }
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment colData rowData
#' @importFrom tibble tibble
#' @importFrom dplyr bind_rows
#' @export
create_full_augmented_data_w_num_sim_genes <- function(sim_prefix, new_id_check, num_sim_genes,
    base_seurat_fn, sims_info_dir, param_file_prefix,
    de_args_save_prfx, condition_save_prfx, sample_ids_save_prfx, gt_clusterings_save_prfx,
    verbose = TRUE, save_down = TRUE) {
  print("Loading seurat null genes with leiden clusters")
  print(paste0("Base seurat file: ", base_seurat_fn))
  base_seurat <- readRDS(base_seurat_fn)
  neg_de_args <- tibble::tibble(name = rownames(base_seurat), is_simulated = FALSE)
  neg_de_counts <- base_seurat@assays$RNA@layers$counts

  print(paste0("Sim info file: ", sims_info_dir, "/", param_file_prefix, sim_prefix, ".Rds"))
  pa <- readRDS(paste0(sims_info_dir, "/", param_file_prefix, sim_prefix, ".Rds"))
  pa$n_de_genes <- num_sim_genes

  print("Preparing augmented data - randomizing sample/cell selections")
  start_time <- Sys.time()
  args_prep <- pa[names(pa) %in% formalArgs(prep_aug_sim)]
  args_prep$processed_sim <- base_seurat
  args_prep$verbose <- verbose
  args_out <- do.call(prep_aug_sim, args_prep)
  # add any parameter in pa that is not in args_out to args_out
  for (param in names(pa)) {
    if (!param %in% names(args_out)) {
      args_out[[param]] <- pa[[param]]
    }
  }
  end_time <- Sys.time()
  print(paste0("Time taken to prepare augmented data: ", end_time - start_time))

  print(paste0("Generating augmented data by generating new counts"))
  start_time <- Sys.time()
  args_aug <- args_out[names(args_out) %in% formalArgs(create_augmented_data)]
  args_aug$verbose <- verbose
  augmented_data <- do.call(create_augmented_data, args_aug)
  end_time <- Sys.time()
  counts_aug_data <- augmented_data$new_counts
  gt_clusterings <- augmented_data$clustering_bank
  de_args <- augmented_data$de_args
  fake_condition <- augmented_data$fake_condition
  print(paste0("Time taken to generate augmented data: ", end_time - start_time))

  print("Creating used sce object by gluing together the negative and positive de counts")
  used_sce <- SingleCellExperiment::SingleCellExperiment(
    list(counts = rbind(unname(neg_de_counts), counts_aug_data)),
    colData = base_seurat@meta.data,
    rowData = dplyr::bind_rows(neg_de_args, de_args))
  rownames(used_sce) <- SummarizedExperiment::rowData(used_sce)$name
  SummarizedExperiment::colData(used_sce)$fake_condition <- fake_condition
  set_de_genes <- SummarizedExperiment::rowData(used_sce)[SummarizedExperiment::rowData(used_sce)$is_simulated, "name"]

  print(paste0("Saving augmented data details for simulation ", new_id_check))
  start_time <- Sys.time()

  if (save_down) {
    de_args$is_de_cell <- NULL
    de_args_fn <- paste0(sims_info_dir, "/", de_args_save_prfx, new_id_check, ".Rds")
    print(paste0("Saving de_args for simulation ", new_id_check, " to ", de_args_fn))
    saveRDS(de_args, de_args_fn)

    fake_condition_fn <- paste0(sims_info_dir, "/", condition_save_prfx, new_id_check, ".Rds")
    print(paste0("Saving fake_condition for simulation ", new_id_check, " to ", fake_condition_fn))
    saveRDS(fake_condition, fake_condition_fn)

    sample_ids_fn <- paste0(sims_info_dir, "/", sample_ids_save_prfx, new_id_check, ".Rds")
    print(paste0("Saving sample ids for simulation ", new_id_check, " to ", sample_ids_fn))
    saveRDS(used_sce$sample, sample_ids_fn)

    gt_clusterings_fn <- paste0(sims_info_dir, "/", gt_clusterings_save_prfx, new_id_check, ".Rds")
    print(paste0("Saving gt clusterings for simulation ", new_id_check, " to ", gt_clusterings_fn))
    saveRDS(gt_clusterings, gt_clusterings_fn)
  } else {
    de_args_fn <- NULL
    fake_condition_fn <- NULL
    sample_ids_fn <- NULL
    gt_clusterings_fn <- NULL
  }
  end_time <- Sys.time()
  print(paste0("Time taken to save augmented data details: ", end_time - start_time))

  pa_de <- list(
    new_id_check = new_id_check,
    num_de_genes = num_sim_genes,
    sim_prefix = sim_prefix,
    include_batch = TRUE,
    de_args_fn = de_args_fn,
    fake_condition_fn = fake_condition_fn,
    sample_ids_fn = sample_ids_fn,
    set_de_genes = set_de_genes,
    leiden_clusters = pa$leiden_clusters
  )

  return(list(used_sce = used_sce, pa_de = pa_de))
}

#' Build parameter list for muscat-style full simulation (mSim)
#'
#' mSim is not implemented in this package; this function is a placeholder.
#' It would have constructed a parameter list for external mSim-style simulation
#' (e.g. cluster count \code{nk}, DE proportion \code{p_de_ratio}, \code{lfc},
#' etc.). Calling it stops with an error. This package does not depend on
#' \code{muscatsc} or \code{muscat} for simulation.
#'
#' @param lfc Numeric. Log fold change for DE genes.
#' @param num_de_genes Integer. Number of DE genes.
#'
#' @return Nothing; stops with an error.
#'
#' @export
gen_muscat_params <- function(lfc, num_de_genes) {
  stop("mSim simulation is not implemented; gen_muscat_params is only used by create_full_mSim_data.")
}

#' Generate full mSim (muscat-style) synthetic data
#'
#' mSim simulation is not implemented in this package. This function exists so
#' that \code{\link{create_simulated_data}} can dispatch on \code{sim_type};
#' calling it directly or with \code{sim_type = "mSim"} will throw an error.
#' This package does not depend on \code{muscatsc}.
#'
#' @param sim_prefix Character. Simulation prefix.
#' @param id_check Integer or character. Run ID.
#' @param lfc Numeric. Log fold change for DE.
#' @param num_de_genes Integer. Number of DE genes.
#' @param base_seurat_fn Unused; mSim not implemented.
#' @param sims_info_dir Unused; mSim not implemented.
#' @param param_file_prefix Unused; mSim not implemented.
#' @param de_args_save_prfx Unused; mSim not implemented.
#' @param condition_save_prfx Unused; mSim not implemented.
#' @param sample_ids_save_prfx Unused; mSim not implemented.
#' @param gt_clusterings_save_prfx Unused; mSim not implemented.
#'
#' @return Nothing; stops with an error.
#'
#' @export
create_full_mSim_data <- function(sim_prefix, id_check, lfc, num_de_genes,
    base_seurat_fn = NULL, sims_info_dir = NULL, param_file_prefix = NULL,
    de_args_save_prfx = NULL, condition_save_prfx = NULL, sample_ids_save_prfx = NULL, gt_clusterings_save_prfx = NULL) {
  stop("mSim simulation is not implemented in this package; use sim_type = \"augData\" instead.")
}