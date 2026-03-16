test_that("Overlap Calculations Run",{
  sim_res_DE <- data.frame(
   gene = c("GENE1", "GENE2", "GENE3"),
   cluster_id = c("C1", "C1", "C2"),
   adj_p_val = c(0.01, 0.15, 0.05)
  )
  overlap_matrix <- list(
   cell_count_assign_de_gene = matrix(c(50, 10, 5, 45), nrow = 2,
                                      dimnames = list(c("GENE1", "GENE2"), c("C1", "C2"))),
   assign_size_clust = list(C1 = 100, C2 = 100),
   gt_size_de_gene = list(GENE1 = 50, GENE2 = 50)
  )
  expect_silent(confusion_from_overlap_matrix(
     sim_res_DE, overlap_matrix,
     cut_off_true = 0.1,
     cut_off_false = 0.05,
     sig_threshold = 0.1,
     overlap_type = "union"
    )
  )
})

# Example with mock data
test_that("Combining Pvalues",{
  de_pvals_by_cluster <- list(
   cluster1 = data.frame(gene = c("A", "B"), cluster_id = "cluster1", p_val = c(0.01, 0.05)),
   cluster2 = data.frame(gene = c("A", "B"), cluster_id = "cluster2", p_val = c(0.02, 0.03))
  )
  expect_silent(result<-combine_de_pvals_by_cluster(de_pvals_by_cluster, method = "adjust_all"))
  
})