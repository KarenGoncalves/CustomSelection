#' Uses average TPM values and the covariance of TPM values to select reference genes from RNAseq data
#'
#' If counts_to_tpm and DAFS functions were already computed, this function will use their results to select the genes with lowest covariance, among those considered as expressed according to DAFS, as references.
#'
#' If there is no interest on keeping the results from DAFS or the TPM values for all genes in each sample, run customReferences instead.
#' @param countsToTpm_result Result of the counts_to_tpm function or data frame with other expression unit (FPKM, RPKM, CPM)
#' @param dafs_result  Result of DAFS fucntion: vector with threshold for noise/true expression for each sample (columns from tpm data frame) in log2
#' @param top_genes Percentage of genes (left after filtering) to be selected as references, default is 0.5\%.
#' @export
#' @return Data frame with genes as rows and two columns: Mean (average TPM) and Covariance.
#' @name gene_selection

gene_selection <- function(countsToTpm_result, dafs_result, top_genes = 0.5) {

  # DAFS returns a vector with cutoffs for each sample (columns from counts) in log2
  # Here we obtain an average cutoff and transform it to use it with the tpm data
  cutoff <- 2 ** mean(dafs_result)

  # Mean and covariance of the data in TPM are calculated
  mean_cov <-
    data.frame(cbind("Mean" = apply(countsToTpm_result, 1, mean),
                     "Covariance" =
                       apply(tpm, 1, sd)/apply(countsToTpm_result, 1, mean)))

  # Exclude genes with average TPM lower than the cutoff
  mean_cov <- mean_cov[mean_cov$Mean > cutoff,]

  # Here we sort the genes by level of covariance and select
  # the top genes with lowest covariance

  top_genes <- length(mean_cov$Mean) * (top_genes/100)
  sorted <- mean_cov[order(mean_cov$Covariance), ]

  customreferences <- as.data.frame(sorted[1:top_genes, ] )
  return(customreferences)
}
