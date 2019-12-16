#' Uses average TPM values and the covariance of TPM values to select reference genes from RNAseq data
#'
#' This function uses the Counts_to_tpm and the DAFS function to select the reference genes.
#'
#' After transforming the counts into TPM values, the tpm data frame is used as input for DAFS function.
#'
#' We then calculate the average cutoff (result from DAFS) for all samples.
#'
#' Genes with average expression level lower than the average cutoff are excluded and the genes with lowest covariance are selected.
#' @param counts Data frame genes/features (in the rows) by samples (in the columns). Rownames must be in the names of the features and in the same format as the names of the featureLength vector
#' @param featureLength Named numeric vector with the length of each gene/feature. Names must be the names of the features and in the same format as the rownames of the data frame counts.
#' @param top_genes Percentage of genes (left after filtering) to be selected as references, default is 0.5\%.
#' @export
#' @return Data frame with genes as rows and two columns: Mean (average TPM) and Covariance.
#'

customReferences <- function(counts, featureLength, top_genes = 0.5) {

  tpm <- Counts_to_tpm(counts, featureLength)

  # DAFS returns a vector with cutoffs for each sample (columns from counts) in log2
  # Here we obtain an average cutoff and transform it to use it with the tpm data
  cutoff <- 2 ** mean(DAFS(tpm))

  # Mean and covariance of the data in TPM are calculated
  mean_cov <-
    data.frame(cbind("Mean" = apply(tpm, 1, mean),
                     "Covariance" =
                       apply(tpm, 1, sd)/apply(tpm, 1, mean)))

  # Exclude genes with average TPM lower than the cutoff
  mean_cov <- mean_cov[mean_cov$Mean > cutoff,]

  # Here we sort the genes by level of covariance and select
  # the top genes with lowest covariance

  top_genes <- length(mean_cov$Mean) * (top_genes/100)
  sorted <- mean_cov[order(mean_cov$Covariance), ]

  customreferences <- as.data.frame(sorted[1:top_genes, ] )
  return(customreferences)
}
