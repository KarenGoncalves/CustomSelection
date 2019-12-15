#' CustomReferences
#'
#' This function uses the Counts_to_tpm and the DAFS function to select the reference genes.
#' After transforming the counts into TPM values, the tpm data frame is used
#' as input for DAFS function.
#' We then calculate the average cutoff (result from DAFS) for all samples
#' Genes with average expression level lower than the average cutoff are excluded
#' and the genes with lowest covariance are selected.
#' @param counts matrix of genes/features (in the rows) by samples (in the columns)
#' @param featureLength a named numeric vector with the length of each gene/feature
#' @param top_genes percentage of genes (left after filtering) to be selected as references, default is 0.5\%
#' @export
#' @return data frame with genes as rows and two columns: Mean (average TPM) and Covariance
#'
customReferences <- function(counts, featureLength, top_genes = 0.5) {

  tpm <- Counts_to_tpm(counts, featureLength)

  # DAFS returns a vector with cutoffs for each sample (columns from counts) in log2
  # Here we obtain an average cutoff and transform it to use it with the tpm data
  cutoff <- 2 ** mean(DAFS(tpm))

  # Mean and covariance of the data in TPM are calculated
  reference_tpm <-
    data.frame(cbind("Mean" = apply(tpm, 1, mean),
                     "Covariance" =
                       apply(tpm, 1, sd)/apply(tpm, 1, mean)))

  # Exclude genes with average TPM lower than the cutoff
  reference_tpm <- reference_tpm[reference_tpm$Mean > cutoff,]

  # Here we sort the genes by level of covariance and select
  # the top genes with lowest covariance

  top_genes <- length(reference_tpm$Mean) * (top_genes/100)
  sorted <- reference_tpm[order(reference_tpm$Covariance), ]

  customreferences <- as.data.frame(sorted[1:top_genes, ] )
  return(customreferences)
}
