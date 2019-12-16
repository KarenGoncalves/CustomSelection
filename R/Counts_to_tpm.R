#' Transforms count data into Transcripts Per Million (TPM) data
#'
#' With the matrix of counts and the size of the genes/transcripts, it calculates the TPM.
#'
#' This function was modified from a gist from Slowkow (https://gist.github.com/slowkow/c6ab0348747f86e2748b).
#'
#' Here, we do not calculate the effective length.
#'
#' @param counts Data frame genes/features (in the rows) by samples (in the columns).
#' Rownames must be in the names of the features and in the same format as the names of the featureLength vector
#' @param featureLength Named numeric vector with the length of each gene/feature.
#' Names must be the names of the features and in the same format as the rownames of the data frame counts.
#' @export
#' @references https://gist.github.com/slowkow/c6ab0348747f86e2748b
#' @return List with: a data frame with TPM of genes (rows) per sample (columns)
#' and a vector with the names of the genes/features for which the TPM was calculated

Counts_to_tpm <- function(counts, featureLength) {

  # Ensure valid arguments.
  counts <- counts[rownames(counts) %in% names(featureLength), ]
  featureLength <- featureLength[names(featureLength) %in% rownames(counts)]

  effLen <- as.data.frame(featureLength)

  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate <- log(counts[, i]) - log(effLen)
    denom <- log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))

  # Copy the row and column names from the counts matrix.
  tpm <- as.data.frame(tpm)
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm, rownames(tpm)) # Returns a data frame
}
