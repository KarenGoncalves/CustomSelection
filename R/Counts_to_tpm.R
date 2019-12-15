#' Counts_to_tpm
#' This function was modified from a gist from Slowkow, we do not calculate the effective length
#' With the matrix of counts and the size of the genes/transcripts, it calculates the TPM
#' @param counts matrix of genes/features (in the rows) by samples (in the columns)
#' @param featureLength a named numeric vector with the length of each gene/feature
#' @export
#' @references https://gist.github.com/slowkow/c6ab0348747f86e2748b
#' @return list with: a data frame with TPM of genes (rows) per sample (columns)
#' and a vector with the names of the genes/features for which the TPM was calculated

Counts_to_tpm <- function(counts, featureLength) {

  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))

  effLen <- as.data.frame(featureLength)

  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  counts <- counts[idx, ]
  effLen <- effLen[idx, ]
  featureLength <- featureLength[idx]

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
