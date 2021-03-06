#' Calculates the threshold for a gene to be considered truly expressed
#'
#' This function calculates the threshold for a gene to be considered truly expressed in each sample (columns of the expression data frame).
#'
#' Modified from George and Chang (2014).
#'
#' The function normally requires two inputs: data and name.
#' We do not export the result, so we do not need a 'name' argument.
#' @param tpm Data frame of log2 transformed expression values (RPKM, FPKM, TPM, CPM, etc.). We recomend the use of the fuction Counts_to_tpm from this package. As the Counts_to_tpm function returns a list, use as input for DAFS the first object of its result, e.g.:
#'
#' \code{DAFS(Counts_to_tpm(...)[[1]])}
#' @author George and Chang (2014)
#' @references doi:10.1186/1471-2105-15-92
#' @export
#' @return A vector with threshold for noise/true expression for each sample (columns from tpm data frame) in log2
#' @name DAFS
#' @examples
#' data("sample_counts"); data("ath_featureLength")
#' tpm = Counts_to_tpm(counts = sample_counts, featureLength = ath_featureLength)
#' DAFS(tpm[[1]])
#'

DAFS <- function(tpm){

  # determine which rows have all 0 counts
  out <- apply(tpm, 1, function(x) all(x == 0))
  if(length(out[out == "TRUE"]) > 0){ data <- tpm[!out,]} else {data = tpm}

  # set vector for cutoff values
  cutv <- rep(0, 0)

  for (i in 1:ncol(data)) {
    # specify array and remove 0 counts
    xx <- data[, i]
    xx <- xx[xx != 0]

    # take log2 of data
    log2xx <- log2(xx)
    dlog2 <- data.frame(LogC = log2xx)

    #vector to store Kolmogorov Smirnov distance statistics
    vv <- rep(0, 0)

    # select start point
    start <- length(log2xx[log2xx == min(log2xx)])/length(log2xx)

    # set sequence
    s <- seq(round(start, 2), 0.5, by = 0.005)

    # loop through cuts of the data to determine targeted K-S statistic
    for(q in s) {

      # select data greater than a quantile and run Mclust on that data to determine theoretical distribution
      d <- log2xx[which(log2xx > quantile(log2xx, q, na.rm = T))]
      out <- mclust::Mclust(d, G = 1)
      ks <- ks.test(d, "pnorm", out$parameter$mean,
                    out$parameter$variance$sigmasq)
      vv <- c(vv, ks$statistic)
    }

    # determine first left-most local minima
    out <- earth::earth(s, vv, thresh = 0.005)

    # save suggested cut
    cutv <- c(cutv, min(out$cuts[out$cuts > 0]))
  }

  names(cutv) <- colnames(data)
  return(cutv)
}
