#' DAFS
#'
#' Modified from George and Chang (2014)
#' The function normally requires two inputs: data and name.
#' We do not export the result, so we do not need a 'name' argument.
#' @param tpm Data frame of TPM of genes (rows) per sample (columns)
#' @author George and Chang (2014)
#' @references doi:10.1186/1471-2105-15-92
#' @export
#' @return a vector with threshold for noise/true expression for each sample (columns from counts) in log2

DAFS <- function(tpm){

  # determine which rows have all 0 counts
  out <- apply(tpm, 1, function(x) all(x == 0))
  if(length(out[out == "TRUE"]) > 0) data <- tpm[-which(out == "TRUE"),]

  # set vector for cutoff values
  cutv <- rep(0, 0)

  for (i in 1:ncol(data)) {
    # specify array and remove 0 counts
    xx <- data[, i]
    xx <- xx[-which(xx == 0)]

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



