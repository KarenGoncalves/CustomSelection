#' Example - vector with gene length
#'
#' Length of Arabidopsis thaliana genes (TAIR10) obtained with the following code:
#' library(biomaRt)
#' \code{ath <- useMart('plants_mart', host = "plants.ensembl.org", dataset = "athaliana_eg_gene")
#' gene_start_end = getBM(attributes = c('ensembl_gene_id', 'start_position', 'end_position'), mart = ath)
#' featureLength <- gene_start_end$end_position - gene_start_end$start_position
#' names(featureLength) <- gene_start_end$ensembl_gene_id}

#' @format Name vector: length (values) of genes (names) of Arabidopsis thaliana
#' @source \url{plants.ensembl.org}
"ath_featureLength"
