

findMut <- function(mut) {
  gtypes <- apply(out$genotypes[,-ncol(out$genotypes)], 1, function(r) any(r %in% mut))
  gtypes <- which(gtypes == TRUE)
  ixs <- which(out$cell_ids$genotype %in% gtypes)
  ixs
}