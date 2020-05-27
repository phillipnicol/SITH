

generateTumor <- function(N = 100000, b = 0.25, d = 0.13, u = 0.01, du = 0.00003, s = 1.1, cutoff = 0.01) {
  input <- list()
  input$params <- c(N, b, d, u, du, s)
  tumor <- simulate_tumor(input)
  out <- list()
  out$cell_ids <- data.frame(tumor[[1]])
  colnames(out$cell_ids) <- c("x", "y", "z", "allele", "nmuts", "distance")
  
  out$species_dict <- data.frame(tumor[[2]]); nc <- ncol(out$species_dict)
  out$species_dict <- out$species_dict[out$species_dict[,nc] > 0,]
  out$species_dict <- out$species_dict[order(-out$species[,nc]),]

  df <-  as.data.frame(tumor[[3]])
  ix <- which(df[,1] > cutoff*N)
  out$muts <- cbind(ix - 1, df[ix,], df[ix,]/N)
  out$muts <- out$muts[order(-out$muts[,3]),]
  colnames(out$muts) <- c("id", "count", "MAF")
  
  return(out)
}