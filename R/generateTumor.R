#' Spatial simulation of tumor growth
#' 
#' @description Simulate the spatial growth of a tumor using a multi-type branching
#' process on the three-dimensional integer lattice.
#' 
#' @param N Number of cells in the tumor
#' @param b Cell division rate
#' @param d Cell death rate 
#' @param u Mutation rate. When a cell divides, both daughter cell acquire $Pois(u)$ genetic alterations
#' @param du The probability that a genetic alteration is a driver mutation
#' @param s The selective advantage conferred to a driver mutation. A cell with k
#' driver mutations is given birth rate $bs^k$
#' @param verbose Whether or not to print simulation details to the R console
#' 
#' @return A list \code{out}. 
generateTumor <- function(N = 100000, b = 0.25, d = 0.13, u = 0.01, du = 0.003, s = 1.1, 
                          seed = -1, verbose = TRUE) {
  if(seed == -1) {seed <- sample((-2^{31}+1):(2^{31}-1), 1)}
  input <- list()
  input$params <- c(N, b, d, u, du, s, verbose, seed)
  tumor <- simulate_tumor(input)
  out <- list()
  out$cell_ids <- data.frame(tumor[[1]])
  colnames(out$cell_ids) <- c("x", "y", "z", "allele", "nmuts", "distance")
  
  out$alleles <- data.frame(tumor[[2]]); nc <- ncol(out$alleles)

  df <-  as.data.frame(tumor[[3]])
  #ix <- which(df[,1] > cutoff*N)
  ix <- as.data.frame(0:(nrow(df)-1))
  out$muts <- cbind(ix, df[,1], df[,1]/N)
  #out$muts <- as.data.frame(cbind(ix - 1, df[ix,], df[ix,]/N))
  #out$muts <- out$muts[order(-out$muts[,3]),]
  colnames(out$muts) <- c("id", "count", "MAF")
  
  out$phylo_tree <- as.data.frame(tumor[[4]])
  colnames(out$phylo_tree) <- c("parent", "child")
  
  color_scheme_mat <- tumor[[5]]
  out$color_scheme <- apply(color_scheme_mat, 2, function(x) {
    return(rgb(x[1],x[2],x[3],1))
  })
  
  out$drivers <- tumor[[7]]
  
  out$params <- input$params
  
  return(out)
}