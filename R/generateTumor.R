#' @title Spatial simulation of tumor growth
#' 
#' @description Simulate the spatial growth of a tumor with a multi-type branching
#' process on the three-dimensional integer lattice.
#' 
#' @param N Number of cells in the tumor.
#' @param b Cell division rate.
#' @param d Cell death rate.
#' @param u Mutation rate. When a cell divides, both daughter cell acquire \eqn{Pois(u)} genetic alterations
#' @param du The probability that a genetic alteration is a driver mutation.
#' @param s The selective advantage conferred to a driver mutation. A cell with k
#' driver mutations is given birth rate \eqn{bs^k}.
#' @param verbose Whether or not to print simulation details to the R console.
#' 
#' @return A list with components 
#' \itemize{
#' \item \code{cell_ids} - A data frame containing the information for the simulated cells. (x,y,z) position, allele ID number
#' (note that 0 is the wild-type allele),
#' number of genetic alterations, and Euclidean distance from origin are included. 
#' \item \code{muts} - A data frame consisting of the mutation ID number, the count of the mutation within the population, 
#' and the mutation allele frequency (which is the count divided by N). 
#' \item \code{phylo_tree} - A data frame giving all of the information necessary to determine the order of mutations. The parent
#' of a mutation is defined to be the most recent mutation that precedes it. Since the ID 0 corresponds to the initial mutation,
#' 0 does not have any parents and is thus the root of the tree. 
#' \item \code{alleles} - A data frame containing the information about the mutations that make up each allele. The \eqn{i}-th 
#' row of this data frame corresponds to the allele ID \eqn{i-1}. The positive numbers in each row correspond to the IDs of the 
#' mutations present in that allele, while a -1 is simply a placeholder and indicates no mutation. The count column gives
#' the number of cells which have the specific allele. 
#' \item \code{color_scheme} - A vector containing an assignment of a color to each allele.
#' \item \code{drivers} - A vector containing the ID numbers for the driver mutations.
#' \item \code{time} - The simulated time (in days). 
#' \item \code{params} - The parameters used for the simulation. 
#' } 
#' 
#' @details The model is based upon Waclaw et. al. (2015), although the simulation algorithm used is different. A growth of a canerous tumor
#' is modeled using an exponential birth-death process on the three-dimensional integer lattice. Each cell is given a birth rate
#' \eqn{b} and a death rate \eqn{d} such that the time until cell division or cell death is exponentially distributed with 
#' parameters \eqn{b} and \eqn{d}, respectively. A cell can replicate if at least one of the six sites adjacent to it is
#' unoccupied. Each time cell replication occurs, both daughter cells recevie \eqn{Pois(u)} genetic alterations. Each 
#' alteration is a driver mutation with some probability \eqn{du}. A cell with k driver mutations is given birth rate 
#' \eqn{bs^k}. The simulation begins with a single cell at the origin at time \eqn{t = 0}. 
#' 
#' The model is simulated using a Gillespie algorithm. See the package vignette for details on how the algorithm is implemented. 
#' 
#' @author Phillip B. Nicol <philnicol740@gmail.com>
#' 
#' @references B. Waclaw, I. Bozic, M. Pittman, R. Hruban, B. Vogelstein and M. Nowak. A spatial model predicts
#' that dispersal and cell turnover limit intratumor heterogeneity. \emph{Nature}, pages 261-264, 2015. 
#' 
#' D. Gillespie. Exact stochastic simulation of coupled chemical reactions. \emph{The Journal of Physical Chemistry},
#' volume 81, pages 2340-2361, 1970.
#' 
simulateTumor <- function(N = 250000, b = 0.25, d = 0.18, u = 0.01, du = 0.003, s = 1.05, verbose = TRUE) {
  #create input list
  input <- list()
  input$params <- c(N, b, d, u, du, s, verbose)
  
  #call Rcpp function
  tumor <- simulate_tumor(input)
  
  out <- list()
  
  #position data for the N cells 
  out$cell_ids <- data.frame(tumor[[1]])
  colnames(out$cell_ids) <- c("x", "y", "z", "allele", "nmuts", "distance")
  
  #record the information for the uniqe alleles
  out$alleles <- data.frame(tumor[[2]]); nc <- ncol(out$alleles)
  colnames(out$alleles)[ncol(out$alleles)] <- "count"

  #get mutation ID and MAF 
  df <-  as.data.frame(tumor[[3]])
  ix <- as.data.frame(0:(nrow(df)-1))
  out$muts <- cbind(ix, df[,1], df[,1]/N)
  colnames(out$muts) <- c("id", "count", "MAF")
  
  #record phylogenetic tree 
  out$phylo_tree <- as.data.frame(tumor[[4]])
  colnames(out$phylo_tree) <- c("parent", "child")
  
  #Record colors for plotting purposes 
  color_scheme_mat <- tumor[[5]]
  out$color_scheme <- apply(color_scheme_mat, 2, function(x) {
    return(rgb(x[1],x[2],x[3],1))
  })
  
  #record a list of drivers and the simulated time (in days)
  out$drivers <- tumor[[7]]
  out$time <- tumor[[8]]
  
  #return parameters used for simulation 
  out$params <- input$params
  
  return(out)
}