
#' Quantify the spatial distribution of mutants
#' 
#' @description Provides a summary the spatial distribution of mutants within the simulated
#' tumor. 
#' 
#' @param tumor A list which is the output of \code{\link{simulateTumor}()}.
#' @param N The number of pairs to sample. 
#' @param cutoff For a plot of clone sizes, all mutations with a MAF below \code{cutoff} are ignored. 
#' @param make.plot Whether or not to make plots. 
#' 
#' @return A list with the following components
#' \itemize{
#' \item \code{mean_mutant} - A data frame with 2 columns giving the mean number of mutants
#' as a function of Euclidean distance from the lattice origin (Euclid. distance rounded to nearest integer). 
#' \item \code{mean_driver} - The same as \code{mean_mutant} except for driver mutations only. Will be \code{NULL} if 
#' no drivers are present in the simulated tumor. 
#' \item \code{jaccard} A data farme with two columns giving mean jaccard index as a function of Euclidean distance 
#' between pairs of cells (rounded to nearest integer). 
#' }
#' 
#' @details The genotype of a cell can be interpreted as a binary vector where the \eqn{i}-th component is 1 if mutation
#' \eqn{i} is present in the cell and is 0 otherwise. Then a natural comparison of the similarity between two cells is the 
#' Jaccard index \eqn{J(A,B) = |I(A,B)|/|U(A,B)|}, where \eqn{I(A,B)} is the intersection of \eqn{A} and \eqn{B} and 
#' \eqn{U(A,B)} is the union. This function estimates the Jaccard index as a function of Euclidean distance between the 
#' cells by randomly sampling \eqn{N} pairs of cells. 
#' 
#' @examples 
#' set.seed(1126490984)
#' out <- simulateTumor(N = 1000, du = 0.1)
#' sp <- spatialDistribution(tumor = out, make.plot = FALSE)
#' 
#' @author Phillip B. Nicol
spatialDistribution <- function(tumor, N = 500, cutoff = 0.01, make.plot = TRUE) {
  out <- list()
  
  #First do mean number of mutants by distance
  #round first to get integers
  tumor$cell_ids$distances <- round(tumor$cell_ids$distance)
  
  vals <- c(0:max(tumor$cell_ids$distance))
  mean_mutant <- as.data.frame(sapply(vals, function(x) {
    return(mean(tumor$cell_ids[tumor$cell_ids$distance == x,5]))
  }))
  
  vals <- as.data.frame(vals)
  out$mean_mutant <- cbind(vals, mean_mutant)
  colnames(out$mean_mutant) <- c("Distance", "Mean # mutations")
  
  #Repeat process for drivers
  if(length(tumor$drivers) > 0) {
    driver_ids <- tumor$cell_ids[which(tumor$cell_ids$allele %in% tumor$drivers),]
    vals <- c(0:max(driver_ids$distance))
    mean_driver <- as.data.frame(sapply(vals, function(x) {
      return(nrow(driver_ids[driver_ids$distance == x,])/nrow(tumor$cell_ids[tumor$cell_ids$distance == x,]) )
    }))
    vals <- as.data.frame(vals)
    out$mean_driver <- cbind(vals, mean_driver)
    colnames(out$mean_driver) <- c("Distance", "Mean # drivers")
  }
  
  #Now do jaccard similarity 
  jaccard_mat <- matrix(0, nrow = N, ncol = 2)
  
  for(i in 1:N) {
    s <- sample(1:nrow(tumor$cell_ids), 2, replace = F)
    dist <- (tumor$cell_ids$x[s[1]] - tumor$cell_ids$x[s[2]])^2 + (tumor$cell_ids$y[s[1]] - tumor$cell_ids$y[s[2]])^2
    dist <- dist + (tumor$cell_ids$z[s[1]] - tumor$cell_ids$z[s[2]])^2
    dist <- sqrt(dist)
    
    allele1 <- tumor$cell_ids$allele[s[1]] + 1
    allele2 <- tumor$cell_ids$allele[s[2]] + 1
    
    row1 <- tumor$alleles[allele1,-ncol(tumor$alleles)]
    row2 <- tumor$alleles[allele2,-ncol(tumor$alleles)]

    v1 <- row1[row1 != -1]
    v2 <- row2[row2 != -1]
    
    jaccard_mat[i,1] <- round(dist)
    jaccard_mat[i,2] <- length(intersect(v1, v2))/length(union(v1,v2))
  }
  vals <- c(0:max(jaccard_mat[,1]))
  mean_jaccard <- as.data.frame(sapply(vals, function(x) {
    return(mean(jaccard_mat[jaccard_mat[,1] == x,2]))
  }))
  
  vals <- as.data.frame(vals)
  out$jaccard <- cbind(vals, mean_jaccard)
  colnames(out$jaccard) <- c("Distance", "Mean jaccard index")
  
  if(make.plot) {make_plot(out, tumor, cutoff)}
  
  return(out)
}

make_plot <- function(out.spatial, tumor, cutoff) {
  par(mfrow=c(2,2))
  
  plot(out.spatial$mean_mutant[,1], out.spatial$mean_mutant[,2], pch = 4, col = "blue", xlab = "Euclid. distance from origin",
       ylab = "Mean # of mutants per cell", main = "Mutations per cell")
  
  hist(tumor$cell_ids$nmuts, breaks = 0:max(tumor$cell_ids$nmuts), main = "Histogram of mutations per cell",
       xlab = "Number of mutations")
  
  plot(out.spatial$jaccard[,1], out.spatial$jaccard[,2], pch = 16, col = "red", main = "Jaccard index comparison",
       xlab = "Euclid. distance between cells", ylab = "Jaccard index")
  
  tumor$muts <- tumor$muts[order(tumor$muts$MAF, decreasing = T),]; tumor$muts <- tumor$muts[tumor$muts$MAF > cutoff,]
  plot(1:length(tumor$muts$MAF), tumor$muts$MAF, pch = 16, col = "green", xlab = "k-th largest clone", 
       ylab = "Mutation allele frequency", main = "Clone sizes")
  
}

