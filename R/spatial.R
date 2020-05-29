
spatialDistribution <- function(tumor, make.plot = TRUE) {
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
  
  #Now do jaccard similarity 
  N <- 500
  jaccard_mat <- matrix(0, nrow = N, ncol = 2)
  
  for(i in 1:N) {
    s <- sample(1:nrow(tumor$cell_ids), 2, replace = F)
    dist <- (tumor$cell_ids$x[s[1]] - tumor$cell_ids$x[s[2]])^2 + (tumor$cell_ids$y[s[1]] - tumor$cell_ids$y[s[2]])^2
    dist <- dist + (tumor$cell_ids$z[s[1]] - tumor$cell_ids$z[s[2]])^2
    dist <- sqrt(dist)
    
    allele1 <- tumor$cell_ids$allele[s[1]] + 1
    allele2 <- tumor$cell_ids$allele[s[2]] + 1
    
    row1 <- which(rownames(tumor$alleles) == sprintf("%d", allele1))
    row2 <- which(rownames(tumor$alleles) == sprintf("%d", allele2))
    
    v1 <- tumor$alleles[row1 ,-ncol(tumor$alleles)]
    v1 <- v1[v1 != -1]
    
    v2 <- tumor$alleles[row2 ,-ncol(tumor$alleles)]
    v2 <- v2[v2 != -1]
    
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
  
  if(make.plot) {make_plot(out, tumor)}
  
  return(out)
}

make_plot <- function(out.spatial, tumor) {
  par(mfrow=c(2,2))
  plot(out.spatial$mean_mutant[,1], out.spatial$mean_mutant[,2], pch = 4, col = "blue", xlab = "Euclid. distance",
       ylab = "Mean # of mutants per cell", main = "Mutations per cell")
  plot(out.spatial$jaccard[,1], out.spatial$jaccard[,2], pch = 4, col = "red", main = "Jaccard index comparison",
       xlab = "Euclid. distance between cells", ylab = "Jaccard index")
  plot(1:length(tumor$muts$MAF), tumor$muts$MAF, pch = 16, col = "green", xlab = "k-th largest clone", 
       ylab = "Mutation allele frequency", main = "Clone sizes")
  
}

