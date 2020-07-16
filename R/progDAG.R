

progressionChain <- function(n) {
  G <- matrix(0, nrow = n, ncol = 4)
  G[,1] <- c(0:(n-1))
  G[,2] <- c(1:n)
  G[,3] <- 10^{-3}
  G[,4] <- 1.0
  
  colnames(G) <- c("Head", "Tail", "Mut. rate", "Selective advantage")
  
  return(G)
}

progressionDAG_from_igraph <- function(iG) {
  if(!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' required to use this function!")
  }  
  if(!is.igraph(iG)) {
    stop("Input must be igraph object!")
  }
  if(!is.dag(iG)) {
    stop("Graph should be directed and acyclic!")
  }
  
  Vcount <- gorder(iG)
  
  Edgelist <- as_edgelist(iG,names = FALSE)
  Edgelist <- Edgelist-1
  
  G <- matrix(0, nrow = nrow(Edgelist), ncol = 4)
  G[,c(1,2)] <- Edgelist
  G[,3] <- 10^{-3}
  G[,4] <- 1
  
  colnames(G) <- c("Head", "Tail", "Mut. rate", "Selective advantage")
  
  return(G)
}