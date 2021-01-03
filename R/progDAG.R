
#' @title Create a linear chain graph to describe the order of mutations
#' 
#' @description A helper function for \code{\link{simulateTumor}()} which returns to the user
#' the edge list for a linear chain. 
#' 
#' @param n Number of vertices in the chain
#' 
#' @return A matrix with 4 columns and n-1 rows which will be accepted as input to
#' \code{\link{simulateTumor}()}.
#' 
#' @author Phillip B. Nicol <philnicol740@gmail.com>
#' 
#' @examples 
#' G <- progressionChain(3)
progressionChain <- function(n) {
  G <- matrix(0, nrow = n, ncol = 4)
  G[,1] <- c(0:(n-1))
  G[,2] <- c(1:n)
  G[,3] <- 10^{-3}
  G[,4] <- 1.0
  
  colnames(G) <- c("Head", "Tail", "Mut. rate", "Selective advantage")
  
  return(G)
}

#' @title Define the progression of mutations from an \code{igraph} object 
#' 
#' @description A helper function for \code{\link{simulateTumor}()} which returns to the user
#' the edge list for a DAG which is defined as an \code{igraph} object.
#' 
#' @param iG An igraph object for a directed acyclic graph. 
#' 
#' @return A matrix with 4 columns which contains the edges of the graph as well as the rate of 
#' crossing each edge and the selective advantage/disadvantage obtained by crossing each edge.
#' 
#' @author Phillip B. Nicol <philnicol740@gmail.com>
#' 
progressionDAG_from_igraph <- function(iG) {
  if(!requireNamespace("igraph", quietly = TRUE)) {
    stop("Package 'igraph' required to use this function!")
  }  
  if(!igraph::is.igraph(iG)) {
    stop("Input must be igraph object!")
  }
  if(!igraph::is.dag(iG)) {
    stop("Graph should be directed and acyclic!")
  }
  
  Vcount <- igraph::gorder(iG)
  
  Edgelist <- igraph::as_edgelist(iG,names = FALSE)
  Edgelist <- Edgelist-1
  
  G <- matrix(0, nrow = nrow(Edgelist), ncol = 4)
  G[,c(1,2)] <- Edgelist
  G[,3] <- 10^{-3}
  G[,4] <- 1
  
  colnames(G) <- c("Head", "Tail", "Mut. rate", "Selective advantage")
  
  return(G)
}

#check the input G for correctness 
checkG <- function(G) {
  if(ncol(G) != 4) {
    stop("G must have 4 columns. See progressionChain() for an example.")
  }
  
  m <- as.vector(G[,c(1,2)])
  if(!is.numeric(m)) {
    stop("Head and tail of G must be numeric. See progressionChain() for an example.")
  }
  
  if(any(G[,3] < 0) | any(G[,3] > 1)) {
    stop("Mutation rate is a probability and should be between 0 and 1.")
  }
  
  if(any(G[,4] < 0)) {
    stop("Selective advantages cannot be negative.")
  }
}