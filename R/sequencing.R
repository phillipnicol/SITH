#' Simulate single cell sequencing data 
#' 
#' @description Simulate single cell sequencing data by random selecting cells from the tumor. 
#' 
#' @param tumor A list which is the output of \code{\link{simulateTumor}}.
#' @param ncells The number of cells to sample.
#' @param noise The false negative rate. 
#' 
#' @return A data frame with sample names on the row and mutation ID on the column. 
#' A 1 indicates that the mutation is present in the cell and a 0 indicates the mutation is not present. 
#' 
#' @details TODO 
#' 
#' @author Phillip B. Nicol 

randomSingleCells <- function(tumor, ncells, noise = 0.0) {
  cells <- sample(1:nrow(tumor$cell_ids), ncells, replace = F)
  
  df <- data.frame(matrix(nrow = ncells, ncol = 0))
  
  counter <- 1
  for(i in cells) {
    allele <- tumor$cell_ids[i,4] + 1
    row <- tumor$alleles[allele,]
    df[counter,] <- rep(0, ncol(df))
    rownames(df)[counter] <- sprintf("SC-%d", counter)
    for(j in row) {
      if(j == -1) {
        break
      }
      
      #put cell into matrix 
      if(sprintf("mutID-%d", j) %in% colnames(df)) {
        key <- which(colnames(df) == sprintf("mutID-%d", j))
        if(rbinom(1,1,noise) == 0) {
          df[counter, key] <- 1
        }
      }
      else {
        df <- cbind(df, 0)
        colnames(df)[ncol(df)] <- sprintf("mutID-%d", j)
        if(rbinom(1,1,noise) == 0) {
          df[counter, ncol(df)] <- 1          
        }
      }
    }
    counter <- counter + 1
  }
  
  return(df)
}

#' Simulate single cell sequencing data 
#' 
#' @description Simulate single cell sequencing data by selecting a cell at a specified position 
#' 
#' @param tumor A list which is the output of \code{\link{simulateTumor}}.
#' @param pos A vector of length 3 giving the (x,y,z) coordinates of the cell to sample. 
#' @param noise The false negative rate. 
#' 
#' @details TODO 
#' 
#' @author Phillip B. Nicol
#' 
#' @references 
#' 
singleCell <- function(tumor, pos, noise = 0.0) {
  if(length(pos) != 3) {
    stop("Position must be a vector of length 3.")
  }
  cell_df <- tumor$cell_ids[tumor$cell_ids$x == pos[1] & tumor$cell_ids$y == pos[2] & tumor$cell_ids$z == pos[3],]
  if(nrow(cell_df) == 0) {
    stop("No cell at specified position.")
  }
  
  df <- data.frame(matrix(nrow = 1, ncol = 0))

  allele <- cell_df[,4] + 1
  row <- tumor$alleles[allele,]
  
  for(j in row) {
    if(j == -1) {
      break
    }
    
    df <- cbind(df, 0)
    colnames(df)[ncol(df)] <- sprintf("mutID-%d", j)
    if(rbinom(1,1,noise) == 0) {df[1, ncol(df)] <- 1}
  }
  
  return(df)
}

#' Simulate bulk sequencing data 
#' 
#' @description Simulate bulk sequencing data by takign a local sample of the tumor
#' and computing the variant allele frequencies of the various mutations. 
#' 
#' @param tumor A list which is the output of \code{\link{simulateTumor}}.
#' @param nsamples The number of samples to take.
#' @param cube_length The side length of the cube in
#' @param threshold TODO 
randomBulkSamples <- function(tumor, nsamples, cube_length = 5, threshold = 0.05) {
  if(cube_length %% 2 == 0 | cube_length < 1) {
    stop("cube_length must be an odd positive integer.")
  }
  
  cells <- sample(1:nrow(tumor$cell_ids), nsamples, replace = F)
  
  df <- data.frame(matrix(nrow = nsamples, ncol = 0))
  
  #Take around 1 percent of the tumor
  bulk_size <- cube_length^3
  cube_length <- cube_length - 1
  
  counter <- 1
  for(i in cells) {
    cx <- tumor$cell_ids[i,1]; cy <- tumor$cell_ids[i,2]; cz <- tumor$cell_ids[i,3]
    cell_subset <-  tumor$cell_ids[((tumor$cell_ids$x >= cx - cube_length/2) & (tumor$cell_ids$x <= cx + cube_length/2) &
                                  (tumor$cell_ids$y >= cy - cube_length/2) & (tumor$cell_ids$y <= cy + cube_length/2) &
                                  (tumor$cell_ids$z >= cz - cube_length/2) & (tumor$cell_ids$z <= cz + cube_length/2)),]
    
    input <- list()
    input$cell_ids <- cell_subset
    input$alleles <- tumor$alleles
    total_sqnc <- as.data.frame(randomSingleCells(input, nrow(cell_subset)))
    total_sqnc <- colSums(total_sqnc)/bulk_size
    total_sqnc <- total_sqnc[total_sqnc > threshold]
    
    total_sqnc <- as.data.frame(t(total_sqnc))

    df[counter,] <- rep(0, ncol(df))
    rownames(df)[counter] <- sprintf("Bulk-%d", counter)    
    for(j in 1:ncol(total_sqnc)) {
      if(colnames(total_sqnc)[j] %in% colnames(df)) {
        key <- which(colnames(df) == colnames(total_sqnc)[j])
        df[counter, key] <- total_sqnc[1,j] 
      }
      else {
        df <- cbind(df, 0)
        colnames(df)[ncol(df)] <- colnames(total_sqnc)[j]
        df[counter, ncol(df)] <- total_sqnc[1,j]
      }
    }
    counter <- counter + 1
  }
  
  return(as.data.frame(df))
}

#' Quantify the spatial distribution of mutants
#' 
#' @param tumor The output of \link[TumorGenerator]{simulateTumor}. 
#' @param pos TODO
#' @param cube_length TODO 
#' @param threshold TODO 
bulkSample <- function(tumor, pos, cube_length = 5, threshold = 0.05) {
  if(length(pos) != 3) {
    stop("Position must be a vector of length 3.")
  }
  
  cx <- pos[1]; cy <- pos[2]; cz <- pos[3]
  bulk_size <- cube_length^3
  cube_length <- cube_length - 1
  
  df <- data.frame(matrix(nrow = 1, ncol = 0))  
  
  cell_subset <-  tumor$cell_ids[((tumor$cell_ids$x >= cx - cube_length/2) & (tumor$cell_ids$x <= cx + cube_length/2) &
                                    (tumor$cell_ids$y >= cy - cube_length/2) & (tumor$cell_ids$y <= cy + cube_length/2) &
                                    (tumor$cell_ids$z >= cz - cube_length/2) & (tumor$cell_ids$z <= cz + cube_length/2)),]
  
  if(nrow(cell_subset) == 0) {
    warning("Sample is empty")
    return(1)
  }
  
  input <- list()
  input$cell_ids <- cell_subset
  input$alleles <- tumor$alleles
  total_sqnc <- as.data.frame(randomSingleCells(input, nrow(cell_subset)))
  total_sqnc <- colSums(total_sqnc)/bulk_size
  total_sqnc <- total_sqnc[total_sqnc > threshold]
  
  total_sqnc <- as.data.frame(t(total_sqnc))
  
  for(j in 1:ncol(total_sqnc)) {
    if(colnames(total_sqnc)[j] %in% colnames(df)) {
      key <- which(colnames(df) == colnames(total_sqnc)[j])
      df[1, key] <- total_sqnc[1,j] 
    }
    else {
      df <- cbind(df, 0)
      colnames(df)[ncol(df)] <- colnames(total_sqnc)[j]
      df[1, ncol(df)] <- total_sqnc[1,j]
    }
  }

return(as.data.frame(df))
}