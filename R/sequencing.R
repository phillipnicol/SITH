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
#' @details The procedure is exactly the same as \code{\link{singleCell}} except that it allows multiple cells
#' to be sequenced at once (chosen randomly throughout the entire tumor). 
#' 
#' @author Phillip B. Nicol <philnicol740@gmail.com>
#' 
#' @examples
#' out <- simulateTumor(N = 1000)
#' df <- randomSingleCells(tumor = out, ncells = 5, noise = 0.10) 

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
#' @return A data frame with 1 row and columns corresponding to the mutations present in the cell. A 1 indicates that
#' the mutation is detected while a 0 indicates the mutation is not detected. 
#' 
#' @details This function selects the cell at \code{pos} (error if no cell at specified position exists) and returns
#' the list of mutations present in the cell. Due to technological artifacts, the false negative rate can be quite higher
#' (10-20 percent). To account for this,
#' the \code{noise} parameter introduces false negatives into the data set at the specified rate. 
#' 
#' @author Phillip B. Nicol <philnicol740@gmail.com>
#' 
#' @examples 
#' set.seed(1126490984)
#' out <- simulateTumor(N = 1000)
#' df <- singleCell(tumor = out, pos = c(0,0,0), noise = 0.1)
#' 
#' @references 
#' K. Jahn, J. Kupiers and N. Beerenwinkel. Tree inference for single-cell data. Genome Biology, volume 17, 2016. 
#' 
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

#' Simulate bulk sampling 
#' 
#' @description Simulate bulk sequencing data by takign a local sample from the tumor
#' and computing the variant allele frequencies of the various mutations. 
#' 
#' @param tumor A list which is the output of \code{\link{simulateTumor}}.
#' @param nsamples The number of bulk samples to take. 
#' @param cube.length The side length of the cube of cells to be sampled. 
#' @param threshold Only mutations with an allele frequency greater than the threshold will be included in the sample.
#' 
#' @return A data frame with \code{nsamples} rows and columns corresponding to the mutations. 
#' The entries are the mutation allele frequency.
#' 
#' @details This is the same as \code{\link{bulkSample}}, except multiple samples are taken 
#' with random center points. 
#' 
#' @examples 
#' out <- simulateTumor(N = 1000)
#' df <- randomBulkSamples(tumor = out, nsamples = 5, cube.length = 5, threshold = 0.05)
#' 
#' @author Phillip B. Nicol 
randomBulkSamples <- function(tumor, nsamples, cube.length = 5, threshold = 0.05) {
  if(cube.length %% 2 == 0 | cube.length < 1) {
    stop("cube.length must be an odd positive integer.")
  }
  
  cells <- sample(1:nrow(tumor$cell_ids), nsamples, replace = F)
  
  df <- data.frame(matrix(nrow = nsamples, ncol = 0))
  
  #Take around 1 percent of the tumor
  bulk_size <- cube.length^3
  cube.length <- cube.length - 1
  
  counter <- 1
  for(i in cells) {
    cx <- tumor$cell_ids[i,1]; cy <- tumor$cell_ids[i,2]; cz <- tumor$cell_ids[i,3]
    cell_subset <-  tumor$cell_ids[((tumor$cell_ids$x >= cx - cube.length/2) & (tumor$cell_ids$x <= cx + cube.length/2) &
                                  (tumor$cell_ids$y >= cy - cube.length/2) & (tumor$cell_ids$y <= cy + cube.length/2) &
                                  (tumor$cell_ids$z >= cz - cube.length/2) & (tumor$cell_ids$z <= cz + cube.length/2)),]
    
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

#' Simulate bulk sampling 
#' 
#' @description Simulate bulk sequencing data by takign a local sample from the tumor
#' and computing the variant allele frequencies of the various mutations. 
#' 
#' @param tumor The output of \link[TumorGenerator]{simulateTumor}. 
#' @param pos The center point of the sample.
#' @param cube.length The side length of the cube of cells to be sampled. 
#' @param threshold Only mutations with an allele frequency greater than the threshold will be included in the sample.
#' 
#' @return A data frame with 1 row and columns corresponding to the mutations. The entries are the mutation allele frequency.
#' 
#' @details A local region of the tumor is sampled by constructing a cube with side length \code{cube.length} around
#' the center point \code{pos}. Each cell within the cube is sampled, and the reported quantity is variant (or mutation) 
#' allele frequency. Lattice sites without cells are assumed to be normal tissue, and thus the reported MAF may be less than
#' 1.0 even if the mutation is present in all cancerous cells. 
#' 
#' Note that \code{cube.length} is required to be an odd integer (in order to have a well-defined center point). 
#' 
#' @author Phillip B. Nicol 
#' 
#' @examples 
#' set.seed(116776544, kind = "Mersenne-Twister", normal.kind = "Inversion")
#' out <- simulateTumor(N = 1000)
#' df <- bulkSample(tumor = out, pos = c(0,0,0))
#' 
#' @references 
#' K. Chkhaidze, T. Heide, B. Werner, M. Williams, W. Huang, G. Caravagna, T. Graham, and 
#' A. Sottoriva. Spatially con- strained tumour growth affects the 
#' patterns of clonal selection and neutral drift in cancer genomic data. PLOS Computational Biology, 2019.
bulkSample <- function(tumor, pos, cube.length = 5, threshold = 0.05) {
  if(length(pos) != 3) {
    stop("Position must be a vector of length 3.")
  }
  
  cx <- pos[1]; cy <- pos[2]; cz <- pos[3]
  bulk_size <- cube.length^3
  cube.length <- cube.length - 1
  
  df <- data.frame(matrix(nrow = 1, ncol = 0))  
  
  cell_subset <-  tumor$cell_ids[((tumor$cell_ids$x >= cx - cube.length/2) & (tumor$cell_ids$x <= cx + cube.length/2) &
                                    (tumor$cell_ids$y >= cy - cube.length/2) & (tumor$cell_ids$y <= cy + cube.length/2) &
                                    (tumor$cell_ids$z >= cz - cube.length/2) & (tumor$cell_ids$z <= cz + cube.length/2)),]
  
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



#' Simulate needle sequencing data
#' 
#' @description Simulate a sampling procedure which takes a fine needle through the simulated tumor and
#' reports the mutation allele frequency of the sampled cells. 
#' 
#' @param tumor A list which is the output of \code{\link{simulateTumor}}.
#' @param nsamples The number of samples to take.
#' @param threshold Only mutations with an allele frequency greater than the threshold will be included in the sample.
#' 
#' @author Phillip B. Nicol
#' 
#' @details This sampling procedure is inspired by Chkhaidze et. al. (2019). A random one-dimensional cross-section 
#' of the tumor is chosen, and the cells within this cross section are sampled, reporting mutation allele frequency. 
#' 
#' @examples 
#' out <- simulateTumor(N = 1000)
#' df <- randomNeedles(tumor = out, nsamples = 5)
#' 
#' @references 
#' K. Chkhaidze, T. Heide, B. Werner, M. Williams, W. Huang, G. Caravagna, T. Graham, and 
#' A. Sottoriva. Spatially con- strained tumour growth affects the 
#' patterns of clonal selection and neutral drift in cancer genomic data. PLOS Computational Biology, 2019.
randomNeedles <- function(tumor, nsamples, threshold = 0.05) {
  cells <- sample(1:nrow(tumor$cell_ids), nsamples, replace = F)
  
  df <- data.frame(matrix(nrow = nsamples, ncol = 0))
  
  counter <- 1
  for(i in cells) {
    dim <- sample(1:3,2,replace = F)
    cell_subset <- tumor$cell_ids[(tumor$cell_ids[,dim[1]] == tumor$cell_ids[i,dim[1]])
                                  & (tumor$cell_ids[,dim[2]] == tumor$cell_ids[i,dim[2]]),]
      
    input <- list()
    input$cell_ids <- cell_subset
    input$alleles <- tumor$alleles
    total_sqnc <- as.data.frame(randomSingleCells(input, nrow(cell_subset)))
    total_sqnc <- colSums(total_sqnc)/nrow(cell_subset)
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

