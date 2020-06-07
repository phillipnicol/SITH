\name{SMITH-package}
\alias{SMITH-package}
\alias{SMITH}
\docType{package}
\title{Visualize and analyze intratumor heterogeneity using a spatial model of tumor growth}
\description{The \pkg{SMITH} (spatial model of intratumor heterogeneity) package implements a lattice based spatial model
of tumor growth and mutation. Interactive 3D visualization of the tumor are possible using \pkg{rgl}. Additional functions for visualization and investigating the spatial distribution of mutants are provided. \pkg{SMITH} also
provides functions to simulate single cell sequencing and bulk sampling data sets from the simulated tumor. 
}
\section{Background}{
  On-lattice models of tumor growth and mutation are computationally efficient and provide a simple setting to study how spatially constrained growth impacts intratumor heterogeneity. While this model has been studied extensively in literature (see Waclaw (2015), Chkhaidze (2019), Opasic (2019)), existing software is either not publically available or inconvenient to use with \R. 
  
  The motivation for creating the \pkg{SMITH} package was to provide a spatial simulator that is both easy to use and can be used entirely with \R. The core function in the package is \code{\link{simulateTumor}}, which wraps a C++ implementation of the model into \R using the \pkg{Rcpp} package. Once the results of the simulation are saved as an \R object, \pkg{SMITH} provides several other useful functions for studying this model. 
  
  See the package vignette for more information on the model and the algorithm used. 
}

\author{Phillip Nicol}
\references{
B. Waclaw, I. Bozic, M. Pittman, R. Hruban, B. Vogelstein and M. Nowak. A spatial model predicts that dispersal and cell turnover limit intratumor heterogeneity. \emph{Nature}, pages 261-264, 2015. 

K. Chkhaidze, T. Heide, B. Werner, M. Williams, W. Huang, G. Caravagna, T. Graham, and A. Sottoriva. Spatially constrained tumour growth affects the patterns of clonal selection and neutral drift in cancer genomic data. \emph{PLOS Computational Biology}, 2019.

L. Opasic, D. Zhou, B. Wener, D. Dingli and A. Traulsen. How many samples are needed to infer truly clonal mutations
from heterogenous tumours? \emph{BMC Cancer}, https://doi.org/10.1186/s12885-019-5597-1. 
}
\examples{
\dontrun{                  
  #Simulate tumor growth
  out <- simulateTumor()
  
  #3d interactive visualization using rgl 
  visualizeTumor(out)
  #or see regions with lots of mutants
  visualizeTumor(out, plot.type = "heat")
  
  #get a summary on the spatial dist. of mutants
  sp <- spatialDistribution(out)
  
  #simulate single cell sequencing
  Scs <- randomSingleCells(tumor = out, ncells = 5, noise = 0.1)
  
  #simulate bulk sampling
  Bulks <- randomBulkSamples(tumor = out, nsamples = 5)
}
}
