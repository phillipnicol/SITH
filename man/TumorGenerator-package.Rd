\name{TumorGenerator-package}
\alias{TumorGenerator-package}
\alias{TumorGenerator}
\docType{package}
\title{
  A short title line describing what the package does
}
\description{
  A more detailed description of what the package does. A length
  of about one to five lines is recommended.
}
\details{
  This section should provide a more detailed overview of how to use the
  package, including the most important functions.
}
\author{
Your Name, email optional.

Maintainer: Your Name <your@email.com>
}
\references{
K. Chkhaidze, T. Heide, B. Werner, M. Williams, W. Huang, G. Caravagna, T. Graham, and A. Sottoriva. 
Spatially con- strained tumour growth affects the patterns of clonal selection and neutral drift in cancer genomic data. 
\emph{PLOS Computational Biology}, 2019.
}
\keyword{ package }
\seealso{
  Optional links to other man pages
}
\examples{
  \dontrun{
    ## Simulate a tumor with 100000 cells
    out <- simulateTumor(N = 100000)
    
    ## Interactive visualization of tumor (if rgl is installed)
    visualizeTumor(tumor_instance = out)
    
    ## View a summary of the spatial distribution of mutations
    sp <- spatialDistribution(tumor = out)
  }
}
