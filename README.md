[![Build Status](https://travis-ci.org/phillipnicol/SMITH.svg?branch=master)](https://travis-ci.org/phillipnicol/SMITH)

# A Spatial model of Intra-Tumor Heterogeneity (SITH) 

**An R package for visualizing and analyzing a spatial model of intratumor heterogeneity**

## Installation 

The following is required for installing the package:

  - `R` version 4.0.0 or newer.
  - Package `Rcpp` which can be downloaded from CRAN by running `install.packages("Rcpp")` from the console.
  - Package `rgl` is **strongly recommended** for visualizations can be downloaded from CRAN by running `install.packages(rgl)` from the console. 
    - macOS users may have to install [Xquartz](https://www.xquartz.org) before installing `rgl`. 
    
See the package documentation (LINK) for a complete list of dependencies.
    
As the package is on CRAN, it can be installed by running `install.packages(SITH)` from the command line. 

Alternatively, the package can be installed directly from this repository by first installing `devtools` (run `install.packages(devtools)`) and then running `install_github("phillipnicol/SITH")` from the console. 

## Features

  - Contains a 3D simulator of spatial tumor growth and mutation, similar to the model described in [this paper](https://www.nature.com/articles/nature14971).
    - The simulator is written in C++ since it is computationally expensive. 
    - A tumor with 1 million cells can be simulated in under a minute on a standard desktop computer.
  - 3D interactive visualizations of the simulated tumor using `rgl`. 
    - Option to color regions with high mutation red and regions with low mutation blue. 
    - 2D cross section
  - Summary of the spatial distribution of mutations throughout the tumor.
    - Creates graphs that show how genetic diversity changes in different spatial regions.
  - Create synthetic bulk and single-cell sequencing data from the simulated tumor.
    - Users can specify what part of the tumor the samples are taken from. 
