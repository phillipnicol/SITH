#ifndef POSTPROC_H
#define POSTPROC_H

#include"simutils.h"

extern int x_dim, y_dim, z_dim; 
#define rgb_lb 0.09
#define rgb_ub 0.91

namespace PostProcessing {
    void write_results(std::vector<cell> &cells, std::vector<specie> &species, 
                    Rcpp::NumericMatrix &cell_coords, Rcpp::IntegerMatrix &species_dict, Rcpp::IntegerVector &muts);
    void write_phylo_tree(std::vector<std::vector<int> > &phylo_tree, Rcpp::IntegerMatrix &rphylo_tree);
    Rcpp::NumericMatrix get_color_scheme(std::vector<specie> &species);
}

#endif 