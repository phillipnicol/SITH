

void write_results(std::vector<cell> &cells, std::vector<specie> &species, Rcpp::NumericMatrix &cell_coords, Rcpp::IntegerMatrix &species_dict, Rcpp::IntegerVector &muts) {
    for(int i = 0; i < cells.size(); ++i) {
        cell_coords(i, 0) = cells[i].x - x_dim/2;
        cell_coords(i, 1) = cells[i].y - y_dim/2;
        cell_coords(i, 2) = cells[i].z - z_dim/2;
        cell_coords(i, 3) = cells[i].species.id;
        cell_coords(i, 4) = cells[i].species.genotype.size();
        cell_coords(i, 5) = sqrt(pow(cells[i].x - x_dim/2, 2) + pow(cells[i].y - y_dim/2, 2) + pow(cells[i].z - z_dim/2,2));
    }    

    //Write species results
    int nmuts;
    for(int i = 0; i < species.size(); ++i) {
        nmuts = species[i].count;
        for(int j = 0; j < species[i].genotype.size(); ++j) {
            species_dict(i, j) = species[i].genotype[j];
            muts[species[i].genotype[j]] += nmuts;
        }
        for(int j = species[i].genotype.size(); j < species_dict.ncol() - 1; ++j) {
            species_dict(i, j) = -1;
        }
        species_dict(i, species_dict.ncol() - 1) = nmuts;
    }
}

void write_phylo_tree(std::vector<std::vector<int> > &phylo_tree, Rcpp::IntegerMatrix &rphylo_tree) {
    for(int i = 0; i < phylo_tree[0].size(); ++i) {
        rphylo_tree(i,0) = phylo_tree[0][i];
        rphylo_tree(i,1) = phylo_tree[1][i];
    }
}

Rcpp::NumericMatrix get_color_scheme(std::vector<specie> &species) {
    Rcpp::NumericMatrix color_scheme(3, species.size());

    color_scheme(0,0) = 0.5; color_scheme(1,0) = 0.5; color_scheme(2,0) = 0.5; 

    for(int i = 0; i < species.size(); ++i) {
        color_scheme(0,i) = R::runif(rgb_lb, rgb_ub);
        color_scheme(1,i) = R::runif(rgb_lb, rgb_ub);
        color_scheme(2,i) = R::runif(rgb_lb, rgb_ub);
    }
    return(color_scheme);
}