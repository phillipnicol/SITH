#include"postproc.h" 


void PostProcessing::write_results(std::vector<cell> &cells, std::vector<specie> &species, 
                                    Rcpp::NumericMatrix &cell_coords, Rcpp::IntegerMatrix &species_dict, Rcpp::IntegerVector &muts) {
    for(int i = 0; i < cells.size(); ++i) {
        cell_coords(i, 0) = cells[i].x - x_dim/2;
        cell_coords(i, 1) = cells[i].y - y_dim/2;
        cell_coords(i, 2) = cells[i].z - z_dim/2;
        cell_coords(i, 3) = cells[i].id;
        cell_coords(i, 4) = species[cells[i].id].genotype.size(); 
        cell_coords(i, 5) = sqrt(cell_coords(i,0)*cell_coords(i,0) + cell_coords(i,1)*cell_coords(i,1) + cell_coords(i,2)*cell_coords(i,2));
        if(species[cells[i].id].treatment_resistance) {
            cell_coords(i,6) = 1;
        } else {
            cell_coords(i,6) = 0;
        }
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

void PostProcessing::write_phylo_tree(std::vector<std::vector<int> > &phylo_tree, Rcpp::IntegerMatrix &rphylo_tree) {
    for(int i = 0; i < phylo_tree[0].size(); ++i) {
        rphylo_tree(i,0) = phylo_tree[0][i];
        rphylo_tree(i,1) = phylo_tree[1][i];
    }
}

Rcpp::NumericMatrix PostProcessing::get_color_scheme(std::vector<specie> &species) {
    Rcpp::NumericMatrix color_scheme(3, species.size());

    for(int i = 0; i < species.size(); ++i) {
        if(species[i].red > 1) {
            color_scheme(0,i) = 1;
        } else if(species[i].red < 0) {
            color_scheme(0,i) = 0;
        } else {
            color_scheme(0,i) = species[i].red;
        }
        if(species[i].green > 1) {
            color_scheme(1,i) = 1;
        } else if(species[i].green < 0) {
            color_scheme(1,i) = 0;
        } else {
            color_scheme(1,i) = species[i].green;
        }
        if(species[i].blue > 1) {
            color_scheme(2,i) = 1;
        } else if(species[i].blue < 0) {
            color_scheme(2,i) = 0;
        } else {
            color_scheme(2,i) = species[i].blue;
        }
    }
    return(color_scheme);
}