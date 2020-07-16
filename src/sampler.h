/*
SITH: An R package for visualizing and analyzing a spatial model of intra-tumor heterogeneity
Author: Phillip Nicol
License: GPL-2 
*/

/***********************************************************
random_index samples a random cell to be used in gillespie_step
Uniform random sampling can not be used because cells with different birth and death rates have
different probabilities of being selected.
Inverse transform sampling is O(n) in our setup, so this is infeasible as well. We want a sampler
that can be relied upon to be O(1). 
We use rejection sampling. Let p_max = max_i (b_i + d_i)
Given candidate cell (sampled uniformly random) and index j with sum of birth and death rate b_j + d_j, we generate
uniform rv u ~ U[0, p_max] and accept the cell if u < b_j + d_j and reject otherwise. 
This obtains the correct sampling probabilities, and for our purposes can be expected
to almost always accept a cell in very few iterations. 
*******************************************************/

#ifndef SAMPLER_H_INCLUDED
#define SAMPLER_H_INCLUDED

extern double p_max;

inline int selectIndexRS(std::vector<cell> &cells, std::vector<specie> &species)
{
    int index, trial;
    double u_trial;

    while(true)
    {
        trial = R::runif(0, cells.size());
        u_trial = R::runif(0, p_max);

        specie cell_species = species[cells[trial].id];
        
        if(u_trial < cell_species.b + cell_species.d)
        {
            index = trial;
            break;
        }
    }
    return index;
}

#endif 