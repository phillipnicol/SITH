/*
SITH: An R package for visualizing and analyzing a spatial model of intra-tumor heterogeneity
Author: Phillip Nicol
License: GPL-2 
*/

#ifndef GILLESPIE_H_INCLUDED
#define GILLESPIE_H_INCLUDED

#include"neighbors.h"

extern double p_max;
extern std::vector<int> drivers; 
extern int total_mutations;
extern int x_dim, y_dim, z_dim; //Size of the lattice (set at the start of simulation to accomodate num of cells)

extern bool*** lattice; 
extern std::vector<std::vector<int> > phylo_tree;
extern std::vector<std::vector<Edge> > G; 

cell birth_cellIA(cell &cell, const int key, const specie cell_species, std::vector<specie> &species, 
                const double wt_dr, const double u, const double du, const double s);
cell birth_cellUDT(cell &cell, const int key, specie cell_species, std::vector<specie> &species);

int find_gtype(std::vector<specie> &species, std::vector<int> gtype);
bool vin(std::vector<int> v, int a);

namespace Gillespie {
    void gillespieIA(std::vector<cell> &cells, std::vector<specie> &species, const int index, double &time,
                double &lambda, const double wt_br, const double wt_dr, const double u, const double du, const double s);
    void gillespieUDT(std::vector<cell> &cells, std::vector<specie> &species, const int index, double &time);
}


#endif 
