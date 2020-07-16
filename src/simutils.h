

#ifndef SIMUTILS_H
#define SIMUTILS_H

#include<iostream> 
#include<vector> 
#include<algorithm> 
#include<cmath> 
#include<float.h> 
#include<string> 
#include<Rcpp.h> 

//How often to print to screen
#define interval 2000000

//**** Structs for simulations ****// 
//Defined types for simulation
//A specie is a unique genotype in the cell population
struct specie {
    int id;
    int count;
    std::vector<int> genotype;
    double d,b;
};

//A cell is specified by its coordinates and specie type
struct cell {
    short int x,y,z;
    int id;
};

//For the multi-type branching process
struct Edge {
    int head; 
    double u, s; 
};

cell initial_cell(std::vector<specie> &species, double wt_br, double wt_dr);

bool*** init_lattice(void); 

void gv_init(const int N, const double wt_br, const double wt_dr, const double u, const double du, const double s);
std::vector<std::vector<int> > get_perms(std::vector<int> v);
std::vector<std::vector<Edge > > processG(Rcpp::NumericMatrix G);

namespace SimUtils {
    void initIA(Rcpp::List input);
    void initMTBP(Rcpp::List input);

    cell initial_cell(std::vector<specie> &species, double wt_br, double wt_dr);    

    void trashcan(bool*** lattice); 
}

inline int max_mut(std::vector<specie> &species) {
    int max = 0;
    for(int i = 0; i < species.size(); ++i) {
        if(species[i].genotype.size() > max) {
            max = species[i].genotype.size(); 
        }
    }
    return max;
}

//obviously C++17 has std::sort, which implements quick or merge sort.
//However, I think it is dangerous to use C++17 functionality since the package
//is tested on a wide range of operating systems... so let's just write bubblesort.
//Probably no difference between bubble and merge for our purposes. 
inline std::vector<int> bubblesort(std::vector<int> a) {
    for(int i = 0; i < a.size(); ++i) {
        for(int j = i+1; j < a.size(); ++j) {
            if(a[j] < a[i]) {
                std::swap(a[i],a[j]);
            }
        }
    }
    return a;
}

#endif 