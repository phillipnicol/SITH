// preprocessing
#include<iostream>
#include<vector>
#include<random>
#include<algorithm>
#include<cmath>
#include<fstream>
#include<set>
#include<float.h>
#include<string>
#include<Rcpp.h>


//define dimensions of lattice
//memory required: (x_dim * y_dim * z_dim) bits
#define x_dim 500
#define y_dim 500
#define z_dim 500

//take data at different time points (1 = true)
#define time_course 0
#define interval 500000

#define rgb_ub 0.91
#define rgb_lb 0.09

//Initialize prng
std::random_device rd;
std::mt19937 generator(rd());

//Make static distributions 
std::poisson_distribution<int> mut(1.0);
std::bernoulli_distribution dmut(1.0);

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
    int x,y,z;
    specie species;
};

//Global variables
double p_max;
std::vector<int> drivers; 
int total_mutations;

bool*** init_lattice(void)
{
    bool*** lattice = new bool**[x_dim];
    for(int i = 0; i < x_dim; ++i)
    {
        lattice[i] = new bool*[y_dim];
        for(int j = 0; j < y_dim; ++j)
        {
            lattice[i][j] = new bool[z_dim];
        }
    }
    //Initialize all zeros
    for(int i = 0; i < x_dim; ++i)
    {
        for(int j = 0; j < y_dim; ++j)
        {
            for(int k = 0; k < z_dim; ++k)
            {
                lattice[i][j][k] = 0;
            }
        }
    }
    lattice[x_dim/2][y_dim/2][z_dim/2] = 1;
    return lattice;
}

void trashcan(bool*** lattice)
{
    // deallocate
    for(int i = 0; i < x_dim; ++i)
    {
        for(int j = 0; j < y_dim; ++j)
        {
            delete[] lattice[i][j];
        } 
        delete[] lattice[i];
    }
    delete[] lattice;
}

cell initial_cell(std::vector<specie> &species, double wt_br, double wt_dr)
{
    //Initial cell lies in the center of the lattice
    cell cell;
    cell.x = x_dim/2;
    cell.y = y_dim/2;
    cell.z = z_dim/2;

    //initial specie type has wild type birth and death rates
    specie initial_type;
    initial_type.b = wt_br;
    initial_type.d = wt_dr;
    initial_type.id = 0;
    initial_type.count = 1;
    initial_type.genotype.push_back(0);

    species.push_back(initial_type);

    //save the species type in vector
    cell.species = initial_type;
    return cell;
}

int max_mut(std::vector<specie> &species) {
    int max = 0;
    for(int i = 0; i < species.size(); ++i) {
        if(species[i].genotype.size() > max) {
            max = species[i].genotype.size(); 
        }
    }
    return max;
}