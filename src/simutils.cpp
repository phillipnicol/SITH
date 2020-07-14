#include"simutils.h"

double p_max;
std::vector<int> drivers; 
int total_mutations;
int x_dim, y_dim, z_dim; //Size of the lattice (set at the start of simulation to accomodate num of cells)

bool*** lattice; 
std::vector<std::vector<int> > phylo_tree(2, std::vector<int>());
std::vector<std::vector<int> > perms;

void SimUtils::initIA(Rcpp::List input) {
    //Read input list (from R)
    std::vector<double> params = input["params"]; 
    int tumor_size = params[0]; 
    double wt_br = params[1]; 
    double wt_dr = params[2]; 
    double u = params[3]; 
    double du = params[4]; 
    double s = params[5]; 
    bool verbose = params[6];

    if(verbose) {Rcpp::Rcout << "Initializing structures ... ...\n";}

    //Init time 
    double time = 0; 

    //Initialize global variables
    gv_init(tumor_size, wt_br, wt_dr, u, du, s);

    //initialize empty lattice 
    lattice = init_lattice();

    std::vector<int> v;
    v.push_back(1); v.push_back(2); v.push_back(3); v.push_back(4); v.push_back(5); v.push_back(6);
    perms = get_perms(v);
}

cell SimUtils::initial_cell(std::vector<specie> &species, double wt_br, double wt_dr)
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

    //Save the species type in the vector
    species.push_back(initial_type);

    //Let the cell id point to this id
    cell.id = initial_type.id;
    return cell;
}

void SimUtils::trashcan(bool*** lattice)
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

bool*** init_lattice()
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

void gv_init(const int N, const double wt_br, const double wt_dr, const double u, const double du, const double s) {
    total_mutations = 0; 
    drivers.clear(); 
    p_max = wt_br + wt_dr;

    //error checking 
    if(N < 1) {Rcpp::stop("N must be at least 2.");}
    if(wt_dr > wt_br) {Rcpp::stop("Death rate can not be greater than birth rate.");}
    if((wt_br < 0) || (wt_dr < 0)) {Rcpp::stop("Birth and death rates must be non-negative.");}
    if(u < 0) {Rcpp::stop("u must be non-negative");}
    if(du < 0.0 || du > 1.0) {Rcpp::stop("du must be in [0,1]");}
    if(s < 0) {Rcpp::stop("s must be non-negative");}

    //set lattice dims 
    if(N > 100000000) {
        //very large lattice dims
        x_dim = 2000; y_dim = 2000; z_dim = 2000;
    }
    else if(N > 10000000) {
        //decently large
        x_dim = 1000; y_dim = 1000; z_dim = 1000;
    }
    else {
        //This should handle most cases
        x_dim = 500; y_dim = 500; z_dim = 500;
    }
}

//Store all permutations 
std::vector<std::vector<int> > get_perms(std::vector<int> v) {
    std::vector<std::vector<int> > perms;

    do {
        std::vector<int> p;
        for(int i = 0; i < v.size(); ++i) {
            p.push_back(v[i]);
        }
        perms.push_back(p);
    } while(std::next_permutation(v.begin(), v.end()));

    return perms;
}


std::vector<std::vector<Edge > > processG(Rcpp::NumericMatrix G) {
    std::vector<std::vector<Edge> > Gvec(G.nrow(), std::vector<Edge> ()); 

    for(int i = 0; i < G.nrow(); ++i) {
        Edge edge;
        edge.head = G[i,1];
        edge.u = G[i,2];
        edge.s = G[i,3]; 

        Gvec[i].push_back(edge); 
    }

    return Gvec; 
}