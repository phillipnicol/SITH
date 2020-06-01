#include"setup.h"
#include"neighbors.h"
#include"sampler.h"
#include"gillespie.h"
#include"time_course.h"
#include"post_processing.h"

// [[Rcpp::export]]
Rcpp::List simulate_tumor(Rcpp::List input) {
    clock_t start, end;
    start = clock();

    std::vector<double> params = input["params"];
    int tumor_size = params[0];
    double wt_br = params[1];
    double wt_dr = params[2];
    double u = params[3];
    double du = params[4];
    double multiplicative_update = params[5];
    bool verbose = params[6];

    double time = 0;
    p_max = wt_br + wt_dr;


    //INIT global vars
    total_mutations = 0;
    drivers.clear(); 
    p_max = 1.0;
    mut = std::poisson_distribution<int>(u);
    dmut = std::bernoulli_distribution(du);

    std::vector<std::vector<int> > phylo_tree(2, std::vector<int>());

    if(verbose) {Rcpp::Rcout << "Initializing structures ... ...\n";}
    //initialize empty lattice
    bool*** lattice = init_lattice();

    //initialize vector of cells:
    std::vector<cell> cells;
    std::vector<specie> species;
    cells.push_back(initial_cell(species, wt_br, wt_dr));

    int index;
    int iteration = 1;
    while(cells.size() < tumor_size)
    {
        index = random_index(cells, species);
        gillespie_step(cells, species, index, lattice, time, wt_dr, u, du, multiplicative_update, phylo_tree);
        ++iteration;
        if(iteration % interval == 0)
        {
            if(verbose) {Rcpp::Rcout << "Simulated time: " << time << " days. Population is " << cells.size() << " cells. \n";}
            iteration = 1;
        }
    }

    //Print summary of simulation
    if(verbose) {Rcpp::Rcout << "Simulation complete. Releasing memory ... ... \n";}
    trashcan(lattice);   

    if(verbose) {Rcpp::Rcout << "Simulated time is " << time << "\n";}

    end = clock();
    if(verbose) {Rcpp::Rcout << "Simulation completed in " << (double)(end - start)/CLOCKS_PER_SEC << " s.\n";}

    //save the data
    Rcpp::NumericMatrix cell_coords(cells.size(), 6);

    int maximum_mut = max_mut(species);
    Rcpp::IntegerMatrix species_dict(species.size(), maximum_mut+1); 

    Rcpp::IntegerVector muts(total_mutations+1);
    write_results(cells, species, cell_coords, species_dict, muts);

    Rcpp::IntegerVector driver_muts = Rcpp::wrap(drivers);

    Rcpp::IntegerMatrix rphylo_tree(phylo_tree[0].size(), 2);
    write_phylo_tree(phylo_tree, rphylo_tree);

    Rcpp::NumericMatrix color_scheme = get_color_scheme(species);

    //create list 
    Rcpp::List out = Rcpp::List::create();
    out.push_back(cell_coords);
    out.push_back(species_dict);
    out.push_back(muts);
    out.push_back(rphylo_tree);
    out.push_back(color_scheme);
    out.push_back(species.size());
    out.push_back(driver_muts);
    return(out);
}