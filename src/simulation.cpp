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

    double time = 0;
    p_max = wt_br + wt_dr;

    std::cout << "Initializing structures..." << std::endl;
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
        gillespie_step(cells, species, index, lattice, time, wt_dr, u, du, multiplicative_update);
        ++iteration;
        if(iteration % interval == 0)
        {
            Rcpp::Rcout << "Simulated time: " << time << " Population: " << cells.size() << std::endl;
            iteration = 1;
        }
    }

    //Print summary of simulation
    std::cout << "Process Complete... Releasing memory..." << std::endl;
    std::cout << "There are " << drivers << " driver mutations" << std::endl;
    std::cout << "Simulated time is " << time << std::endl;
    end = clock();
    std::cout << "Simulation completed in " << (double)(end - start)/CLOCKS_PER_SEC << " s." << std::endl;
    trashcan(lattice);
    

    //save the data
    Rcpp::NumericMatrix cell_coords(cells.size(), 6);

    int maximum_mut = max_mut(species);
    Rcpp::IntegerMatrix species_dict(species.size(), maximum_mut+1); 

    Rcpp::IntegerVector muts(species.size());
    write_results(cells, species, cell_coords, species_dict, muts);

    //create list 
    Rcpp::List out = Rcpp::List::create();
    out.push_back(cell_coords);
    out.push_back(species_dict);
    out.push_back(muts);
    out.push_back(species.size());
    return(out);
}