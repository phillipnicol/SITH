#include"simulations.h"

Rcpp::List Sims::simulateIA(Rcpp::List input) {
    //Read input list (from R)
    std::vector<double> params = input["params"]; 
    int tumor_size = params[0]; 
    double wt_br = params[1]; 
    double wt_dr = params[2]; 
    double u = params[3]; 
    double du = params[4]; 
    double s = params[5]; 
    bool verbose = params[6];

    //Init time 
    double time = 0; 

    //initialize vector of cells:
    std::vector<cell> cells; 
    std::vector<specie> species;
    cells.push_back(SimUtils::initial_cell(species, wt_br, wt_dr));

    int index;
    int iteration = 1;

    //start clock
    clock_t start, end; 
    start = clock();

    //main simulation loop
    while(cells.size() < tumor_size)
    {
        index = selectIndexRS(cells, species);
        Gillespie::gillespieIA(cells, species, index, time, wt_dr, u, du, s);
        ++iteration;
        if(iteration % interval == 0)
        {
            if(verbose) {Rcpp::Rcout << "Simulated time: " << time << " days. Population is " << cells.size() << " cells. \n";}
            iteration = 1;
        }
    }

    //Print summary of simulation
    if(verbose) {Rcpp::Rcout << "Simulation complete. Releasing memory ... ... \n";}
    SimUtils::trashcan(lattice);   

    if(verbose) {Rcpp::Rcout << "Simulated time is " << time << " days \n";}

    end = clock();
    if(verbose) {Rcpp::Rcout << "Simulation completed in " << (double)(end - start)/CLOCKS_PER_SEC << " s.\n";}
    if(verbose) {Rcpp::Rcout << "Writing results ... ... \n";}

    //save the data and write them to R objects 
    Rcpp::NumericMatrix cell_coords(cells.size(), 6);

    int maximum_mut = max_mut(species);
    Rcpp::IntegerMatrix species_dict(species.size(), maximum_mut+1); 

    Rcpp::IntegerVector muts(total_mutations+1);
    PostProcessing::write_results(cells, species, cell_coords, species_dict, muts);

    Rcpp::IntegerVector driver_muts = Rcpp::wrap(drivers);

    Rcpp::IntegerMatrix rphylo_tree(phylo_tree[0].size(), 2);
    PostProcessing::write_phylo_tree(phylo_tree, rphylo_tree);

    Rcpp::NumericMatrix color_scheme = PostProcessing::get_color_scheme(species);

    //create list 
    Rcpp::List out = Rcpp::List::create();
    out.push_back(cell_coords);
    out.push_back(species_dict);
    out.push_back(muts);
    out.push_back(rphylo_tree);
    out.push_back(color_scheme);
    out.push_back(species.size());
    out.push_back(driver_muts);
    out.push_back(time);
    return(out);
}

Rcpp::List Sims::simulateMTBP(Rcpp::List input) {
    //Read input list (from R)
    std::vector<double> params = input["params"]; 
    int tumor_size = params[0]; 
    double wt_br = params[1]; 
    double wt_dr = params[2]; 
    bool verbose = params[3];

    //Init time 
    double time = 0; 

    //initialize vector of cells:
    std::vector<cell> cells; 
    std::vector<specie> species;
    cells.push_back(SimUtils::initial_cell(species, wt_br, wt_dr));

    int index;
    int iteration = 1;

    //start clock
    clock_t start, end; 
    start = clock();

    //main simulation loop
    while(cells.size() < tumor_size)
    {
        index = selectIndexRS(cells, species);
        Gillespie::gillespieMTBP(cells, species, index, time);
        ++iteration;
        if(iteration % interval == 0)
        {
            if(verbose) {Rcpp::Rcout << "Simulated time: " << time << " days. Population is " << cells.size() << " cells. \n";}
            iteration = 1;
        }
    }

    for(auto sp : species) {
        for(int j = 0; j < sp.genotype.size(); ++j) {
            std::cout << sp.genotype[j] << " ";
        }
        std::cout << std::endl;
    }

    //Print summary of simulation
    if(verbose) {Rcpp::Rcout << "Simulation complete. Releasing memory ... ... \n";}
    SimUtils::trashcan(lattice);   

    if(verbose) {Rcpp::Rcout << "Simulated time is " << time << " days \n";}

    end = clock();
    if(verbose) {Rcpp::Rcout << "Simulation completed in " << (double)(end - start)/CLOCKS_PER_SEC << " s.\n";}
    if(verbose) {Rcpp::Rcout << "Writing results ... ... \n";}

    //save the data and write them to R objects 
    Rcpp::NumericMatrix cell_coords(cells.size(), 6);

    int maximum_mut = max_mut(species);
    Rcpp::IntegerMatrix species_dict(species.size(), maximum_mut+1); 

    Rcpp::IntegerVector muts(maximum_mut+1);
    PostProcessing::write_results(cells, species, cell_coords, species_dict, muts);

    Rcpp::IntegerVector driver_muts = Rcpp::wrap(drivers);

    Rcpp::IntegerMatrix rphylo_tree(1,2);

    Rcpp::NumericMatrix color_scheme = PostProcessing::get_color_scheme(species);

    //create list 
    Rcpp::List out = Rcpp::List::create();
    out.push_back(cell_coords);
    out.push_back(species_dict);
    out.push_back(muts);
    out.push_back(rphylo_tree);
    out.push_back(color_scheme);
    out.push_back(species.size());
    out.push_back(driver_muts);
    out.push_back(time);
    return(out);    
}