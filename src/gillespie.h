/*
filename: gillespie.h
project: TumorGenerator R package
author: Phillip B. Nicol
license: GPL v3 
date: June 4, 2020

summary: contains functions that perform cell division, cell death, 
mutation, and updates time. 
*/

cell birth_cell(cell &cell, const int key, const specie cell_species, std::vector<specie> &species, const double wt_dr, const double u, const double du, 
                const double s, std::vector<std::vector<int> > &phylo_tree) {
    //form new cell 
    struct cell new_cell;
    new_cell.x = cell.x;
    new_cell.y = cell.y;
    new_cell.z = cell.z;

    //Check for which coordinate to update 
    switch(key) {
        case 1: ++new_cell.x; break;
        case 2: --new_cell.x; break; 
        case 3: ++new_cell.y; break; 
        case 4: --new_cell.y; break; 
        case 5: ++new_cell.z; break; 
        case 6: --new_cell.z; break; 
        default: Rcpp::stop("Algorithm internal error.");
    }

    //daughter cell 
    //Receieve a poisson number of genetic alterations 
    int nmuts = R::rpois(u);
    if(nmuts > 0) {
        specie new_species;
        new_species.id = species.size(); 
        std::vector<int> new_gt = cell_species.genotype;
        double br = cell_species.b;
        for(int i = 0; i < nmuts; ++i) {
            ++total_mutations;
            new_gt.push_back(total_mutations); 

            phylo_tree[0].push_back(cell_species.genotype.back());
            phylo_tree[1].push_back(total_mutations); 

            if((int)R::rbinom(1,du) == 1) {
                //apply multiplicative update
                br *= s;
                drivers.push_back(total_mutations);          
            }
        }
        new_species.b = br;
        new_species.d = wt_dr;
        new_species.genotype = new_gt;
        new_species.count = 1;

        if(new_species.b + wt_dr > p_max)
        {
            p_max = new_species.b + wt_dr;
        }

        species.push_back(new_species);
        new_cell.id = new_species.id;
            
    }
    else {
        //If no mutation, the cell inherits all of the qualities of its parent
        new_cell.id = cell_species.id; 
        ++species[cell.id].count;
    }

    //original cell may mutate as well
    nmuts = R::rpois(u);
    if(nmuts > 0) {
        specie new_species;
        new_species.id = species.size(); 
        std::vector<int> new_gt = cell_species.genotype;
            
        double br = cell_species.b;
        for(int i = 0; i < nmuts; ++i) {
            ++total_mutations;
            new_gt.push_back(total_mutations); 

            phylo_tree[0].push_back(cell_species.genotype.back());
            phylo_tree[1].push_back(total_mutations); 

            if((int)R::rbinom(1,du) == 1) {
                //apply multiplicative update
                br *= s;
                drivers.push_back(total_mutations);          
            }
        }
        new_species.b = br;
        new_species.d = wt_dr;
        new_species.genotype = new_gt;
        new_species.count = 1;

        if(new_species.b + wt_dr > p_max)
        {
            p_max = new_species.b + wt_dr;
        }

        species[cell_species.id].count--;

        species.push_back(new_species);
        cell.id = new_species.id;
    }

    return new_cell;
}

void gillespie_step(std::vector<cell> &cells, std::vector<specie> &species, const int index, bool*** lattice, double &time,
              const double wt_dr, const double u, const double du, const double s, std::vector<std::vector<int> > &phylo_tree,
              std::vector<std::vector<int> > &perms) {

    //Update time--approximate 
    double lambda = 1/(cells.size()*p_max);
    time += R::rexp(lambda);

    //Randomly selected cell (from previous step)
    cell cell = cells[index];

    //Find the allele of the chosen cell 
    specie cell_species = species[cell.id];

    //Look at the neighbors of the cell and find a (random) neighbor
    //If no random neighbor, key = 0
    int key = random_neighbor(cell, lattice, perms);
    if(key != 0)
    {
        //key != 0 so there is at least one free neighbor
        //Probability of birth event
        int bd = R::rbinom(1,cell_species.b/(cell_species.b + cell_species.d));
        if(bd == 1) {
            //Birth
            update_lattice(cell, key, lattice);
            struct cell new_cell = birth_cell(cells[index], key, cell_species, species, wt_dr, u, du, s, phylo_tree);
            cells.push_back(new_cell);         
        }
        else
        {
            //Death
            if(cells.size() > 1) {
                //Free up the space in the lattice
                lattice[cell.x][cell.y][cell.z] = 0;
                //remove cell from list
                std::swap(cells[index], cells.back());
                cells.pop_back();   
                --species[cell.id].count;    
            }
        }
    }
    else
    {
        //Key = 0 and the cell has no free neighbors
        //No birth can occur, but the cell could still die
        int bd = R::rbinom(1,cell_species.b/(cell_species.b + cell_species.d));
        if(bd == 0)
        {
            //Death
            if(cells.size() > 1)
            {
                //procedure same as above
                lattice[cell.x][cell.y][cell.z] = 0;
                std::swap(cells[index], cells.back());
                cells.pop_back();
                --species[cell.id].count;
            }        
        }
    }
}
