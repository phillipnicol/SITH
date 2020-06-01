
cell birth_cell(cell &cell, const int key, std::vector<specie> &species, const double wt_dr, const double u, const double du, 
                const double multiplicative_udpate, std::vector<std::vector<int> > &phylo_tree) {
    struct cell new_cell;
    new_cell.x = cell.x;
    new_cell.y = cell.y;
    new_cell.z = cell.z;

    //based on key value, coordinates of new cell are chosen
    if(key == 1)
    {
        ++new_cell.x;
    }
    else if(key == 2)
    {
        --new_cell.x;
    }
    else if(key == 3)
    {
        ++new_cell.y;
    }
    else if(key == 4)
    {
        --new_cell.y;
    }
    else if(key == 5)
    {
        ++new_cell.z;
    }
    else
    {
        --new_cell.z;
    } 

    //daughter cell 
    int nmuts = R::rpois(u);
    if(nmuts > 0) {
        specie new_species;
        new_species.id = species.size(); 
        std::vector<int> new_gt = cell.species.genotype;
        double br = cell.species.b;
        for(int i = 0; i < nmuts; ++i) {
            ++total_mutations;
            new_gt.push_back(total_mutations); 

            phylo_tree[0].push_back(cell.species.genotype.back());
            phylo_tree[1].push_back(total_mutations); 

            if((int)R::rnbinom(1,du) == 1) {
                //apply multiplicative update
                br *= multiplicative_udpate;
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
        new_cell.species = new_species;
            
    }
    else {
        //If no mutation, the cell inherits all of the qualities of its parent
        new_cell.species.b = cell.species.b;
        new_cell.species.d = cell.species.d;
        new_cell.species.id = cell.species.id;
        new_cell.species.genotype = cell.species.genotype;
        ++species[cell.species.id].count;
    }

    //original cell may mutate as well
    nmuts = R::rpois(u);
    if(nmuts > 0) {
        specie new_species;
        new_species.id = species.size(); 
        std::vector<int> new_gt = cell.species.genotype;
            
        double br = cell.species.b;
        for(int i = 0; i < nmuts; ++i) {
            ++total_mutations;
            new_gt.push_back(total_mutations); 

            phylo_tree[0].push_back(cell.species.genotype.back());
            phylo_tree[1].push_back(total_mutations); 

            if((int)R::rnbinom(1,du) == 1) {
                //apply multiplicative update
                br *= multiplicative_udpate;
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

        species[cell.species.id].count--;

        species.push_back(new_species);
        cell.species = new_species;
    }

    return new_cell;
}


void gillespie_step(std::vector<cell> &cells, std::vector<specie> &species, const int index, bool*** lattice, double &time,
              const double wt_dr, const double u, const double du, const double multiplicative_update, std::vector<std::vector<int> > &phylo_tree) {

    //Randomly selected cell (from previous step)
    cell cell = cells[index];
    //Look at the neighbors of the cell and find a (random) neighbor
    //If no random neighbor, key = 0
    int key = random_neighbor(cell, lattice);
    if(key != 0)
    {
        //key != 0 so there is at least one free neighbor
        //Probability of birth event
        int bd = R::rbinom(1,cell.species.b/(cell.species.b + cell.species.d));
        //std::bernoulli_distribution dist(cell.species.b/(cell.species.b + cell.species.d));
        if(bd == 1)
        {
            //Birth
            update_lattice(cell, key, lattice);
            struct cell new_cell = birth_cell(cells[index], key, species, wt_dr, u, du, multiplicative_update, phylo_tree);
            cells.push_back(new_cell);
            //Update time (approximate)
            double lambda = 1/(cells.size()*p_max);
            time += R::rexp(lambda);            
        }
        else
        {
            //Death
            if(cells.size() > 1)
            {
                //Free up the space in the lattice
                lattice[cell.x][cell.y][cell.z] = 0;
                //remove cell from list
                std::swap(cells[index], cells.back());
                cells.pop_back();   
                --species[cell.species.id].count;
                //update time(approximate)
                double lambda = 1/(cells.size()*p_max);
                time += R::rexp(lambda);           
            }
        }
    }
    else
    {
        //Key = 0 and the cell has no free neighbors
        //No birth can occur, but the cell could still die
        int bd = R::rbinom(1,cell.species.b/(cell.species.b + cell.species.d));
        //std::bernoulli_distribution dist(cell.species.b/(cell.species.b + cell.species.d));
        if(bd == 0)
        {
            //Death
            if(cells.size() > 1)
            {
                //procedure same as above
                lattice[cell.x][cell.y][cell.z] = 0;
                std::swap(cells[index], cells.back());
                cells.pop_back();
                --species[cell.species.id].count;
                double lambda = 1/(cells.size()*p_max);
                time += R::rexp(lambda);
            }        
        }
    }
}