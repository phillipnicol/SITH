
cell birth_cell(cell cell, int key, std::vector<specie> &species, double wt_dr, double u, double du, double multiplicative_udpate)
{
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

    static std::bernoulli_distribution dist(u);
    static std::bernoulli_distribution dist2(du/u);
    if(dist(generator) == 1)
    {
        //Mutation has occurred 
        int mutation_count = species.size();

        specie new_species;
        new_species.id = mutation_count;
        //du gives unconditional probability of driver, so du/u gives probability of driver
        //Assume du < u on input, otherwise this will throw an error
        if(dist2(generator) == 1)
        {
            //apply multiplicative update
            new_species.b = multiplicative_udpate*cell.species.b;
            ++drivers;
            if(new_species.b + wt_dr > p_max)
            {
                p_max = new_species.b + wt_dr;
            }
        }
        else
        {
            //If no driver, birth rate remains the same
            new_species.b = cell.species.b;
        }
        //Death rate stays the same in mutation
        new_species.d = wt_dr;
        new_species.count = 1;
        new_species.genotype = cell.species.genotype;
        new_species.genotype.push_back(mutation_count);

        //Save the new species
        species.push_back(new_species);

        new_cell.species = new_species;
    }
    else
    {
        //If no mutation, the cell inherits all of the qualities of its parent
        new_cell.species.b = cell.species.b;
        new_cell.species.d = cell.species.d;
        new_cell.species.id = cell.species.id;
        new_cell.species.genotype = cell.species.genotype;
        ++species[cell.species.id].count;
    }
    return new_cell;
}


void gillespie_step(std::vector<cell> &cells, std::vector<specie> &species, int index, bool*** lattice, double &time,
              double wt_dr, double u, double du, double multiplicative_update)
{
    //Randomly selected cell (from previous step)
    cell cell = cells[index];
    //Look at the neighbors of the cell and find a (random) neighbor
    //If no random neighbor, key = 0
    int key = random_neighbor(cell, lattice);
    if(key != 0)
    {
        //key != 0 so there is at least one free neighbor
        //Probability of birth event
        std::bernoulli_distribution dist(cell.species.b/(cell.species.b + cell.species.d));
        if(dist(generator) == 1)
        {
            //Birth
            update_lattice(cell, key, lattice);
            struct cell new_cell = birth_cell(cell, key, species, wt_dr, u, du, multiplicative_update);
            cells.push_back(new_cell);
            //Update time (approximate)
            double lambda = cells.size()*p_max;
            std::exponential_distribution<double> dist(lambda);
            time += dist(generator);            
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
                double lambda = cells.size()*p_max;
                std::exponential_distribution<double> dist(lambda);
                time += dist(generator);             
            }
        }
    }
    else
    {
        //Key = 0 and the cell has no free neighbors
        //No birth can occur, but the cell could still die
        std::bernoulli_distribution dist(cell.species.b/(cell.species.b + cell.species.d));
        if(dist(generator) == 0)
        {
            //Death
            if(cells.size() > 1)
            {
                //procedure same as above
                lattice[cell.x][cell.y][cell.z] = 0;
                std::swap(cells[index], cells.back());
                cells.pop_back();
                --species[cell.species.id].count;
                double lambda = cells.size()*p_max;
                std::exponential_distribution<double> dist(lambda);
                time += dist(generator);    
            }        
        }
    }
}
