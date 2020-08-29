#include"gillespie.h"

void Gillespie::gillespieIA(std::vector<cell> &cells, std::vector<specie> &species, const int index, double &time,
                const double wt_dr, const double u, const double du, const double s) {

    //Update time--approximate 
    double lambda = 1/(cells.size()*p_max);
    time += R::rexp(lambda);

    //Randomly selected cell (from previous step)
    cell cell = cells[index];

    //Find the allele of the chosen cell 
    specie cell_species = species[cell.id];

    //Look at the neighbors of the cell and find a (random) neighbor
    //If no random neighbor, key = 0
    int key = random_neighbor(cell);
    if(key != 0)
    {
        //key != 0 so there is at least one free neighbor
        //Probability of birth event
        int bd = R::rbinom(1,cell_species.b/(cell_species.b + cell_species.d));
        if(bd == 1) {
            //Birth
            update_lattice(cell, key, lattice);
            struct cell new_cell = birth_cellIA(cells[index], key, cell_species, species, wt_dr, u, du, s);
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

void Gillespie::gillespieUDT(std::vector<cell> &cells, std::vector<specie> &species, const int index, double &time) {

    //Update time--approximate 
    double lambda = 1/(cells.size()*p_max);
    time += R::rexp(lambda);

    //Randomly selected cell (from previous step)
    cell cell = cells[index];

    //Find the type of the chosen cell 
    specie cell_species = species[cell.id];

    //Look at the neighbors of the cell and find a (random) neighbor
    //If no random neighbor, key = 0
    int key = random_neighbor(cell);
    if(key != 0)
    {
        //key != 0 so there is at least one free neighbor
        //Probability of birth event
        int bd = R::rbinom(1,cell_species.b/(cell_species.b + cell_species.d));
        if(bd == 1) {
            //Birth
            update_lattice(cell, key, lattice);
            struct cell new_cell = birth_cellUDT(cells[index], key, cell_species, species);
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

cell birth_cellIA(cell &cell, const int key, const specie cell_species, std::vector<specie> &species, 
                            const double wt_dr, const double u, const double du, const double s) {
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

cell birth_cellUDT(cell &cell, const int key, specie cell_species, std::vector<specie> &species) {
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

    bool mutflag = false; 

    for(std::vector<int>::iterator it = cell_species.genotype.begin(); it != cell_species.genotype.end(); ++it) {
        for(int j = 0; j < G[*it].size(); ++j) {
            if((int)R::rbinom(1,G[*it][j].u)) {
                //Mutation event 
                if((int)R::rbinom(1,0.5)) {
                    //New daughter cell mutates
                    std::vector<int> gtype = cell_species.genotype;
                    gtype = bubblesort(gtype);
                    if(vin(gtype, G[*it][j].head)) {
                        continue; 
                    }
                    mutflag = true; 

                    gtype.push_back(G[*it][j].head);
                    gtype = bubblesort(gtype); 

                    int sp_id = find_gtype(species, gtype);
                    if(sp_id == -1) {
                        //New genotype
                        specie new_specie;
                        new_specie.b = species[cell.id].b*G[*it][j].s;
                        new_specie.d = species[cell.id].d;
                        if(new_specie.b + new_specie.d > p_max) {p_max = new_specie.b + new_specie.d;}
                        new_specie.id = species.size(); 
                        new_specie.genotype = gtype;
                        new_specie.count = 1; 

                        new_cell.id = new_specie.id; 

                        species.push_back(new_specie); 
                    }
                    else {
                        new_cell.id = sp_id; 
                        ++species[sp_id].count; 
                    }
                }
                else {
                    //Original cell mutates 
                    std::vector<int> gtype = cell_species.genotype;
                    gtype = bubblesort(gtype);

                    if(vin(gtype, G[*it][j].head)) {
                        continue; 
                    }
                    mutflag = true; 

                    //the daughter gets what the original cell had
                    new_cell.id = cell.id;
                    //no need to update count since we will subtract it from original cell 

                    gtype.push_back(G[*it][j].head);
                    gtype = bubblesort(gtype); 

                    int sp_id = find_gtype(species, gtype);
                    if(sp_id == -1) {
                        //New genotype
                        specie new_specie;
                        new_specie.b = species[cell.id].b*G[*it][j].s;
                        new_specie.d = species[cell.id].d;
                        if(new_specie.b + new_specie.d > p_max) {p_max = new_specie.b + new_specie.d;}
                        new_specie.id = species.size(); 
                        new_specie.genotype = gtype;
                        new_specie.count = 1; 

                        cell.id = new_specie.id; 

                        species.push_back(new_specie); 
                    }
                    else {
                        cell.id = sp_id; 
                        ++species[cell.id].count; 
                    }
                }
            }
            if(mutflag) {
                return new_cell; 
            }
        }
    }
    //No mutation occurred
    new_cell.id = cell.id;
    ++species[cell.id].count;

    return new_cell;
}

int find_gtype(std::vector<specie> &species, std::vector<int> gtype) {
    for(int i = 0; i < species.size(); ++i) {
        if(species[i].genotype == gtype) {
            return i;
        }
    }
    return -1; 
}

bool vin(std::vector<int> v, int a) {
    for(int i = 0; i < v.size(); ++i) {
        if(v[i] == a) {
            return true;
        }
    }
    return false; 
}