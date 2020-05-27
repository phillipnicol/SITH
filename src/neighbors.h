
//free_neighbor returns true if there is an open space at the neighbor position specified by key
bool free_neighbor(cell cell, bool*** lattice, int key)
{
    if(key == 1)
    {
        if(cell.x < x_dim - 1)
        {
            if(lattice[cell.x+1][cell.y][cell.z] == 0)
            {
                return true;
            }
        }
    }
    else if(key == 2)
    {
        if(cell.x > 0)
        {
            if(lattice[cell.x-1][cell.y][cell.z] == 0)
            {
                return true;
            }
        }
    }
    else if(key == 3)
    {
        if(cell.y < y_dim - 1)
        {
            if(lattice[cell.x][cell.y+1][cell.z] == 0)
            {
                return true;
            }
        }
    }
    else if(key == 4)
    {
        if(cell.y > 0)
        {
            if(lattice[cell.x][cell.y-1][cell.z] == 0)
            {
                return true;
            }
        }
    }
    else if(key == 5)
    {
        if(cell.z < z_dim - 1)
        {
            if(lattice[cell.x][cell.y][cell.z+1] == 0)
            {
                return true;
            }
        }
    }
    else
    {
        if(cell.z > 0)
        {
            if(lattice[cell.x][cell.y][cell.z-1] == 0)
            {
                return true;
            }
        }
    }
    return false;
}

int random_neighbor(cell cell, bool*** lattice)
{
    //randomly permute an array of keys
    static int nbhd[6] = {1,2,3,4,5,6};
    std::shuffle(nbhd, nbhd+6, generator);
    for(int i = 0; i < 6; ++i)
    {
        if(free_neighbor(cell, lattice, nbhd[i]) == true)
        {
            //if a free neighbor is found, this is the key
            return nbhd[i];
        }
    }
    //otherwise, no free neighbors, and key is 0
    return 0;
}

//update_lattice puts a 1 in the location where a new cell is born
void update_lattice(cell cell, int key, bool*** lattice)
{
    if(key == 1)
    {
        lattice[cell.x+1][cell.y][cell.z] = 1;
    }
    else if(key == 2)
    {
        lattice[cell.x-1][cell.y][cell.z] = 1;
    }
    else if(key == 3)
    {
        lattice[cell.x][cell.y+1][cell.z] = 1;
    }
    else if(key == 4)
    {
        lattice[cell.x][cell.y-1][cell.z] = 1;
    }
    else if(key == 5)
    {
        lattice[cell.x][cell.y][cell.z+1] = 1;
    }
    else
    {
        lattice[cell.x][cell.y][cell.z-1] = 1;
    }
}
