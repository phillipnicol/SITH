
double jaccard_similarity(std::vector<int> a, std::vector<int> b)
{
    std::set<int> set_union;
    int intersection = 0;

    int ita = 0;
    int itb = 0;
    while(ita < a.size() && itb < b.size())
    {
        set_union.insert(a[ita]);
        set_union.insert(a[itb]);
        if(a[ita] == b[itb])
        {
            ++intersection;
            ++ita;
            ++itb;
        }
        else if(a[ita] < b[itb])
        {
            ++ita;
        }
        else
        {
            ++itb;
        }
    }

    return (double)intersection/(double)set_union.size();
}

double avg_mutants(std::vector<cell> &cells)
{
    int n_mutants = 0;
    for(int i = 0; i < cells.size(); ++i)
    {
        n_mutants += cells[i].species.genotype.size();
    }
    double avg_mutants = (double)n_mutants/(double)cells.size();
    return avg_mutants;
}

double avg_jaccard(std::vector<cell> &cells, int N)
{
    double sum = 0.0;
    std::uniform_int_distribution<int> d(0, cells.size() - 1);
    for(int i = 0; i < N; ++i)
    {
        int xi = d(generator);
        int xj = d(generator);
        sum += jaccard_similarity(cells[xi].species.genotype, cells[xj].species.genotype);
    }
    return sum/(double)N;
}

void time_course_data(std::vector<cell> &cells, double time, std::ofstream &avg_mutants_time, std::ofstream &avg_jaccard_time)
{
    avg_mutants_time << time << ", " << avg_mutants(cells) << std::endl;
    avg_jaccard_time << time << ", " << avg_jaccard(cells, 5000) << std::endl;
}