#ifndef SIMULATIONS_H
#define SIMULATIONS_H

#include"simutils.h"
#include"sampler.h"
#include"gillespie.h"
#include"postproc.h"

extern std::vector<int> drivers; 

namespace Sims {
    Rcpp::List simulateIA(Rcpp::List input); 
    Rcpp::List simulateMTBP(Rcpp::List input);
}







#endif 