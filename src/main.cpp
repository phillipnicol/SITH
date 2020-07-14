#include"simulations.h"

// [[Rcpp::export]]
Rcpp::List simulate_tumor(Rcpp::List input) {
    SimUtils::initIA(input); 
    Rcpp::List out = Sims::simulateIA(input);
    return out; 
}

