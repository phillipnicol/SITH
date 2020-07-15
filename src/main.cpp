#include"simulations.h"

// [[Rcpp::export]]
Rcpp::List simulate_tumor(Rcpp::List input) {
    SimUtils::initIA(input); 
    Rcpp::List out = Sims::simulateIA(input);
    return out; 
}

// [[Rcpp::export]] 
Rcpp::List simulateTumorMTBPcpp(Rcpp::List input) {
    SimUtils::initMTBP(input); 
    Rcpp::List out = Sims::simulateMTBP(input);
    return out; 
}