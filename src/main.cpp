#include"simulations.h"

// [[Rcpp::export]]
Rcpp::List simulateTumorcpp(Rcpp::List input) {
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