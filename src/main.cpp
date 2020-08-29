#include"simulations.h"

// [[Rcpp::export]]
Rcpp::List simulateTumorcpp(Rcpp::List input) {
    SimUtils::initIA(input); 
    Rcpp::List out = Sims::simulateIA(input);
    return out; 
}

// [[Rcpp::export]] 
Rcpp::List simulateTumorUDTcpp(Rcpp::List input) {
    SimUtils::initUDT(input); 
    Rcpp::List out = Sims::simulateUDT(input);
    return out; 
}