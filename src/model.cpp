#include <RcppArmadillo.h>
#include <vector>
#include "model.h"
// [[Rcpp::depends(RcppArmadillo)]]


MatricesGroup *createMatricesGroup(Rcpp::List submodel) {
  Rcpp::List matrices = submodel["matrices"];
  MatricesGroup *mg = new MatricesGroup();
 
  // Initialize
  mg->BStar     = new arma::mat(Rcpp::as<arma::mat>(matrices["BStar"]));
  mg->GammaStar = new arma::mat(Rcpp::as<arma::mat>(matrices["GammaStar"]));
  mg->Phi       = new arma::mat(Rcpp::as<arma::mat>(matrices["Phi"]));
  mg->BStarInv  = new arma::mat(Rcpp::as<arma::mat>(matrices["BStarInv"]));
  mg->G         = new arma::mat(Rcpp::as<arma::mat>(matrices["G"]));
  mg->S         = new arma::mat(Rcpp::as<arma::mat>(matrices["S"]));
  mg->Sigma     = new arma::mat(Rcpp::as<arma::mat>(matrices["Sigma"]));
  
  return mg;
}


Model *createModel(Rcpp::List model) {
  Model *m = new Model();

  m->ngroups = Rcpp::as<Rcpp::NumericVector>(model["groups"]).length();

  Rcpp::List models = model["models"];
  for (int i = 0; i < m->ngroups; i++) {
    MatricesGroup *mg = createMatricesGroup(models[i]);
    m->models.push_back(mg);
  }

  return m;
}


// [[Rcpp::export]]
Rcpp::NumericVector ViewModelCreation(Rcpp::List model) {
  Model *m = createModel(model);
  Rcpp::Rcout << *(m->models[0]->GammaStar) << '\n';

  return m->ngroups;  
}
