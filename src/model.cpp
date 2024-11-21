#include <RcppArmadillo.h>
#include <vector>
#include "model.h"
// [[Rcpp::depends(RcppArmadillo)]]


MatricesGroup *createMatricesGroup(Rcpp::List matrices) {
  MatricesGroup *mg = new MatricesGroup();

  return mg;
}


Model *createModel(Rcpp::List model) {
  Model *m = new Model();

  Rcpp::NumericVector groups = model["groups"];
  m->groups = groups;

  MatricesGroup *mg = createMatricesGroup(model);
  m->models.push_back(mg);

  return m;
}


// [[Rcpp::export]]
Rcpp::NumericVector ViewModelCreation(Rcpp::List model) {
  Model *m = createModel(model);
  Rcpp::Rcout << m->models[0]->GammaStar << '\n';

  return m->groups;  
}
