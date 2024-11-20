#include "model.h"


MatricesGroup *createMatricesGroup(Rcpp::List matrices) {
  MatricesGroup *mg = new MatricesGroup();

  return mg;
}


Model *createModel(Rcpp::List model) {
  Model *m = new Model();
  MatricesGroup *mg = createMatricesGroup(model);
  m->models.push_back(mg);
  return m;
}


// [[Rcpp::export]]
Rcpp::NumericVector addOne(Rcpp::NumericVector x) {
  return x + 1;
}
