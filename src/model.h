#ifndef model_h
#define model_h


// Matrices for a single group in a SEM
typedef struct {
  arma::mat *GammaStar;
  arma::mat *BetaStar;
  arma::mat *BetaStarInv;
  arma::mat *Phi;
  arma::mat *G;
  arma::mat *S;
  arma::mat *Sigma;
} MatricesGroup;


// Model
typedef struct {
  std::vector<MatricesGroup*> models;
  Rcpp::NumericVector groups;
  int p;
} Model;


MatricesGroup *createMatricesGroup(Rcpp::List matrices);
Model *createModel(Rcpp::List model);


#endif
