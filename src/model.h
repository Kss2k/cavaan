#ifndef model_h
#define model_h


#include <RcppArmadillo.h>
#include <vector>
// [[Rcpp::depends(RcppArmadillo)]]


#define BETA_STAR 0
#define GAMMA_STAR 1
#define PHI 2


// Matrices for a single group in a SEM
typedef struct {
  arma::mat *BStar;
  arma::mat *GammaStar;
  arma::mat *Phi;
  arma::mat *BStarInv;
  arma::mat *G;
  arma::mat *S;
  arma::mat *Sigma;
  int p;
} MatricesGroup;


// ParTable
typedef struct {
  arma::vec *est;
  std::vector<int>  *matrix;
  std::vector<int>  *group;
  std::vector<int>  *row;
  std::vector<int>  *col;
  std::vector<bool> *free;
  std::vector<bool> *fill;
  std::vector<bool> *continueFromLast;
} ParTable;


// Model
typedef struct {
  std::vector<MatricesGroup*> models;
  ParTable *parTable;
  int ngroups;
  int p;
} Model;


MatricesGroup *createMatricesGroup(Rcpp::List matrices);
Model *createModel(Rcpp::List model);
ParTable *createParTable(Rcpp::List model);
void fillModel(Model *model, arma::vec theta, bool replace);


#endif
