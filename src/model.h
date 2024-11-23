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
  arma::mat BStar;
  arma::mat GammaStar;
  arma::mat Phi;
  arma::mat BStarInv;
  arma::mat G;
  arma::mat S;
  arma::mat Sigma;
  int p;
} MatricesGroup;


// ParTable
typedef struct {
  arma::vec est;
  std::vector<int>  matrix;
  std::vector<int>  group;
  std::vector<int>  row;
  std::vector<int>  col;
  std::vector<bool> free;
  std::vector<bool> fill;
  std::vector<bool> continueFromLast;
  int nfree;
} ParTable;


// GradientMatricesParam
typedef struct {
  std::vector<MatricesGroup*> matricesGroups; // vector with elem for each group
} GradientMatricesParam;


typedef struct {
  int iterations;
  int maxIterations;
  double threshold;
  bool convergence;
  arma::vec par;
} OptimizerInfo;


// Model
typedef struct {
  std::vector<MatricesGroup*> matricesGroups;
  std::vector<GradientMatricesParam*> gradientMatricesParams; // vector with elem for each param
  ParTable *parTable;
  int ngroups;
  int p;
  OptimizerInfo optimizerInfo;
} Model;


// model.cpp
MatricesGroup *createMatricesGroup(Rcpp::List matrices);
MatricesGroup *copyMatricesGroup(MatricesGroup *matricesGroup, bool fillZero);
Model *createModel(Rcpp::List model);
ParTable *createParTable(Rcpp::List model);
void fillModel(Model *model, arma::vec &theta, bool replace, bool calcSigma);
void fillMatricesGroups(std::vector<MatricesGroup*> matricesGroups, ParTable *parTable, 
    arma::vec &theta, bool calcSigma, bool fillConst);
int countFree(std::vector<bool> free);


// gradient.cpp
void getBaseGradients(Model *model);
arma::vec getGradientModel(arma::vec theta, Model *model);
arma::vec normalize(arma::vec x);
double getLogLikModel(arma::vec theta, Model *model);


// optimize.cpp
arma::vec optim(arma::vec theta, Model* model, double (*objective)(arma::vec, Model*),
    arma::vec (*gradient)(arma::vec, Model*));
arma::vec optimBFGS(arma::vec theta, Model* model, double (*objective)(arma::vec, Model*),
    arma::vec (*gradient)(arma::vec, Model*));

#endif
