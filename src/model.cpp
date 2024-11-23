#include "model.h"
#include <cstdio>


MatricesGroup *createMatricesGroup(Rcpp::List submodel) {
  Rcpp::List matrices = submodel["matrices"];
  MatricesGroup *matricesGroup = new MatricesGroup();
 
  // Initialize
  matricesGroup->BStar     = Rcpp::as<arma::mat>(matrices["BStar"]);
  matricesGroup->GammaStar = Rcpp::as<arma::mat>(matrices["GammaStar"]);
  matricesGroup->Phi       = Rcpp::as<arma::mat>(matrices["Phi"]);
  matricesGroup->BStarInv  = Rcpp::as<arma::mat>(matrices["BStarInv"]);
  matricesGroup->G         = Rcpp::as<arma::mat>(matrices["G"]);
  matricesGroup->S         = Rcpp::as<arma::mat>(matrices["S"]);
  matricesGroup->Sigma     = Rcpp::as<arma::mat>(matrices["Sigma"]);

  matricesGroup->p =   Rcpp::as<int>(matrices["p"]);

  return matricesGroup;
}


MatricesGroup *copyMatricesGroup(MatricesGroup *matricesGroup, bool fillZero = true) {
  MatricesGroup *copyMatricesGroup = new MatricesGroup();
  
  copyMatricesGroup->BStar     = matricesGroup->BStar;
  copyMatricesGroup->GammaStar = matricesGroup->GammaStar;
  copyMatricesGroup->Phi       = matricesGroup->Phi;
  copyMatricesGroup->BStarInv  = matricesGroup->BStarInv;
  copyMatricesGroup->G         = matricesGroup->G;
  copyMatricesGroup->S         = matricesGroup->S;
  copyMatricesGroup->Sigma     = matricesGroup->Sigma;


  if (fillZero) {
    copyMatricesGroup->BStar.fill(0);
    copyMatricesGroup->GammaStar.fill(0);
    copyMatricesGroup->Phi.fill(0);
    copyMatricesGroup->BStarInv.fill(0);
    copyMatricesGroup->G.fill(0);
    copyMatricesGroup->S.fill(0);
    copyMatricesGroup->Sigma.fill(0);
  }

  return copyMatricesGroup;
}


std::vector<int> reindexToZero(std::vector<int> x) {
  for (int i = 0; i < (int)(x.size()); i++) {
    x[i] = x[i] - 1;
  }
  return x;
}


int countFree(std::vector<bool> free) {
  int n = 0;
  for (int i = 0; i < (int)(free.size()); i++) {
    if (free[i]) n++;
  }

  return n;
}


ParTable *createParTable(Rcpp::List rParTable) {
  ParTable *parTable = new ParTable();

  parTable->est = Rcpp::as<arma::vec>(rParTable["est"]);
  parTable->col = Rcpp::as<std::vector<int>>(rParTable["col"]);
  parTable->row = Rcpp::as<std::vector<int>>(rParTable["row"]);
  parTable->matrix = Rcpp::as<std::vector<int>>(rParTable["matrix"]);
  parTable->group  = Rcpp::as<std::vector<int>>(rParTable["group"]);
  parTable->free   = Rcpp::as<std::vector<bool>>(rParTable["free"]);
  parTable->fill   = Rcpp::as<std::vector<bool>>(rParTable["fill"]);
  parTable->continueFromLast = Rcpp::as<std::vector<bool>>(rParTable["continue"]);

  parTable->row   = reindexToZero(parTable->row);
  parTable->col   = reindexToZero(parTable->col);
  parTable->group = reindexToZero(parTable->group);

  parTable->nfree = countFree(parTable->free);

  return parTable;
}


Model *createModel(Rcpp::List model) {
  Model *m = new Model();

  m->ngroups = Rcpp::as<Rcpp::NumericVector>(model["groups"]).length();

  Rcpp::List matricesGroups = model["models"];
  for (int i = 0; i < m->ngroups; i++) {
    MatricesGroup *matricesGroup = createMatricesGroup(matricesGroups[i]);
    m->matricesGroups.push_back(matricesGroup);
  }

  Rcpp::List parTable = model["parTable.d"];
  m->parTable = createParTable(parTable);

  getBaseGradients(m);

  return m;
}


void fillMatricesGroups(std::vector<MatricesGroup*> matricesGroups, ParTable *parTable, 
    arma::vec &theta, bool calcSigma = true, bool fillConst = true) {

  MatricesGroup *matrices;
  int t = 0, row, col, group;
  double tp;

  for (int i = 0; i < (int)(parTable->fill.size()); i++) {
    if (!parTable->fill[i] || (!parTable->continueFromLast[i] && 
        !parTable->free[i] && !fillConst)) continue;

    if      (!parTable->continueFromLast[i] && parTable->free[i]) tp = theta[t++];
    else if (!parTable->continueFromLast[i]) tp = parTable->est[i];
  
    row   = parTable->row[i];
    col   = parTable->col[i];
    group = parTable->group[i];

    matrices = matricesGroups[group];

    switch (parTable->matrix[i]) {
      case BETA_STAR:  {matrices->BStar.at(row, col)     = tp; break;}
      case GAMMA_STAR: {matrices->GammaStar.at(row, col) = tp; break;}
      case PHI:        {matrices->Phi.at(row, col)       = tp; break;}
      default: Rcpp::stop("Unrecognized matrix index");
    }
  } 
  
  if (calcSigma) {
    arma::mat BStarInv, Sigma;

    for (int i = 0; i < (int)(matricesGroups.size()); i++) {
      matrices = matricesGroups[i];
      BStarInv = arma::inv(matrices->BStar);
      
      Sigma = matrices->G * BStarInv * matrices->GammaStar * 
        matrices->Phi * matrices->GammaStar.t() * BStarInv.t() *
        matrices->G.t();

      matrices->BStarInv = BStarInv;
      matrices->Sigma = Sigma;
    } 
  }

}


void fillModel(Model *model, arma::vec &theta, bool replace = false,
               bool calcSigma = true) {
  fillMatricesGroups(model->matricesGroups, model->parTable, theta, true);

  // replace is not implemented yet!!
}


// [[Rcpp::export]]
Rcpp::NumericVector ViewModelCreation(Rcpp::List RModel, arma::vec theta) {
  Model *model = createModel(RModel);
  
  fillModel(model, theta, false, true);
  Rcpp::Rcout << model->matricesGroups[0]->Sigma << '\n';
  
  getBaseGradients(model);
  Rcpp::Rcout << model->gradientMatricesParams[0]->matricesGroups[0]->GammaStar << '\n';

  return model->ngroups;  
}


// [[Rcpp::export]]
RcppExport SEXP createRcppModel(Rcpp::List RModel) {
  Model *model = createModel(RModel);
  Rcpp::XPtr xptr(model);
  return xptr;
}


// [[Rcpp::export]]
RcppExport SEXP fillRcppModel(SEXP xptr, arma::vec theta) {
  Rcpp::XPtr<Model> model(xptr);
  fillModel(model, theta, false, true);
  Rcpp::XPtr<Model> rmodel(model);
  return rmodel;
}


// [[Rcpp::export]]
arma::vec gradLogLikCpp(arma::vec theta, SEXP xptr) {
  Rcpp::XPtr<Model> model(xptr);
  return getGradientModel(theta, model);
}


// [[Rcpp::export]]
Rcpp::NumericVector logLikCpp(arma::vec theta, SEXP xptr) {
  Rcpp::XPtr<Model> model(xptr);
  fillModel(model, theta, false, true);

  double logLik = 0, valLDSigma, signLDSigma, valLDS, signLDS;
  bool okLDSigma, okLDS, okSigmaInv;
  arma::mat SigmaInv, S, Sigma;
  int p;

  MatricesGroup* matrices;

  for (int i = 0; i < model->ngroups; i++) {
    matrices = model->matricesGroups[i];

    S     = matrices->S;
    Sigma = matrices->Sigma;
    p     = matrices->p;

    okLDS      = arma::log_det(valLDS, signLDS, S);
    okLDSigma  = arma::log_det(valLDSigma, signLDSigma, Sigma);
    okSigmaInv = arma::inv(SigmaInv, Sigma);

    if (!okLDS || !okLDSigma || !okSigmaInv || signLDSigma < 0 || signLDS < 0) {
      return Rcpp::NumericVector::create(Rcpp::NumericVector::get_na());
    }
    
    logLik += valLDSigma + arma::trace(S * SigmaInv) - valLDS - p;
  }

  return Rcpp::NumericVector::create(logLik);
}


// // [[Rcpp::export]]
// Rcpp::NumericVector logLikR2Cpp(arma::vec theta, Rcpp::List model) {
//   // Create the model and fill it
//   Model* m = createModel(model);
//   fillModel(m, theta, false, true);
// 
//   double logLik = 0;
//   double log_det_Sigma;
//   double sign_Sigma;
//   double log_det_S;
//   double sign_S;
//   arma::mat Sigma_inv;
//   double trace_value;
// 
//   MatricesGroup* matrices;
// 
//   for (int i = 0; i < m->ngroups; i++) {
//     matrices = m->matricesGroups[i];
// 
//     // Compute log(det(Sigma))
//     arma::log_det(log_det_Sigma, sign_Sigma, matrices->Sigma);
//     if (sign_Sigma <= 0) {
//       return Rcpp::NumericVector::create(Rcpp::NumericVector::get_na());  // Return NaN
//     }
// 
//     // Compute log(det(S))
//     arma::log_det(log_det_S, sign_S, matrices->S);
//     if (sign_S <= 0) {
//       return Rcpp::NumericVector::create(Rcpp::NumericVector::get_na());  // Return NaN
//     }
// 
//     // Compute trace of solve(Sigma) %*% S
//     Sigma_inv = arma::inv(matrices->Sigma);
//     trace_value = arma::trace(Sigma_inv * matrices->S);
// 
//     // Accumulate the log-likelihood
//     logLik += log_det_Sigma + trace_value - log_det_S - matrices->p;
//   }
// 
//   return Rcpp::NumericVector::create(logLik);  // Return log-likelihood as a single-element vector
// }
// 
// 
// // [[Rcpp::export]]
// arma::vec gradLogLikR2Cpp(arma::vec theta, Rcpp::List model) {
//   Model* m = createModel(model);
//   fillModel(m, theta, false, true);
// 
//   return getGradientModel(theta, m);
// }
// 
// 
// // [[Rcpp::export]]
// arma::vec semR2Cpp(arma::vec theta, Rcpp::List rmodel) {
//   Model *model = createModel(rmodel);
//   
//   return optimBFGS(theta, model, getLogLikModel, getGradientModel);
// }
