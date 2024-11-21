#include "model.h"
#include <cstdio>
#include <vector>


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

  mg->p =   Rcpp::as<int>(matrices["p"]);

  return mg;
}


void reindexToZero(std::vector<int> *x) {
  for (int i = 0; i < (int)(x->size()); i++) {
    (*x)[i] = (*x)[i] - 1;
  }
}


ParTable *createParTable(Rcpp::List rParTable) {
  ParTable *parTable = new ParTable();

  parTable->est = new arma::vec(Rcpp::as<arma::vec>(rParTable["est"]));
  parTable->col = new std::vector<int>(Rcpp::as<std::vector<int>>(rParTable["col"]));
  parTable->row = new std::vector<int>(Rcpp::as<std::vector<int>>(rParTable["row"]));
  parTable->matrix = new std::vector<int>(Rcpp::as<std::vector<int>>(rParTable["matrix"]));
  parTable->group = new std::vector<int>(Rcpp::as<std::vector<int>>(rParTable["group"]));
  parTable->free = new std::vector<bool>(Rcpp::as<std::vector<bool>>(rParTable["free"]));
  parTable->fill = new std::vector<bool>(Rcpp::as<std::vector<bool>>(rParTable["fill"]));
  parTable->continueFromLast = new std::vector<bool>(Rcpp::as<std::vector<bool>>(rParTable["continue"]));

  reindexToZero(parTable->row);
  reindexToZero(parTable->col);
  reindexToZero(parTable->group);

  return parTable;
}


Model *createModel(Rcpp::List model) {
  Model *m = new Model();

  m->ngroups = Rcpp::as<Rcpp::NumericVector>(model["groups"]).length();

  Rcpp::List models = model["models"];
  for (int i = 0; i < m->ngroups; i++) {
    MatricesGroup *mg = createMatricesGroup(models[i]);
    m->models.push_back(mg);
  }

  Rcpp::List parTable = model["parTable.d"];
  m->parTable = createParTable(parTable);

  return m;
}


void fillModel(Model *model, arma::vec &theta, bool replace = false,
               bool calcSigma = true) {
  ParTable *parTable = model->parTable;

  int t = 0, row, col, group;
  double tp = -999.0; 
  for (int i = 0; i < parTable->fill->size(); i++) {
    if (!(*parTable->fill)[i]) continue;

    if (!(*parTable->continueFromLast)[i] && (*parTable->free)[i]) {
      tp = theta[t++];
    } else if (!(*parTable->continueFromLast)[i]) {
      tp = (*parTable->est)[i];
    }
  
    row = (*parTable->row)[i];
    col = (*parTable->col)[i];
    group = (*parTable->group)[i];

    switch ((*parTable->matrix)[i]) {
      case BETA_STAR: {
        model->models[group]->BStar->at(row, col) = tp;
        break;
      }
      case GAMMA_STAR: {
        model->models[group]->GammaStar->at(row, col) = tp;
        break;
      }
      case PHI: {
        model->models[group]->Phi->at(row, col) = tp;
        break;
      }
      default:
        Rcpp::stop("Unrecognized matrix index");

    }
  } 

  if (calcSigma) {
    MatricesGroup *matrices;
    arma::mat BStarInv;
    arma::mat Sigma;

    for (int i = 0; i < model->ngroups; i++) {
      matrices = model->models[i];
      BStarInv = arma::inv(*matrices->BStar);
      
      delete matrices->BStarInv;
      delete matrices->Sigma;

      Sigma = (*matrices->G) * BStarInv * (*matrices->GammaStar) * 
        (*matrices->Phi) * (*matrices->GammaStar).t() * BStarInv.t() *
        (*matrices->G).t();

      matrices->BStarInv = new arma::mat(BStarInv);
      matrices->Sigma = new arma::mat(Sigma);
    } 
  }
}


// [[Rcpp::export]]
Rcpp::NumericVector ViewModelCreation(Rcpp::List model, arma::vec theta) {
  Model *m = createModel(model);
  
  fillModel(m, theta, false, true);
  Rcpp::Rcout << *(m->models[0]->Sigma) << '\n';

  return m->ngroups;  
}


// [[Rcpp::export]]
Rcpp::NumericVector logLikR2Cpp(arma::vec theta, Rcpp::List model) {
  // Create the model and fill it
  Model* m = createModel(model);
  fillModel(m, theta, false, true);

  double logLik = 0;
  double log_det_Sigma;
  double sign_Sigma;
  double log_det_S;
  double sign_S;
  arma::mat Sigma_inv;
  double trace_value;

  MatricesGroup* matrices;

  for (int i = 0; i < m->ngroups; i++) {
    matrices = m->models[i];

    // Compute log(det(Sigma))
    arma::log_det(log_det_Sigma, sign_Sigma, *matrices->Sigma);
    if (sign_Sigma <= 0) {
      return Rcpp::NumericVector::create(Rcpp::NumericVector::get_na());  // Return NaN
    }

    // Compute log(det(S))
    arma::log_det(log_det_S, sign_S, *matrices->S);
    if (sign_S <= 0) {
      return Rcpp::NumericVector::create(Rcpp::NumericVector::get_na());  // Return NaN
    }

    // Compute trace of solve(Sigma) %*% S
    Sigma_inv = arma::inv(*matrices->Sigma);
    trace_value = arma::trace(Sigma_inv * (*matrices->S));

    // Accumulate the log-likelihood
    logLik += log_det_Sigma + trace_value - log_det_S - matrices->p;
  }

  return Rcpp::NumericVector::create(logLik);  // Return log-likelihood as a single-element vector
}

// old version
// double logLikR2Cpp(arma::vec theta, Rcpp::List model) {
//   Model *m = createModel(model);
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
//   MatricesGroup *matrices;
//   for (int i = 0; i < m->ngroups; i++) {
//     matrices = m->models[i];
//     
//     arma::log_det(log_det_Sigma, sign_Sigma, *matrices->Sigma);
//     arma::log_det(log_det_S, sign_S, *matrices->S);
//     
//     Sigma_inv = arma::inv(*matrices->Sigma);
//     trace_value = arma::trace(Sigma_inv * (*matrices->S));
//     
//     logLik += log_det_Sigma + trace_value - log_det_S - matrices->p;
//   }
// 
//   return logLik;
// }
// 
