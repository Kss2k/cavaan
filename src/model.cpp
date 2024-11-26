#include "cavaan.h"
#include <cstdio>


MatricesGroup *createMatricesGroup(Rcpp::List submodel) {
  Rcpp::List matrices = submodel["matrices"];
  MatricesGroup *matricesGroup = new MatricesGroup();
 
  // Initialize
  matricesGroup->BStar     = Rcpp::as<arma::mat>(matrices["BStar"]);
  matricesGroup->GammaStar = Rcpp::as<arma::mat>(matrices["GammaStar"]);
  matricesGroup->Phi       = Rcpp::as<arma::mat>(matrices["Phi"]);
  matricesGroup->Tau       = Rcpp::as<arma::mat>(matrices["Tau"]);
  matricesGroup->BStarInv  = Rcpp::as<arma::mat>(matrices["BStarInv"]);
  matricesGroup->G         = Rcpp::as<arma::mat>(matrices["G"]);
  matricesGroup->S         = Rcpp::as<arma::mat>(matrices["S"]);
  matricesGroup->Mu        = Rcpp::as<arma::mat>(matrices["Mu"]);
  matricesGroup->Nu        = Rcpp::as<arma::mat>(matrices["Nu"]);

  matricesGroup->p =   Rcpp::as<int>(matrices["p"]);

  return matricesGroup;
}


MatricesGroup *copyMatricesGroup(MatricesGroup *matricesGroup, bool fillZero = true) {
  MatricesGroup *copyMatricesGroup = new MatricesGroup();
  
  copyMatricesGroup->BStar     = matricesGroup->BStar;
  copyMatricesGroup->GammaStar = matricesGroup->GammaStar;
  copyMatricesGroup->Phi       = matricesGroup->Phi;
  copyMatricesGroup->Tau       = matricesGroup->Tau;
  copyMatricesGroup->BStarInv  = matricesGroup->BStarInv;
  copyMatricesGroup->G         = matricesGroup->G;
  copyMatricesGroup->S         = matricesGroup->S;
  copyMatricesGroup->Sigma     = matricesGroup->Sigma;
  copyMatricesGroup->Mu        = matricesGroup->Mu;
  copyMatricesGroup->Nu        = matricesGroup->Nu;


  if (fillZero) {
    copyMatricesGroup->BStar.fill(0);
    copyMatricesGroup->GammaStar.fill(0);
    copyMatricesGroup->Phi.fill(0);
    copyMatricesGroup->Tau.fill(0);
    copyMatricesGroup->BStarInv.fill(0);
    copyMatricesGroup->G.fill(0);
    copyMatricesGroup->S.fill(0);
    copyMatricesGroup->Sigma.fill(0);
    copyMatricesGroup->Nu.fill(0);
    copyMatricesGroup->Mu.fill(0);
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

  parTable->lhs    = Rcpp::as<std::vector<std::string>>(rParTable["lhs"]);
  parTable->op     = Rcpp::as<std::vector<std::string>>(rParTable["op"]);
  parTable->rhs    = Rcpp::as<std::vector<std::string>>(rParTable["rhs"]);
  parTable->est    = Rcpp::as<arma::vec>(rParTable["est"]);
  parTable->label  = Rcpp::as<std::vector<std::string>>(rParTable["label"]);
  parTable->col    = Rcpp::as<std::vector<int>>(rParTable["col"]);
  parTable->row    = Rcpp::as<std::vector<int>>(rParTable["row"]);
  parTable->matrix = Rcpp::as<std::vector<int>>(rParTable["matrix"]);
  parTable->group  = Rcpp::as<std::vector<int>>(rParTable["group"]);
  parTable->free   = Rcpp::as<std::vector<bool>>(rParTable["free"]);
  parTable->fill   = Rcpp::as<std::vector<bool>>(rParTable["fill"]);
  parTable->continueFromLast = Rcpp::as<std::vector<bool>>(rParTable["continue"]);
  parTable->isEquation = Rcpp::as<std::vector<bool>>(rParTable["isEquation"]);

  parTable->row   = reindexToZero(parTable->row);
  parTable->col   = reindexToZero(parTable->col);
  parTable->group = reindexToZero(parTable->group);

  parTable->nfree = countFree(parTable->free);

  parTable->expressions = std::vector<Expression*>();
  for (int i = 0; i < (int)(parTable->free.size()); i++) {
    if (!parTable->isEquation[i]) continue;
    parTable->expressions.push_back(createExpression(parTable->rhs[i]));
  }
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

  // getBaseGradients(m);

  return m;
}


void fillMatricesGroups(std::vector<MatricesGroup*> matricesGroups, ParTable *parTable, 
    const arma::vec &theta, bool calcSigma = true, bool fillConst = true) {

  MatricesGroup *matrices;
  int t = 0, e = 0, row, col, group;
  double tp;
  Rcpp::List evaluatedParams;

  for (int i = 0; i < (int)(parTable->fill.size()); i++) {
    if ((!parTable->fill[i] && !parTable->isEquation[i]) || 
        (!parTable->continueFromLast[i] && !parTable->free[i] && !fillConst)) continue;
    if (parTable->isEquation[i]) {
      tp = evaluateExpression(parTable->expressions[e++], evaluatedParams);
      evaluatedParams[parTable->label[i]] = tp;
      continue;
    } else if (!parTable->continueFromLast[i] && parTable->free[i]) {
      tp = theta[t++];
      evaluatedParams[parTable->label[i]] = tp;
    } else if (!parTable->continueFromLast[i]) {
      tp = parTable->est[i];
      evaluatedParams[parTable->label[i]] = tp;
    }

    // if (!arma::is_finite(tp)) Rcpp::Rcout << "label: " << parTable->label[i] << ", val: " << tp << '\n';
    row   = parTable->row[i];
    col   = parTable->col[i];
    group = parTable->group[i];

    matrices = matricesGroups[group];

    switch (parTable->matrix[i]) {
      case BETA_STAR:  {matrices->BStar.at(row, col)     = tp; break;}
      case GAMMA_STAR: {matrices->GammaStar.at(row, col) = tp; break;}
      case PHI:        {matrices->Phi.at(row, col)       = tp; break;}
      case TAU:        {matrices->Tau.at(row, col)       = tp; break;}
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
      matrices->Mu = matrices->G * BStarInv * matrices->GammaStar * matrices->Tau;
    } 
  }

}


void fillModel(Model *model, const arma::vec &theta, bool replace = false,
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
  Rcpp::XPtr<Model> xptr(model);
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
  return getGradientModelSimple(theta, model);
}


// [[Rcpp::export]]
Rcpp::NumericVector logLikCpp(const arma::vec &theta, SEXP xptr) {
  Rcpp::XPtr<Model> model(xptr);
  fillModel(model, theta, false, true);

  double logLik = 0, valLDSigma, signLDSigma, valLDS, signLDS;
  bool okLDSigma, okLDS, okSigmaInv;
  arma::mat SigmaInv, S, Sigma, Nu, Mu;
  int p;

  MatricesGroup* matrices;

  for (int i = 0; i < model->ngroups; i++) {
    matrices = model->matricesGroups[i];

    S     = matrices->S;
    Sigma = matrices->Sigma;
    Nu    = matrices->Nu;
    Mu    = matrices->Mu;
    p     = matrices->p;

    okLDS      = arma::log_det(valLDS, signLDS, S);
    okLDSigma  = arma::log_det(valLDSigma, signLDSigma, Sigma);
    okSigmaInv = arma::inv(SigmaInv, Sigma);

    if (!okLDS || !okLDSigma || !okSigmaInv) {
      return Rcpp::NumericVector::create(Rcpp::NumericVector::get_na());
    }
   
    logLik += valLDSigma + arma::trace(S * SigmaInv) - valLDS - p + 
      arma::as_scalar((Nu - Mu).t() * SigmaInv * (Nu - Mu));
  }

  return Rcpp::NumericVector::create(logLik);
}


// [[Rcpp::export]]
void debugCppModel(SEXP xptr, arma::vec theta) {
  Rcpp::XPtr<Model> model(xptr);
  fillModel(model, theta, false, true);
  Rcpp::Rcout << "Sigma: \n";
  Rcpp::Rcout << model->matricesGroups[0]->Sigma << '\n';
  Rcpp::Rcout << "SigmaInv: \n";
  Rcpp::Rcout << arma::inv(model->matricesGroups[0]->Sigma) << '\n';
  Rcpp::Rcout << "GammaStar: \n";
  Rcpp::Rcout << model->matricesGroups[0]->GammaStar << '\n';
  Rcpp::Rcout << "BStar: \n";
  Rcpp::Rcout << model->matricesGroups[0]->BStar << '\n';
  Rcpp::Rcout << "BStarInv: \n";
  Rcpp::Rcout << model->matricesGroups[0]->BStarInv << '\n';
  Rcpp::Rcout << "Phi: \n";
  Rcpp::Rcout << model->matricesGroups[0]->Phi << '\n';
  Rcpp::Rcout << "Mu: \n";
  Rcpp::Rcout << model->matricesGroups[0]->Mu << '\n';
  Rcpp::Rcout << "Nu: \n";
  Rcpp::Rcout << model->matricesGroups[0]->Nu << '\n';
}
