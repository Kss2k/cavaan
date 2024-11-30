#include "cavaan.h"
#include <cmath>


arma::vec getGradientModelNoMeans(const arma::vec& theta, Model* model) {
  fillModel(model, theta, false, true);

  int npar = model->parTable->nfree;
  int ngroups = model->ngroups;

  arma::vec grad(npar, arma::fill::zeros);

  ParTable *parTable = model->parTable;
  for (int g = 0; g < ngroups; g++) {
    MatricesGroup* matricesGroup = model->matricesGroups[g];

    const arma::mat& Sigma     = matricesGroup->Sigma;
    const arma::mat& G         = matricesGroup->G;
    const arma::mat& S         = matricesGroup->S;
    const arma::mat& Phi       = matricesGroup->Phi;
    const arma::mat& GammaStar = matricesGroup->GammaStar;
    const arma::mat& BStarInv  = matricesGroup->BStarInv;
    const arma::mat  I = arma::eye(Phi.n_rows, Phi.n_cols);   

    arma::mat SigmaInv;
    if (!arma::inv(SigmaInv, Sigma)) Rcpp::stop("Sigma not invertible");

    const arma::mat BStarInv_T  = BStarInv.t();
    const arma::mat GammaStar_T = GammaStar.t();
    const arma::mat G_T         = G.t();
    const arma::mat Q           = SigmaInv - SigmaInv * S * SigmaInv;

    const arma::mat DerivPhi = (GammaStar_T * BStarInv_T * G_T * Q * G * BStarInv * GammaStar);
    const arma::mat DerivGammaStar = (BStarInv_T * G_T * Q * G * BStarInv * GammaStar * Phi);
    const arma::mat DerivBStar = (BStarInv_T * G_T * Q * G * BStarInv * 
        GammaStar * Phi * GammaStar_T * BStarInv_T);
    
    int t = 0;
    for (int i = 0; i < (int)(parTable->free.size()) && t < npar; i++) {
      if (!parTable->free[i]) continue;

      int row = parTable->row[i];
      int col = parTable->col[i];

      switch (parTable->matrix[i]) {
        case BETA_STAR:  {grad[t++] += -2 * DerivBStar.at(row, col); break;}
        case GAMMA_STAR: {grad[t++] +=  2 * DerivGammaStar.at(row, col); break;}
        case PHI:        {grad[t++] += (2 - I.at(row, col)) * DerivPhi.at(row, col); break;}
        default: Rcpp::stop("Unrecognized matrix index");
      }
    }
  }


  return grad;
}


arma::vec getGradientModelOVMeans(const arma::vec& theta, Model* model) {
  fillModel(model, theta, false, true);

  int npar = model->parTable->nfree;
  int ngroups = model->ngroups;

  arma::vec grad(npar, arma::fill::zeros);

  ParTable *parTable = model->parTable;
  for (int g = 0; g < ngroups; g++) {
    MatricesGroup* matricesGroup = model->matricesGroups[g];

    const arma::mat& Sigma     = matricesGroup->Sigma;
    const arma::mat& G         = matricesGroup->G;
    const arma::mat& S         = matricesGroup->S;
    const arma::mat& Phi       = matricesGroup->Phi;
    const arma::mat& GammaStar = matricesGroup->GammaStar;
    const arma::mat& BStarInv  = matricesGroup->BStarInv;
    const arma::mat& Nu        = matricesGroup->Nu;
    const arma::mat& Mu        = matricesGroup->Mu;
 
    const arma::mat I = arma::eye(Phi.n_rows, Phi.n_cols);   
    arma::mat SigmaInv;
    if (!arma::inv(SigmaInv, Sigma)) Rcpp::stop("Sigma not invertible");

    const arma::mat BStarInv_T  = BStarInv.t();
    const arma::mat GammaStar_T = GammaStar.t();
    const arma::mat G_T         = G.t();
    const arma::mat Q           = SigmaInv - SigmaInv * S * SigmaInv;

    const arma::mat DerivPhi = (GammaStar_T * BStarInv_T * G_T * Q * G * BStarInv * GammaStar);
    const arma::mat DerivGammaStar = (BStarInv_T * G_T * Q * G * BStarInv * GammaStar * Phi);
    const arma::mat DerivBStar = (BStarInv_T * G_T * Q * G * BStarInv * 
        GammaStar * Phi * GammaStar_T * BStarInv_T);
    const arma::mat DerivTau = (G * BStarInv * GammaStar).t() * SigmaInv * (Nu - Mu);
    
    int t = 0;
    for (int i = 0; i < (int)(parTable->free.size()) && t < npar; i++) {
      if (!parTable->free[i]) continue;

      int row = parTable->row[i];
      int col = parTable->col[i];

      switch (parTable->matrix[i]) {
        case BETA_STAR:  {grad[t++] += -2 * DerivBStar.at(row, col); break;}
        case GAMMA_STAR: {grad[t++] +=  2 * DerivGammaStar.at(row, col); break;}
        case PHI:        {grad[t++] += (2 - I.at(row, col)) * DerivPhi.at(row, col); break;}
        case TAU:        {grad[t++] += -2 * DerivTau.at(row, 0); break;}
        default: Rcpp::stop("Unrecognized matrix index");
      }
    }
  }


  return grad;
}


arma::vec getGradientModelLVMeans(const arma::vec& theta, Model* model) {
  fillModel(model, theta, false, true);

  int npar = model->parTable->nfree;
  int ngroups = model->ngroups;

  arma::vec grad(npar, arma::fill::zeros);

  ParTable *parTable = model->parTable;
  for (int g = 0; g < ngroups; g++) {
    MatricesGroup* matricesGroup = model->matricesGroups[g];

    const arma::mat& Sigma     = matricesGroup->Sigma;
    const arma::mat& G         = matricesGroup->G;
    const arma::mat& S         = matricesGroup->S;
    const arma::mat& Phi       = matricesGroup->Phi;
    const arma::mat& GammaStar = matricesGroup->GammaStar;
    const arma::mat& BStarInv  = matricesGroup->BStarInv;
    const arma::mat& Tau       = matricesGroup->Tau;
    const arma::mat& Nu        = matricesGroup->Nu;
    const arma::mat& Mu        = matricesGroup->Mu;
 
    const arma::mat I = arma::eye(Phi.n_rows, Phi.n_cols);   
    arma::mat SigmaInv;
    if (!arma::inv(SigmaInv, Sigma)) Rcpp::stop("Sigma not invertible");

    const arma::mat BStarInv_T  = BStarInv.t();
    const arma::mat GammaStar_T = GammaStar.t();
    const arma::mat G_T         = G.t();
    const arma::mat Q           = SigmaInv - SigmaInv * S * SigmaInv;
   
    const arma::mat DerivTau = (G * BStarInv * GammaStar).t() * SigmaInv * (Nu - Mu);
    const arma::mat M = Sigma;
    const arma::mat Snm = Nu - Mu;
    const arma::mat K = G * BStarInv;
    const arma::mat K_T = K.t();
    const arma::mat H = G * BStarInv;
    const arma::mat H_T = H.t();

    const arma::mat DerivBStar1 = - 2 * (BStarInv_T * G_T * Q * G * BStarInv * 
        GammaStar * Phi * GammaStar_T * BStarInv_T);
    const arma::mat DerivBStar2 = - 2 * BStarInv * (
        GammaStar * 
        (Tau - Phi * GammaStar_T * H_T * SigmaInv * Snm) * (G_T * SigmaInv * Snm).t()
        ) * BStarInv;
    const arma::mat DerivBStar = DerivBStar1 + DerivBStar2;

    const arma::mat DerivGammaStar1 = 2 * (BStarInv_T * G_T * Q * G * BStarInv * GammaStar * Phi);
    const arma::mat DerivGammaStar2 = - 2 * H.t() * SigmaInv * Snm * Tau.t() - 2 * H.t() * SigmaInv *
      Snm * Snm.t() * SigmaInv * H * GammaStar * Phi;

    const arma::mat DerivGammaStar = DerivGammaStar1 + DerivGammaStar2;
    
    const arma::mat DerivPhi1 = (GammaStar_T * BStarInv_T * G_T * Q * G * BStarInv * GammaStar);
    const arma::mat DerivPhi2 = - GammaStar_T * H.t() * SigmaInv * Snm * Snm.t() * 
      SigmaInv * H * GammaStar;
    const arma::mat DerivPhi = DerivPhi1 + DerivPhi2;

    
    int t = 0;
    for (int i = 0; i < (int)(parTable->free.size()) && t < npar; i++) {
      if (!parTable->free[i]) continue;

      int row = parTable->row[i];
      int col = parTable->col[i];

      switch (parTable->matrix[i]) {
        case BETA_STAR:  {grad[t++] += DerivBStar.at(row, col); break;}
        case GAMMA_STAR: {grad[t++] += DerivGammaStar.at(row, col); break;}
        case PHI:        {grad[t++] += (2 - I.at(row, col)) * DerivPhi.at(row, col); break;}
        case TAU:        {grad[t++] += -2 * DerivTau.at(row, 0); break;}
        default: Rcpp::stop("Unrecognized matrix index");
      }
    }
  }


  return grad;
}


// [[Rcpp::export]]
arma::vec gradLogLikCppNoMeans(arma::vec theta, SEXP xptr) {
  Rcpp::XPtr<Model> model(xptr);
  return getGradientModelOVMeans(theta, model);
}


// [[Rcpp::export]]
arma::vec gradLogLikCppOVMeans(arma::vec theta, SEXP xptr) {
  Rcpp::XPtr<Model> model(xptr);
  return getGradientModelOVMeans(theta, model);
}


// [[Rcpp::export]]
arma::vec gradLogLikCppLVMeans(arma::vec theta, SEXP xptr) {
  Rcpp::XPtr<Model> model(xptr);
  return getGradientModelLVMeans(theta, model);
}


// [[Rcpp::export]]
Rcpp::NumericVector gradLogLikNumericCpp(const arma::vec& theta, SEXP xptr, double h=1e-6) {
  Rcpp::XPtr<Model> model(xptr);
  Rcpp::NumericVector grad(theta.size());
  Rcpp::NumericVector logLikBase = logLikCpp(theta, model);
  Rcpp::NumericVector logLikPlus;

  for (int i = 0; i < (int)theta.size(); i++) {
    arma::vec thetaPlus = theta;
    thetaPlus[i] += h;
    logLikPlus = logLikCpp(thetaPlus, model);
    grad[i] = (logLikPlus[0] - logLikBase[0]) / h; 
  }

  return grad;
}
