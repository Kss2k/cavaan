#include "cavaan.h"
#include <cmath>


arma::vec getGradientModelSimple(const arma::vec& theta, Model* model) {
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
 
    arma::mat I = arma::eye(Phi.n_rows, Phi.n_cols);   
    arma::mat SigmaInv;
    if (!arma::inv(SigmaInv, Sigma)) Rcpp::stop("Sigma not invertible");

    const arma::mat BStarInv_T  = BStarInv.t();
    const arma::mat GammaStar_T = GammaStar.t();
    const arma::mat G_T         = G.t();
    const arma::mat Q           = SigmaInv - SigmaInv * S * SigmaInv;
    

    // THIS IS NOT ENTIERLY CORRECT! -------------------------------------------
    // The derivatives fro Phi, Gamma and B assume that they do not affect
    // the meanstructure of the model! In practice, it seems like this is fine!
    // but it might affect parameter estimates for some more complicated examples
    // and more likely the std.error estimates obtained using the derivatives. 
    // DerivTau is correct however, as the mean-structure does not affect the 
    // other parts of the model.

    arma::mat DerivPhi = (GammaStar_T * BStarInv_T * G_T * Q * G * BStarInv * GammaStar);
    arma::mat DerivGammaStar = (BStarInv_T * G_T * Q * G * BStarInv * GammaStar * Phi);
    arma::mat DerivBStar = (BStarInv_T * G_T * Q * G * BStarInv * 
        GammaStar * Phi * GammaStar_T * BStarInv_T);
    arma::mat DerivTau = (G * BStarInv * GammaStar).t() * SigmaInv * (Nu - Mu);
    
    // -------------------------------------------------------------------------
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
 
    arma::mat I = arma::eye(Phi.n_rows, Phi.n_cols);   
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

    const arma::mat DerivBStar1 = - 2 * (BStarInv_T * G_T * Q * G * BStarInv * 
        GammaStar * Phi * GammaStar_T * BStarInv_T);


   // const arma::mat DerivBStar2 = BStarInv * G_T * SigmaInv * Snm * Tau.t() * 
   //   GammaStar_T * BStarInv + 2 * BStarInv * G_T * SigmaInv * Sigma * SigmaInv *
   //   Snm * Snm.t() * SigmaInv * G * BStarInv;
     // 2B−1Γτ(H′M−1S)′−[M−1SS′M−1HΓΦΓ′B−1+(M−1SS′M−1HΓΦΓ′B−1)′]
    const arma::mat DerivBStar2 = 2 * BStarInv * GammaStar * Tau.t() * SigmaInv * H.t() - 
      (SigmaInv * S * S.t() * SigmaInv * H * GammaStar * Phi * GammaStar.t() * BStarInv + 
       SigmaInv * S * S.t() * SigmaInv * H * GammaStar * Phi * GammaStar.t() * BStarInv).t();

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
