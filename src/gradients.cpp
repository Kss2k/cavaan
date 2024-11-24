#include "cavaan.h"
#include <cmath>
#include <stdexcept>


void getBaseGradients(Model *model) {
  GradientMatricesParam *gradientMatricesParam; 
  MatricesGroup *matricesGroupParam;
  arma::vec theta = arma::vec(model->parTable->nfree).fill(0);

  for (int i = 0; i < model->parTable->nfree; i++) {
    gradientMatricesParam = new GradientMatricesParam();

    theta[i] = 1;

    for (int g = 0; g < (int)(model->matricesGroups.size()); g++) {
      matricesGroupParam = copyMatricesGroup(model->matricesGroups[g], true);
      gradientMatricesParam->matricesGroups.push_back(matricesGroupParam);
    }
    
    fillMatricesGroups(gradientMatricesParam->matricesGroups, model->parTable, theta, false, false);
    theta[i] = 0; // reset before next iteration
    
    model->gradientMatricesParams.push_back(gradientMatricesParam);
  }
}


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
 
    arma::mat I = arma::eye(Phi.n_rows, Phi.n_cols);   
    arma::mat SigmaInv;
    if (!arma::inv(SigmaInv, Sigma)) Rcpp::stop("Sigma not invertible");

    const arma::mat BStarInv_T  = BStarInv.t();
    const arma::mat GammaStar_T = GammaStar.t();
    const arma::mat G_T         = G.t();
    const arma::mat Q           = SigmaInv - SigmaInv * S * SigmaInv;


    arma::mat DerivPhi = (GammaStar_T * BStarInv_T * G_T * Q * G * BStarInv * GammaStar);
    arma::mat DerivGammaStar = (BStarInv_T * G_T * Q * G * BStarInv * GammaStar * Phi);
    arma::mat DerivBStar = (BStarInv_T * G_T * Q * G * BStarInv * 
        GammaStar * Phi * GammaStar_T * BStarInv_T);

    int t = 0;
    for (int i = 0; i < parTable->free.size() && t < npar; i++) {
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


arma::vec getGradientModelGeneral(const arma::vec& theta, Model* model) {
  fillModel(model, theta, false, true);

  int npar = model->parTable->nfree;
  int ngroups = model->ngroups;

  arma::vec grad(npar, arma::fill::zeros);

  for (int g = 0; g < ngroups; g++) {
    MatricesGroup* matricesGroup = model->matricesGroups[g];

    const arma::mat& Sigma     = matricesGroup->Sigma;
    const arma::mat& G         = matricesGroup->G;
    const arma::mat& S         = matricesGroup->S;
    const arma::mat& Phi       = matricesGroup->Phi;
    const arma::mat& GammaStar = matricesGroup->GammaStar;
    const arma::mat& BStarInv  = matricesGroup->BStarInv;
    
    arma::mat SigmaInv;
    if (!arma::inv(SigmaInv, Sigma)) Rcpp::stop("Sigma not invertible");

    const arma::mat BStarInv_T  = BStarInv.t();
    const arma::mat GammaStar_T = GammaStar.t();
    const arma::mat Q           = SigmaInv - SigmaInv * S * SigmaInv;

    for (int t = 0; t < npar; t++) {
      GradientMatricesParam* gradientMatricesParam = model->gradientMatricesParams[t];
      MatricesGroup* gradientMatricesParamGroup = gradientMatricesParam->matricesGroups[g];

      const arma::mat& DerivGammaStar  = gradientMatricesParamGroup->GammaStar;
      const arma::mat& DerivBStar      = gradientMatricesParamGroup->BStar;
      const arma::mat& DerivPhi        = gradientMatricesParamGroup->Phi;
      const arma::mat DerivGammaStar_T = DerivGammaStar.t();
      const arma::mat DerivBStarInv    = -BStarInv * DerivBStar * BStarInv;

      arma::mat temp = BStarInv * GammaStar;
      arma::mat temp_T = temp.t();
      arma::mat DerivSigma = G * (
          DerivBStarInv * GammaStar * Phi * GammaStar_T * BStarInv_T +
          BStarInv * DerivGammaStar * Phi * GammaStar_T * BStarInv_T +
          temp * DerivPhi * GammaStar_T * BStarInv_T +
          temp * Phi * DerivGammaStar_T * BStarInv_T +
          temp * Phi * GammaStar_T * DerivBStarInv.t()
          ) * G.t();

      grad[t] += arma::trace(Q * DerivSigma);
    }
  }

  return grad;
}

arma::vec normalize(arma::vec x) {
  return x / std::sqrt(arma::dot(x, x));
}
