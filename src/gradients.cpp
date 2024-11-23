#include "model.h"
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


arma::vec getGradientModel(arma::vec theta, Model* model) {
  fillModel(model, theta, false, true);

  int npar = model->parTable->nfree, ngroups = model->ngroups;

  arma::vec grad = arma::vec(npar).fill(0);

  // Precompute Sigma inverses for all groups
  std::vector<arma::mat> SigmaInvs(ngroups);
  for (int g = 0; g < ngroups; g++) {
    arma::mat Sigma = model->matricesGroups[g]->Sigma;
    if (!arma::inv(SigmaInvs[g], Sigma)) {
      Rcpp::stop("Sigma not invertible");
    }
  }

  // Iterate over free parameters
  for (int t = 0; t < npar; t++) {
    GradientMatricesParam* gradientMatricesParam = model->gradientMatricesParams[t];

    // Iterate over groups
    for (int g = 0; g < ngroups; g++) {
      MatricesGroup* gradientMatricesParamGroup = gradientMatricesParam->matricesGroups[g];
      MatricesGroup* matricesGroup = model->matricesGroups[g];

      // Load derivative matrices
      arma::mat DerivGammaStar = gradientMatricesParamGroup->GammaStar;
      arma::mat DerivBStar = gradientMatricesParamGroup->BStar;
      arma::mat DerivPhi = gradientMatricesParamGroup->Phi;

      // Load base matrices
      arma::mat G = matricesGroup->G;
      arma::mat S = matricesGroup->S;
      arma::mat Phi = matricesGroup->Phi;
      arma::mat GammaStar = matricesGroup->GammaStar;
      arma::mat BStar = matricesGroup->BStar;
      arma::mat BStarInv = matricesGroup->BStarInv;
      arma::mat Sigma = matricesGroup->Sigma;
      arma::mat SigmaInv = SigmaInvs[g];

      // Precompute reusable terms
      arma::mat BStarInv_T = BStarInv.t();
      arma::mat GammaStar_T = GammaStar.t();

      // Derivative of BStarInv
      arma::mat DerivBStarInv = -BStarInv * DerivBStar * BStarInv;

      // Derivative of Sigma
      arma::mat DerivSigma = G * (
          DerivBStarInv * GammaStar * Phi * GammaStar_T * BStarInv_T +
          BStarInv * DerivGammaStar * Phi * GammaStar_T * BStarInv_T +
          BStarInv * GammaStar * DerivPhi * GammaStar_T * BStarInv_T +
          BStarInv * GammaStar * Phi * DerivGammaStar.t() * BStarInv_T +
          BStarInv * GammaStar * Phi * GammaStar_T * DerivBStarInv.t()
      ) * G.t();

      // Update gradient
      arma::mat diff = SigmaInv - SigmaInv * S * SigmaInv;
      grad[t] += arma::trace(diff * DerivSigma);
    }
  }

  return grad;
}


arma::vec normalize(arma::vec x) {
  return x / std::sqrt(arma::dot(x, x));
}
