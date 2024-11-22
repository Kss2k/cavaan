#include "model.h"


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


arma::vec getGradientModel(arma::vec theta, Model *model) {
  fillModel(model, theta, false, true);
  
  int npar = model->parTable->nfree, ngroups = model->ngroups;

  GradientMatricesParam *gradientMatricesParam;
  MatricesGroup *gradientMatricesParamGroup;
 
  arma::mat DerivSigma, DerivGammaStar, DerivBStar, DerivPhi, 
    DerivBStarInv, G, Phi, GammaStar, BStar, BStarInv, Sigma, SigmaInv, S;

  arma::vec grad = arma::vec(npar).fill(0);

  std::vector<arma::mat> SigmaInvs;
  for (int g = 0; g < ngroups; g++) {
    SigmaInvs.push_back(arma::inv(model->matricesGroups[g]->Sigma[0]));
  }

  for (int t = 0; t < npar; t++) {
    gradientMatricesParam = model->gradientMatricesParams[t];

    for (int g = 0; g < ngroups; g++) {
      gradientMatricesParamGroup = gradientMatricesParam->matricesGroups[g];
      
      DerivGammaStar = gradientMatricesParamGroup->GammaStar[0];
      DerivBStar     = gradientMatricesParamGroup->BStar[0];
      DerivPhi       = gradientMatricesParamGroup->Phi[0];

      G         = model->matricesGroups[g]->G[0];
      S         = model->matricesGroups[g]->S[0];
      Phi       = model->matricesGroups[g]->Phi[0];
      GammaStar = model->matricesGroups[g]->GammaStar[0];
      BStar     = model->matricesGroups[g]->BStar[0];
      BStarInv  = model->matricesGroups[g]->BStarInv[0];
      Sigma     = model->matricesGroups[g]->Sigma[0];

      SigmaInv  = SigmaInvs[g];

      DerivBStarInv = -BStarInv * DerivBStar * BStarInv;

      DerivSigma = G * (
        DerivBStarInv * GammaStar      * Phi      * GammaStar.t()      * BStarInv.t() +
        BStarInv      * DerivGammaStar * Phi      * GammaStar.t()      * BStarInv.t() +
        BStarInv      * GammaStar      * DerivPhi * GammaStar.t()      * BStarInv.t() +
        BStarInv      * GammaStar      * Phi      * DerivGammaStar.t() * BStarInv.t() +
        BStarInv      * GammaStar      * Phi      * GammaStar.t()      * DerivBStarInv.t()
          ) * G.t();
     
      grad[t] += arma::trace((SigmaInv - SigmaInv * S * SigmaInv) * DerivSigma);
    }
  }

  return grad;
}
