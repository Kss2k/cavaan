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


arma::vec getGradientModel(const arma::vec& theta, Model* model) {
    fillModel(model, theta, false, true);

    int npar = model->parTable->nfree;
    int ngroups = model->ngroups;

    arma::vec grad(npar, arma::fill::zeros);

    // Precompute Sigma inverses for all groups
    std::vector<arma::mat> SigmaInvs(ngroups);
    for (int g = 0; g < ngroups; ++g) {
        const arma::mat& Sigma = model->matricesGroups[g]->Sigma;
        if (!arma::inv(SigmaInvs[g], Sigma)) {
            Rcpp::stop("Sigma not invertible");
        }
    }

    // Iterate over free parameters
    for (int t = 0; t < npar; ++t) {
        GradientMatricesParam* gradientMatricesParam = model->gradientMatricesParams[t];

        // Iterate over groups
        for (int g = 0; g < ngroups; ++g) {
            MatricesGroup* gradientMatricesParamGroup = gradientMatricesParam->matricesGroups[g];
            MatricesGroup* matricesGroup = model->matricesGroups[g];

            // Load derivative matrices using references to avoid copies
            const arma::mat& DerivGammaStar = gradientMatricesParamGroup->GammaStar;
            const arma::mat& DerivBStar = gradientMatricesParamGroup->BStar;
            const arma::mat& DerivPhi = gradientMatricesParamGroup->Phi;

            // Load base matrices using references
            const arma::mat& G = matricesGroup->G;
            const arma::mat& S = matricesGroup->S;
            const arma::mat& Phi = matricesGroup->Phi;
            const arma::mat& GammaStar = matricesGroup->GammaStar;
            const arma::mat& BStar = matricesGroup->BStar;
            const arma::mat& BStarInv = matricesGroup->BStarInv;
            const arma::mat& SigmaInv = SigmaInvs[g];

            // Precompute reusable terms
            const arma::mat BStarInv_T = BStarInv.t();
            const arma::mat GammaStar_T = GammaStar.t();
            const arma::mat DerivGammaStar_T = DerivGammaStar.t();
            const arma::mat DerivBStarInv = -BStarInv * DerivBStar * BStarInv;

            // Compute Derivative of Sigma more efficiently
            arma::mat temp = BStarInv * GammaStar;
            arma::mat temp_T = temp.t();
            arma::mat DerivSigmaInner = (
                DerivBStarInv * GammaStar * Phi * GammaStar_T * BStarInv_T +
                BStarInv * DerivGammaStar * Phi * GammaStar_T * BStarInv_T +
                temp * DerivPhi * GammaStar_T * BStarInv_T +
                temp * Phi * DerivGammaStar_T * BStarInv_T +
                temp * Phi * GammaStar_T * DerivBStarInv.t()
            );
            arma::mat DerivSigma = G * DerivSigmaInner * G.t();

            // Update gradient using element-wise multiplication for efficiency
            arma::mat diff = SigmaInv - SigmaInv * S * SigmaInv;
            grad[t] += arma::accu(diff % DerivSigma);
        }
    }

    return grad;
}

arma::vec normalize(arma::vec x) {
  return x / std::sqrt(arma::dot(x, x));
}
