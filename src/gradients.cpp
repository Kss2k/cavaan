#include "cavaan.h"
#include <cmath>


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


arma::vec getGradientModelOVMeanStructure(const arma::vec& theta, Model* model) {
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
    Rcpp::Rcout << "DerivTau:" << "\n";    
    Rcpp::Rcout << DerivTau << "\n";    
    Rcpp::Rcout << "Nu:" << "\n";    
    Rcpp::Rcout << Nu << "\n";    
    Rcpp::Rcout << "Mu:" << "\n";    
    Rcpp::Rcout << Mu << "\n";    
    Rcpp::Rcout << "Tau" << "\n";    
    Rcpp::Rcout << Tau << "\n";    
    Rcpp::Rcout << "Mu calculated directly" << "\n";    
    Rcpp::Rcout << G * BStarInv * GammaStar * Tau << "\n";
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


arma::vec getGradientModelLVMeanStructure(const arma::vec& theta, Model* model) {
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
    const arma::vec& Tau       = matricesGroup->Tau;
    const arma::vec& Nu        = matricesGroup->Nu;
    const arma::vec& Mu        = matricesGroup->Mu;

    arma::mat I = arma::eye(Phi.n_rows, Phi.n_cols);
    arma::mat SigmaInv;
    if (!arma::inv(SigmaInv, Sigma)) Rcpp::stop("Sigma not invertible");

    const arma::mat BStarInv_T  = BStarInv.t();
    const arma::mat GammaStar_T = GammaStar.t();
    const arma::mat G_T         = G.t();
    const arma::mat Q           = SigmaInv - SigmaInv * S * SigmaInv;

    // Compute residuals and weighted residuals
    arma::vec delta = Nu - Mu;          // Residuals
    arma::vec m = SigmaInv * delta;     // Weighted residuals

    arma::mat A = G * BStarInv;         // A = G * B^{-1}

    // Compute derivative matrices without scaling factors
    arma::mat DerivPhi = (GammaStar_T * BStarInv_T * G_T * Q * G * BStarInv * GammaStar);

    // Derivative with respect to GammaStar
    arma::mat DerivGammaStar_cov = (BStarInv_T * G_T * Q * G * BStarInv * GammaStar * Phi);
    arma::mat DerivGammaStar_mean = - (A.t() * m) * Tau.t();
    arma::mat DerivGammaStar = DerivGammaStar_cov + DerivGammaStar_mean;

    // Derivative with respect to BStar
    arma::mat DerivBStar_cov = (BStarInv_T * G_T * Q * G * BStarInv * GammaStar * Phi * GammaStar_T * BStarInv_T);

    // Mean part calculations
    arma::vec C = GammaStar * Tau;       // C = GammaStar * Tau
    arma::vec S_vec = BStarInv * C;      // S = B^{-1} * C
    arma::mat V = G * BStarInv;          // V = G * B^{-1}
    arma::vec W = V.t() * m;             // W = V' * m

    arma::mat DerivBStar_mean = - S_vec * W.t();  // Outer product
    arma::mat DerivBStar = DerivBStar_cov + DerivBStar_mean;

    // Derivative with respect to Tau
    arma::vec DerivTau = - (A.t() * m);

    int t = 0;
    for (int i = 0; i < (int)(parTable->free.size()) && t < npar; i++) {
      if (!parTable->free[i]) continue;

      int row = parTable->row[i];
      int col = parTable->col[i];

      switch (parTable->matrix[i]) {
        case BETA_STAR:  {
          grad[t++] += -2 * DerivBStar.at(row, col);
          break;
        }
        case GAMMA_STAR: {
          grad[t++] += -2 * DerivGammaStar.at(row, col);
          break;
        }
        case PHI:        {
          // Adjust for symmetry
          if (row == col) {
            grad[t++] += 2 * DerivPhi.at(row, col);
          } else {
            grad[t++] += 4 * DerivPhi.at(row, col);
          }
          break;
        }
        case TAU:        {
          grad[t++] += -2 * DerivTau.at(row);
          break;
        }
        default: Rcpp::stop("Unrecognized matrix index");
      }
    }
  }

  return grad;
}


arma::vec getGradientModelGeneralNoMean(const arma::vec& theta, Model* model) {
  // THIS IS NOT IMPLEMENTED FOR MODELS WITH MEANSTRUCTURE!
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
      // THIS IS NOT IMPLEMENTED FOR MODELS WITH MEANSTRUCTURE!
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
