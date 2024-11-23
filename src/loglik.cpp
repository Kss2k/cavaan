#include "model.h"
#include <float.h>


double getLogLikModel(arma::vec theta, Model *model) {
  fillModel(model, theta, false, true);

  double logLik = 0, valLDSigma, signLDSigma, valLDS, signLDS;
  bool okLDSigma, okLDS, okSigmaInv;
  arma::mat SigmaInv, S, Sigma;
  int p;

  MatricesGroup* matrices;

  for (int i = 0; i < model->ngroups; i++) {
    matrices = model->matricesGroups[i];

    S     = matrices->S[0];
    Sigma = matrices->Sigma[0];
    p     = matrices->p;

    okLDS      = arma::log_det(valLDS, signLDS, S);
    okLDSigma  = arma::log_det(valLDSigma, signLDSigma, Sigma);
    okSigmaInv = arma::inv(SigmaInv, Sigma);

    if (!okLDS || !okLDSigma || !okSigmaInv) return(DBL_MAX); // Not sure if this is the best approach...
    
    logLik += valLDSigma + arma::trace(S * SigmaInv) - valLDS - p;
  }

  return logLik;
}
