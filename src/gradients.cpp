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


arma::mat getGradientBStar(const arma::mat& BStar,
                           const arma::mat& G,
                           const arma::mat& GammaStar,
                           const arma::vec& Nu,
                           const arma::mat& Phi,
                           const arma::mat& S,
                           const arma::vec& Tau,
                           double          p) {

    const arma::mat T0 = arma::inv(BStar);                                    // B*⁻¹
    const arma::mat T1 = G * T0;
    const arma::mat T2 = T1 * GammaStar * Phi * GammaStar.t() * T0.t() * G.t();
    const arma::mat T3 = arma::inv(T2);                                       // T2⁻¹

    const arma::vec t4 = T0 * (GammaStar * Tau);
    const arma::vec t5 = Nu - G * t4;
    const arma::mat T6 = T1.t() * T3;
    const arma::vec t7 = T3 * t5;
    const arma::vec t8 = T0.t() * (G.t() * t7);

    double logDetT2,  signT2;
    double logDetS,   signS;
    if (!arma::log_det(logDetT2, signT2, T2))               // works for any nonsingular
      Rcpp::stop("log_det(T2) failed!");
    if (!arma::log_det(logDetS,  signS,  S))                // (S expected SPD)
      Rcpp::stop("log_det(S) failed!");

    // const arma::double functionValue =
    //     (logDetT2 + arma::trace(S * T3) - logDetS - p)
    //     + arma::as_scalar(t5.t() * t7);                // dot via *

    const arma::mat gradCore =
          T6 * S * T3 * G * T0 * GammaStar * Phi * GammaStar.t() * T0.t()
        - T6 *     G * T0 * GammaStar * Phi * GammaStar.t() * T0.t();

    const arma::vec inner =
        (t5.t() * T3 * G * T0 * GammaStar * Phi * GammaStar.t() * T0.t()).t();  // column vec

    const arma::mat gradient =
        2.0 * (gradCore + t8 * t4.t() + t8 * inner.t());

    return gradient;
}


arma::mat getGradientGammaStar(const arma::mat& BStar,
                                const arma::mat& G,
                                const arma::mat& GammaStar,
                                const arma::vec& Nu,
                                const arma::mat& Phi,
                                const arma::mat& S,
                                const arma::vec& Tau,
                                double          p) {

    arma::mat T0 = arma::inv(BStar);
    arma::mat T1 = G * T0;
    arma::mat T2 = T1 * GammaStar * Phi * GammaStar.t() * T0.t() * G.t();
    arma::mat T3 = arma::inv(T2);

    arma::vec t4 = Nu - G * ( T0 * ( GammaStar * Tau ) );          // (r)
    arma::mat T5 = T1.t() * T3;                                    // (r×r)
    arma::vec t6 = T3 * t4;                                        // (r)
    arma::vec t7 = T0.t() * ( G.t() * t6 );                        // (r)

    double logDetT2, signT2;
    if (!arma::log_det(logDetT2, signT2, T2)) 
      Rcpp::stop("log_det(T2) failed!");
    double logDetS,  signS;
    if (!arma::log_det(logDetS,  signS,  S))
      Rcpp::stop("log_det(S) failed!");

    // double functionValue =
    //     (logDetT2 + arma::trace(S * T3) - logDetS - p)
    //     + arma::as_scalar(t4.t() * t6);

    arma::mat gradLeft  = ( (T5 * G) * T0 ) * GammaStar * Phi;               // r×k
    arma::mat gradRight = ( ( (T5 * S) * T3 * G ) * T0 ) * GammaStar * Phi;  // r×k

    arma::rowvec inner = ( ( ( t4.t() * T3 ) * G ) * T0 * GammaStar * Phi ); // 1×k

    arma::mat gradient =
        2.0 * ( gradLeft - gradRight
              - t7 * Tau.t()              // outer(t7, Tau)
              - t7 * inner );             // outer(t7, inner)

    return gradient;
}

arma::mat getGradientPhi( const arma::mat& BStar,
                          const arma::mat& G,
                          const arma::mat& GammaStar,
                          const arma::vec& Nu,
                          const arma::mat& Phi,
                          const arma::mat& S,
                          const arma::vec& Tau,
                          double          p)
                    {
    /* ---------- core algebra (names exactly follow Python) ------------- */
    arma::mat T_0 = arma::inv(BStar);

    arma::mat T_1 = (G * T_0) * GammaStar;                                    //  r×k
    arma::mat T_2 = (((T_1 * Phi) * GammaStar.t()) * T_0.t()) * G.t();        //  r×r
    arma::mat T_3 = arma::inv(T_2);

    arma::vec t_4 = Nu - G * ( T_0 * ( GammaStar * Tau ) );                   //  r

    arma::mat T_5 = G * ( T_0 * GammaStar );                                  //  r×k
    arma::mat T_6 = arma::inv( (((T_1 * Phi.t()) * GammaStar.t()) * T_0.t()) * G.t() ); // r×r

    arma::mat T_7 = T_5.t() * T_6;                                            //  k×r
    arma::mat T_8 = T_5.t() * T_3;                                            //  k×r
    arma::vec t_9 = T_3 * t_4;                                                //  r

    double logDetT2, signT2;
    if (!arma::log_det(logDetT2, signT2, T_2))
        Rcpp::stop("arma::log_det failed for T_2 (singular or non-finite).");

    double logDetS,  signS;
    if (!arma::log_det(logDetS,  signS,  S))
        Rcpp::stop("arma::log_det failed for S (singular or non-finite).");

    // double functionValue =
    //       ( logDetT2 + arma::trace(S * T_3) - logDetS - p )
    //     + arma::as_scalar(t_4.t() * t_9);

    /* ---------- gradient (k × k) --------------------------------------- */
    arma::mat termA = ( (T_7 * G) * T_0 ) * GammaStar;                        // k×k
    arma::mat termB = ( (((T_7 * S.t()) * T_6) * G) * T_0 ) * GammaStar;      // k×k

    arma::vec  v1 = GammaStar.t() * ( T_0.t() * ( G.t() * ( T_6 * t_4 ) ) );  // k
    arma::rowvec w1 = ( ( ( t_4.t() * T_6 ) * G ) * T_0 ) * GammaStar;        // 1×k
    arma::mat outer1 = v1 * w1;                                               // k×k

    arma::mat termC = ( (T_8 * G) * T_0 ) * GammaStar;                        // k×k
    arma::mat termD = ( (((T_8 * S) * T_3) * G) * T_0 ) * GammaStar;          // k×k

    arma::vec  v2 = GammaStar.t() * ( T_0.t() * ( G.t() * t_9 ) );            // k
    arma::rowvec w2 = ( ( ( t_4.t() * T_3 ) * G ) * T_0 ) * GammaStar;        // 1×k
    arma::mat outer2 = v2 * w2;                                               // k×k

    arma::mat gradient = 0.5 *
        (  termA
         - termB
         - outer1
         + termC
         - termD
         - outer2 );

    return gradient;
}


arma::vec getGradientTau(const arma::mat& BStar,
                         const arma::mat& G,
                         const arma::mat& GammaStar,
                         const arma::vec& Nu,
                         const arma::mat& Phi,
                         const arma::mat& S,
                         const arma::vec& Tau,
                         double          p) {

    arma::mat T0 = arma::inv(BStar);                                         // B*⁻¹
    arma::mat T1 = G * T0 * GammaStar * Phi * GammaStar.t() * T0.t() * G.t();
    arma::mat T2 = arma::inv(T1);

    arma::vec t3 = Nu - G * (T0 * (GammaStar * Tau));                       // (r)
    arma::vec t4 = T2 * t3;                                                 // (r)

    double logDetT1 = 0.0, signT1 = 0.0;
    if (!arma::log_det(logDetT1, signT1, T1))
        Rcpp::stop("arma::log_det failed for T1 – matrix not invertible?");

    double logDetS = 0.0,  signS  = 0.0;
    if (!arma::log_det(logDetS,  signS,  S))
        Rcpp::stop("arma::log_det failed for S – matrix not SPD?");

    double functionValue =
        (logDetT1 + arma::trace(S * T2) - logDetS - p)
        + arma::as_scalar(t3.t() * t4);

    arma::vec gradient = -2.0 * (GammaStar.t() * (T0.t() * (G.t() * t4)));

    return gradient;
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
    const double p = matricesGroup->p;
 
    const arma::mat I = arma::eye(Phi.n_rows, Phi.n_cols);   
    arma::mat SigmaInv;
    if (!arma::inv(SigmaInv, Sigma)) Rcpp::stop("Sigma not invertible");

    const arma::mat DerivBStar = getGradientBStar(BStarInv, G, GammaStar, Nu, Phi, S, Tau, p);
    const arma::mat DerivGammaStar = getGradientGammaStar(BStarInv, G, GammaStar, Nu, Phi, S, Tau, p);
    const arma::mat DerivPhi = getGradientPhi(BStarInv, G, GammaStar, Nu, Phi, S, Tau, p);
    const arma::mat DerivTau = getGradientTau(BStarInv, G, GammaStar, Nu, Phi, S, Tau, p);
    
    int t = 0;
    for (int i = 0; i < (int)(parTable->free.size()) && t < npar; i++) {
      if (!parTable->free[i]) continue;

      int row = parTable->row[i];
      int col = parTable->col[i];

      switch (parTable->matrix[i]) {
        case BETA_STAR:  {grad[t++] += DerivBStar.at(row, col); break;}
        case GAMMA_STAR: {grad[t++] += DerivGammaStar.at(row, col); break;}
        // case PHI:        {grad[t++] += (2 - I.at(row, col)) * DerivPhi.at(row, col); break;}
        case PHI:        {grad[t++] += DerivPhi.at(row, col); break;}
        case TAU:        {grad[t++] += DerivTau.at(row, 0); break;}
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
