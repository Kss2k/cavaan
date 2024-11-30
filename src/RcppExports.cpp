// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// gradLogLikCppNoMeans
arma::vec gradLogLikCppNoMeans(arma::vec theta, SEXP xptr);
RcppExport SEXP _cavaan_gradLogLikCppNoMeans(SEXP thetaSEXP, SEXP xptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type xptr(xptrSEXP);
    rcpp_result_gen = Rcpp::wrap(gradLogLikCppNoMeans(theta, xptr));
    return rcpp_result_gen;
END_RCPP
}
// gradLogLikCppOVMeans
arma::vec gradLogLikCppOVMeans(arma::vec theta, SEXP xptr);
RcppExport SEXP _cavaan_gradLogLikCppOVMeans(SEXP thetaSEXP, SEXP xptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type xptr(xptrSEXP);
    rcpp_result_gen = Rcpp::wrap(gradLogLikCppOVMeans(theta, xptr));
    return rcpp_result_gen;
END_RCPP
}
// gradLogLikCppLVMeans
arma::vec gradLogLikCppLVMeans(arma::vec theta, SEXP xptr);
RcppExport SEXP _cavaan_gradLogLikCppLVMeans(SEXP thetaSEXP, SEXP xptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type xptr(xptrSEXP);
    rcpp_result_gen = Rcpp::wrap(gradLogLikCppLVMeans(theta, xptr));
    return rcpp_result_gen;
END_RCPP
}
// gradLogLikNumericCpp
Rcpp::NumericVector gradLogLikNumericCpp(const arma::vec& theta, SEXP xptr, double h);
RcppExport SEXP _cavaan_gradLogLikNumericCpp(SEXP thetaSEXP, SEXP xptrSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type xptr(xptrSEXP);
    Rcpp::traits::input_parameter< double >::type h(hSEXP);
    rcpp_result_gen = Rcpp::wrap(gradLogLikNumericCpp(theta, xptr, h));
    return rcpp_result_gen;
END_RCPP
}
// ViewModelCreation
Rcpp::NumericVector ViewModelCreation(Rcpp::List RModel, arma::vec theta);
RcppExport SEXP _cavaan_ViewModelCreation(SEXP RModelSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type RModel(RModelSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(ViewModelCreation(RModel, theta));
    return rcpp_result_gen;
END_RCPP
}
// createRcppModel
RcppExport SEXP createRcppModel(Rcpp::List RModel);
RcppExport SEXP _cavaan_createRcppModel(SEXP RModelSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type RModel(RModelSEXP);
    rcpp_result_gen = Rcpp::wrap(createRcppModel(RModel));
    return rcpp_result_gen;
END_RCPP
}
// fillRcppModel
RcppExport SEXP fillRcppModel(SEXP xptr, arma::vec theta);
RcppExport SEXP _cavaan_fillRcppModel(SEXP xptrSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xptr(xptrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(fillRcppModel(xptr, theta));
    return rcpp_result_gen;
END_RCPP
}
// logLikCpp
Rcpp::NumericVector logLikCpp(const arma::vec& theta, SEXP xptr);
RcppExport SEXP _cavaan_logLikCpp(SEXP thetaSEXP, SEXP xptrSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< SEXP >::type xptr(xptrSEXP);
    rcpp_result_gen = Rcpp::wrap(logLikCpp(theta, xptr));
    return rcpp_result_gen;
END_RCPP
}
// debugCppModel
void debugCppModel(SEXP xptr, arma::vec theta);
RcppExport SEXP _cavaan_debugCppModel(SEXP xptrSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type xptr(xptrSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    debugCppModel(xptr, theta);
    return R_NilValue;
END_RCPP
}
// getVariablesEquation
Rcpp::StringVector getVariablesEquation(std::string expr);
RcppExport SEXP _cavaan_getVariablesEquation(SEXP exprSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type expr(exprSEXP);
    rcpp_result_gen = Rcpp::wrap(getVariablesEquation(expr));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_cavaan_gradLogLikCppNoMeans", (DL_FUNC) &_cavaan_gradLogLikCppNoMeans, 2},
    {"_cavaan_gradLogLikCppOVMeans", (DL_FUNC) &_cavaan_gradLogLikCppOVMeans, 2},
    {"_cavaan_gradLogLikCppLVMeans", (DL_FUNC) &_cavaan_gradLogLikCppLVMeans, 2},
    {"_cavaan_gradLogLikNumericCpp", (DL_FUNC) &_cavaan_gradLogLikNumericCpp, 3},
    {"_cavaan_ViewModelCreation", (DL_FUNC) &_cavaan_ViewModelCreation, 2},
    {"_cavaan_createRcppModel", (DL_FUNC) &_cavaan_createRcppModel, 1},
    {"_cavaan_fillRcppModel", (DL_FUNC) &_cavaan_fillRcppModel, 2},
    {"_cavaan_logLikCpp", (DL_FUNC) &_cavaan_logLikCpp, 2},
    {"_cavaan_debugCppModel", (DL_FUNC) &_cavaan_debugCppModel, 2},
    {"_cavaan_getVariablesEquation", (DL_FUNC) &_cavaan_getVariablesEquation, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_cavaan(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
