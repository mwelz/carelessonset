// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// LSP
Rcpp::IntegerVector LSP(Rcpp::IntegerVector x, int J);
RcppExport SEXP _carelessonset_LSP(SEXP xSEXP, SEXP JSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type J(JSEXP);
    rcpp_result_gen = Rcpp::wrap(LSP(x, J));
    return rcpp_result_gen;
END_RCPP
}
// Tn_multivariate
Rcpp::NumericVector Tn_multivariate(Rcpp::NumericMatrix x, Rcpp::StringVector theta);
RcppExport SEXP _carelessonset_Tn_multivariate(SEXP xSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(Tn_multivariate(x, theta));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _carelessonset_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// mean_ab_cpp
double mean_ab_cpp(std::vector<double> x, int start, int stop);
RcppExport SEXP _carelessonset_mean_ab_cpp(SEXP xSEXP, SEXP startSEXP, SEXP stopSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<double> >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type start(startSEXP);
    Rcpp::traits::input_parameter< int >::type stop(stopSEXP);
    rcpp_result_gen = Rcpp::wrap(mean_ab_cpp(x, start, stop));
    return rcpp_result_gen;
END_RCPP
}
// Tn_univariate
Rcpp::NumericVector Tn_univariate(Rcpp::NumericVector x, Rcpp::StringVector theta);
RcppExport SEXP _carelessonset_Tn_univariate(SEXP xSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    Rcpp::traits::input_parameter< Rcpp::StringVector >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(Tn_univariate(x, theta));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_carelessonset_LSP", (DL_FUNC) &_carelessonset_LSP, 2},
    {"_carelessonset_Tn_multivariate", (DL_FUNC) &_carelessonset_Tn_multivariate, 2},
    {"_carelessonset_rcpp_hello_world", (DL_FUNC) &_carelessonset_rcpp_hello_world, 0},
    {"_carelessonset_mean_ab_cpp", (DL_FUNC) &_carelessonset_mean_ab_cpp, 3},
    {"_carelessonset_Tn_univariate", (DL_FUNC) &_carelessonset_Tn_univariate, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_carelessonset(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
