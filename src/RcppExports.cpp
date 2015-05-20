// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// computeCoherency
ComplexVector computeCoherency(ComplexVector V, NumericVector d1, NumericVector d2);
RcppExport SEXP quantspec_computeCoherency(SEXP VSEXP, SEXP d1SEXP, SEXP d2SEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< ComplexVector >::type V(VSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d1(d1SEXP);
    Rcpp::traits::input_parameter< NumericVector >::type d2(d2SEXP);
    __result = Rcpp::wrap(computeCoherency(V, d1, d2));
    return __result;
END_RCPP
}
// computeSdNaive
ComplexVector computeSdNaive(ComplexVector V, NumericVector W);
RcppExport SEXP quantspec_computeSdNaive(SEXP VSEXP, SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< ComplexVector >::type V(VSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type W(WSEXP);
    __result = Rcpp::wrap(computeSdNaive(V, W));
    return __result;
END_RCPP
}
