// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// IsingProcess
IntegerMatrix IsingProcess(int nSample, NumericMatrix graph, NumericVector thresholds, double beta, IntegerVector responses);
RcppExport SEXP IsingSampler_IsingProcess(SEXP nSampleSEXP, SEXP graphSEXP, SEXP thresholdsSEXP, SEXP betaSEXP, SEXP responsesSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type nSample(nSampleSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type graph(graphSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thresholds(thresholdsSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type responses(responsesSEXP);
    __result = Rcpp::wrap(IsingProcess(nSample, graph, thresholds, beta, responses));
    return __result;
END_RCPP
}
// IsingSamplerCpp
IntegerMatrix IsingSamplerCpp(int n, NumericMatrix graph, NumericVector thresholds, double beta, int nIter, IntegerVector responses, bool exact, IntegerMatrix constrain);
RcppExport SEXP IsingSampler_IsingSamplerCpp(SEXP nSEXP, SEXP graphSEXP, SEXP thresholdsSEXP, SEXP betaSEXP, SEXP nIterSEXP, SEXP responsesSEXP, SEXP exactSEXP, SEXP constrainSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type graph(graphSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thresholds(thresholdsSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< int >::type nIter(nIterSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type responses(responsesSEXP);
    Rcpp::traits::input_parameter< bool >::type exact(exactSEXP);
    Rcpp::traits::input_parameter< IntegerMatrix >::type constrain(constrainSEXP);
    __result = Rcpp::wrap(IsingSamplerCpp(n, graph, thresholds, beta, nIter, responses, exact, constrain));
    return __result;
END_RCPP
}
// H
double H(NumericMatrix J, IntegerVector s, NumericVector h);
RcppExport SEXP IsingSampler_H(SEXP JSEXP, SEXP sSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type J(JSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h(hSEXP);
    __result = Rcpp::wrap(H(J, s, h));
    return __result;
END_RCPP
}
// f
double f(IntegerMatrix Y, NumericMatrix J, NumericVector h);
RcppExport SEXP IsingSampler_f(SEXP YSEXP, SEXP JSEXP, SEXP hSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerMatrix >::type Y(YSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type J(JSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type h(hSEXP);
    __result = Rcpp::wrap(f(Y, J, h));
    return __result;
END_RCPP
}
// Hvec
double Hvec(IntegerVector s, NumericVector Theta, int N);
RcppExport SEXP IsingSampler_Hvec(SEXP sSEXP, SEXP ThetaSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerVector >::type s(sSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Theta(ThetaSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    __result = Rcpp::wrap(Hvec(s, Theta, N));
    return __result;
END_RCPP
}
// expvalues
NumericVector expvalues(IntegerMatrix x);
RcppExport SEXP IsingSampler_expvalues(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerMatrix >::type x(xSEXP);
    __result = Rcpp::wrap(expvalues(x));
    return __result;
END_RCPP
}
// vec2Thresh
NumericVector vec2Thresh(NumericVector vec, int P);
RcppExport SEXP IsingSampler_vec2Thresh(SEXP vecSEXP, SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< int >::type P(PSEXP);
    __result = Rcpp::wrap(vec2Thresh(vec, P));
    return __result;
END_RCPP
}
// vec2Graph
NumericMatrix vec2Graph(NumericVector vec, int P);
RcppExport SEXP IsingSampler_vec2Graph(SEXP vecSEXP, SEXP PSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericVector >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< int >::type P(PSEXP);
    __result = Rcpp::wrap(vec2Graph(vec, P));
    return __result;
END_RCPP
}
// Broderick2013
NumericVector Broderick2013(IntegerMatrix x, int M, int T, int nIter, IntegerVector responses);
RcppExport SEXP IsingSampler_Broderick2013(SEXP xSEXP, SEXP MSEXP, SEXP TSEXP, SEXP nIterSEXP, SEXP responsesSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< IntegerMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type nIter(nIterSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type responses(responsesSEXP);
    __result = Rcpp::wrap(Broderick2013(x, M, T, nIter, responses));
    return __result;
END_RCPP
}
// PseudoLikelihood
double PseudoLikelihood(NumericMatrix x, NumericMatrix graph, NumericVector thresholds, double beta, IntegerVector responses, bool logis);
RcppExport SEXP IsingSampler_PseudoLikelihood(SEXP xSEXP, SEXP graphSEXP, SEXP thresholdsSEXP, SEXP betaSEXP, SEXP responsesSEXP, SEXP logisSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type graph(graphSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type thresholds(thresholdsSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type responses(responsesSEXP);
    Rcpp::traits::input_parameter< bool >::type logis(logisSEXP);
    __result = Rcpp::wrap(PseudoLikelihood(x, graph, thresholds, beta, responses, logis));
    return __result;
END_RCPP
}
