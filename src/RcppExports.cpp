// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// ibm
Rcpp::List ibm(const int nbInd, const int length, const double prob_mut, const int loci, const Eigen::MatrixXf& mutationMatrix, const Eigen::MatrixXf& selectionMatrix);
RcppExport SEXP _IBMselection_ibm(SEXP nbIndSEXP, SEXP lengthSEXP, SEXP prob_mutSEXP, SEXP lociSEXP, SEXP mutationMatrixSEXP, SEXP selectionMatrixSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int >::type nbInd(nbIndSEXP);
    Rcpp::traits::input_parameter< const int >::type length(lengthSEXP);
    Rcpp::traits::input_parameter< const double >::type prob_mut(prob_mutSEXP);
    Rcpp::traits::input_parameter< const int >::type loci(lociSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXf& >::type mutationMatrix(mutationMatrixSEXP);
    Rcpp::traits::input_parameter< const Eigen::MatrixXf& >::type selectionMatrix(selectionMatrixSEXP);
    rcpp_result_gen = Rcpp::wrap(ibm(nbInd, length, prob_mut, loci, mutationMatrix, selectionMatrix));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_IBMselection_ibm", (DL_FUNC) &_IBMselection_ibm, 6},
    {NULL, NULL, 0}
};

RcppExport void R_init_IBMselection(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}