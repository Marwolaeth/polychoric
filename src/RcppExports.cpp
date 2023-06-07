// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// poly_tab
double poly_tab(Eigen::MatrixXd G, const double& correct);
RcppExport SEXP _polychoric_poly_tab(SEXP GSEXP, SEXP correctSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type G(GSEXP);
    Rcpp::traits::input_parameter< const double& >::type correct(correctSEXP);
    rcpp_result_gen = Rcpp::wrap(poly_tab(G, correct));
    return rcpp_result_gen;
END_RCPP
}
// poly_xy
double poly_xy(const Eigen::VectorXd& X, const Eigen::VectorXd& Y, const double& correct);
RcppExport SEXP _polychoric_poly_xy(SEXP XSEXP, SEXP YSEXP, SEXP correctSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const double& >::type correct(correctSEXP);
    rcpp_result_gen = Rcpp::wrap(poly_xy(X, Y, correct));
    return rcpp_result_gen;
END_RCPP
}
// poly_df
Eigen::MatrixXd poly_df(Rcpp::List X, double correct);
RcppExport SEXP _polychoric_poly_df(SEXP XSEXP, SEXP correctSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type correct(correctSEXP);
    rcpp_result_gen = Rcpp::wrap(poly_df(X, correct));
    return rcpp_result_gen;
END_RCPP
}
// poly_tab_full
Rcpp::List poly_tab_full(Eigen::MatrixXd G, const double& correct);
RcppExport SEXP _polychoric_poly_tab_full(SEXP GSEXP, SEXP correctSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Eigen::MatrixXd >::type G(GSEXP);
    Rcpp::traits::input_parameter< const double& >::type correct(correctSEXP);
    rcpp_result_gen = Rcpp::wrap(poly_tab_full(G, correct));
    return rcpp_result_gen;
END_RCPP
}
// poly_xy_full
Rcpp::List poly_xy_full(const Eigen::VectorXd& x, const Eigen::VectorXd& y, const double& correct);
RcppExport SEXP _polychoric_poly_xy_full(SEXP xSEXP, SEXP ySEXP, SEXP correctSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type y(ySEXP);
    Rcpp::traits::input_parameter< const double& >::type correct(correctSEXP);
    rcpp_result_gen = Rcpp::wrap(poly_xy_full(x, y, correct));
    return rcpp_result_gen;
END_RCPP
}
// poly_df_full
Rcpp::List poly_df_full(Rcpp::List X, double correct);
RcppExport SEXP _polychoric_poly_df_full(SEXP XSEXP, SEXP correctSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type X(XSEXP);
    Rcpp::traits::input_parameter< double >::type correct(correctSEXP);
    rcpp_result_gen = Rcpp::wrap(poly_df_full(X, correct));
    return rcpp_result_gen;
END_RCPP
}
// cor_polyserial
double cor_polyserial(const Eigen::VectorXd& x, const Eigen::VectorXd& d, const double& correct);
RcppExport SEXP _polychoric_cor_polyserial(SEXP xSEXP, SEXP dSEXP, SEXP correctSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const double& >::type correct(correctSEXP);
    rcpp_result_gen = Rcpp::wrap(cor_polyserial(x, d, correct));
    return rcpp_result_gen;
END_RCPP
}
// cor_polyserial_full
Rcpp::List cor_polyserial_full(const Eigen::VectorXd& x, const Eigen::VectorXd& d, const double& correct);
RcppExport SEXP _polychoric_cor_polyserial_full(SEXP xSEXP, SEXP dSEXP, SEXP correctSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const Eigen::VectorXd& >::type d(dSEXP);
    Rcpp::traits::input_parameter< const double& >::type correct(correctSEXP);
    rcpp_result_gen = Rcpp::wrap(cor_polyserial_full(x, d, correct));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_polychoric_poly_tab", (DL_FUNC) &_polychoric_poly_tab, 2},
    {"_polychoric_poly_xy", (DL_FUNC) &_polychoric_poly_xy, 3},
    {"_polychoric_poly_df", (DL_FUNC) &_polychoric_poly_df, 2},
    {"_polychoric_poly_tab_full", (DL_FUNC) &_polychoric_poly_tab_full, 2},
    {"_polychoric_poly_xy_full", (DL_FUNC) &_polychoric_poly_xy_full, 3},
    {"_polychoric_poly_df_full", (DL_FUNC) &_polychoric_poly_df_full, 2},
    {"_polychoric_cor_polyserial", (DL_FUNC) &_polychoric_cor_polyserial, 3},
    {"_polychoric_cor_polyserial_full", (DL_FUNC) &_polychoric_cor_polyserial_full, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_polychoric(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
