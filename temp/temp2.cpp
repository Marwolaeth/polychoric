// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <cmath>
#include <numeric>
#include <vector>
#include <unordered_set>
// #include <chrono>
// #include <ctime>
#include "toms462.cpp"
#include "LBFGSB.h"

// Logical vectors and matrices — very useful
typedef Eigen::Array<bool, Eigen::Dynamic, 1> VectorXl;
typedef Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> MatrixXl;

// ## UTILS
// [[Rcpp::export]]
VectorXl is_na(const Eigen::VectorXd& x) {
  return x.array().isNaN();
}

// [[Rcpp::export]]
MatrixXl shadow_matrix(
    const Eigen::MatrixXd& M
) {
  MatrixXl shadow = !(M.array().isNaN());
  
  return shadow.matrix();
}