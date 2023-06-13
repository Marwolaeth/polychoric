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

// [[Rcpp::export]]
Eigen::MatrixXd correct_table_loop(Eigen::MatrixXd G, const double& correct = .1) {
  int r = G.rows();
  int s = G.cols();
  
  // Correct for continuity
  for (int i = 0; i < r; i++) {
    for (int j = 0; j < s; j++) {
      double x = G(i,j);
      double v = (x == 0.0) ? (x + correct) : x;
      G(i,j) = v;
    }
  }
  
  return G;
}

// [[Rcpp::export]]
Eigen::MatrixXd correct_table_mask(Eigen::MatrixXd G, const double& correct = .1) {
  Eigen::MatrixXd zero_entries = (G.array() == .0).cast<double>();
  
  G = G.array() + (G.array() == .0).cast<double>() * correct;
  
  return G;
}