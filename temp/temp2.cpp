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

// Function to calculate rank of vector elements with ties resolved by mean rank
// [[Rcpp::export]]
Eigen::VectorXd rank_(const Eigen::VectorXd& x) {
  int n = x.size();
  
  std::vector<std::size_t> w(n);
  std::iota(begin(w), end(w), 0);
  std::sort(begin(w), end(w), 
            [&x](std::size_t i, std::size_t j) { return x(i) < x(j); });
  
  Eigen::VectorXd r(w.size());
  for (std::size_t n, i = 0; i < w.size(); i += n)
  {
    n = 1;
    while (i + n < w.size() && x(w[i]) == x(w[i+n])) ++n;
    for (std::size_t k = 0; k < n; ++k)
    {
      r(w[i+k]) = i + (n + 1) / 2.0; // average rank of n tied values
      // r(w[i+k]) = i + 1;          // min 
      // r(w[i+k]) = i + n;          // max
      // r(w[i+k]) = i + k + 1;      // random order
    }
  }
  return r;
}