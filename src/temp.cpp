#include <Rcpp.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]

#include <cmath>
#include <numeric>
#include <vector>
#include <unordered_set>

const double pi = 3.1415926535; // set pi

Eigen::Array<bool, Eigen::Dynamic, 1> is_na(const Eigen::VectorXd& x) {
  return x.array().isNaN();
}

// Find pairwise complete observations for two vectors
Eigen::Array<bool, Eigen::Dynamic, 1> pairwise_complete(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& y
) {
  Eigen::Array<bool, Eigen::Dynamic, 1> miss_x = x.array().isNaN();
  Eigen::Array<bool, Eigen::Dynamic, 1> miss_y = y.array().isNaN();
  return !(miss_x || miss_y);
}

// Filter vector based on the boolean mask
Eigen::VectorXd filter(
    const Eigen::VectorXd& x,
    const Eigen::Array<bool, Eigen::Dynamic, 1>& mask
) {
  // Vreate empty vector to store selected values
  Eigen::VectorXd xc = {};
  
  // Filter the array based on the mask
  for (int i = 0; i < x.size(); i++) {
    // Analogous to std::push_back()
    if (mask(i)) {
      Eigen::VectorXd tmp = xc;
      xc.resize(xc.size()+1);
      xc.head(tmp.size()) = tmp;
      xc(xc.size()-1) = x(i);
    }
  }
  
  return xc;
}

// Function to count number of unique values in a vector
int n_unique(const Eigen::VectorXd& X) {
  // Filter out NA values
  Eigen::VectorXd x = filter(X, !is_na(X));
  
  int n = x.size();
  std::unordered_set<double> seen;
  for (int i = 0; i < n; i++) {
    seen.insert(x(i));
  }
  return seen.size();
}

// Function to calculate rank of vector elements with ties resolved by mean rank
Eigen::VectorXd rank_vec(const Eigen::VectorXd& v) {
  // Filter out NA values
  Eigen::VectorXd x = filter(v, !is_na(v));
  int n = x.size();
  
  std::vector<std::size_t> w(n);
  std::iota(begin(w), end(w), 0);
  std::sort(begin(w), end(w), 
            [&v](std::size_t i, std::size_t j) { return v(i) < v(j); });
  
  Eigen::VectorXd r(w.size());
  for (std::size_t n, i = 0; i < w.size(); i += n)
  {
    n = 1;
    while (i + n < w.size() && v(w[i]) == v(w[i+n])) ++n;
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

double cor_spearman(const Eigen::VectorXd& X, const Eigen::VectorXd& Y) {
  // Filter for pairwise complete observations
  Eigen::Array<bool, Eigen::Dynamic, 1> pco = pairwise_complete(X, Y);
  Eigen::VectorXd x = filter(X, pco);
  Eigen::VectorXd y = filter(Y, pco);
  
  long long int n = x.size();
  Eigen::VectorXd d(n);
  double rho;
  // check for equal length
  // if (y.size() != n) {
  //   Rcpp::stop("Vectors are of different dimensionality");
  // }
  d = rank_vec(x).array() - rank_vec(y).array();
  rho = 1.0 - (6.0 * (d.dot(d))) / (n*(n*n - 1));
  rho = (rho < -1.0) ? -1.0 : rho;
  rho = (rho >  1.0) ?  1.0 : rho;
  return rho;
}

Eigen::VectorXd correct_data(Eigen::VectorXd x) {
  double minx = x.minCoeff();
  if (minx <= 0) {
    Rcpp::warning("Variable is scaled to have a minimum of at least 1");
    double correct = 1.0 - minx;
    x = x.array() + correct;
  }
  return x;
}

double Phi(const double& x, const double& e = 0.0, const double& s = 1.0) {
  double z = (x - e) / (s * std::sqrt(2)); // calculate z-score
  double cdf = 0.5 * (1 + std::erf(z)); // calculate CDF using erf function
  return cdf; // return CDF value
}

double Phi_inv(const double& p, const double& e = 0.0, const double& s = 1.0) {
  // check for valid input probability
  // if (p < 0.0 || p > 1.0) {
  //   Rcpp::stop("Invalid input probability");
  // }
  // initial guess for inverse CDF value
  double z = 0.0;
  // use Newton's method to iteratively solve for inverse CDF value
  while (true) {
    double cdf = Phi(z);
    double pdf = std::exp(-0.5 * z * z) / std::sqrt(2 * M_PI);
    double z_new = z - (cdf - p) / pdf;
    if (std::abs(z_new - z) < 1e-8) {
      break;
    }
    z = z_new;
  }
  double x = z*s + e; 
  return x;
}

Eigen::VectorXd Phi_inv_vec(
    const Eigen::VectorXd& p,
    const double& e = 0.0,
    const double& s = 1.0
) {
  // check for valid input standard deviations
  // if (s <= 0) {
  //   Rcpp::stop("Standard deviation must be positive");
  // }
  int n = p.size();
  Eigen::VectorXd x(n);
  for (int i = 0; i < n; i++) {
    x(i) = Phi_inv(p(i), e, s);
  }
  return x;
}

// normal probability density function of x given e and s
// [[Rcpp::export]]
double phi(const double& x, const double& e = 0.0, const double& s = 1.0) {
  return std::exp(-0.5 * std::pow((x - e) / s, 2)) / (s * std::sqrt(2.0 * pi));
}

// # UNIVARIATE ----
// [[Rcpp::export]]
double Var(const Eigen::VectorXd& x) {
  Eigen::ArrayXd x_arr = x.array().isNaN().select(0, x);
  Eigen::Index n = x_arr.size();
  double var = (x_arr - x_arr.mean()).square().sum() / (n - 1);
  return var;
}

double SD(const Eigen::VectorXd& x) {
  double v = Var(x);
  return std::sqrt(v);
}

double skewness(const Eigen::VectorXd& x) {
  Eigen::ArrayXd x_arr = x.array();
  Eigen::Index n = x.size();
  double mean = x.mean();
  double var = (x_arr - mean).square().sum() / (n - 1);
  double std_dev = std::sqrt(var);
  double skew = ((x_arr - mean).pow(3.0)).sum() /
    ((n - 1) * std_dev * std_dev * std_dev);
  return skew;
}

double Cov(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
  Eigen::Index n = x.size();
  Eigen::VectorXd x_centered = x.array() - x.mean();
  Eigen::VectorXd y_centered = y.array() - y.mean();
  double covar = (x_centered.dot(y_centered)) / (n - 1);
  return covar;
}

// [[Rcpp::export]]
double cor_pearson(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
  return Cov(x, y) / (SD(x) * SD(y));
}

Eigen::VectorXd cumsum(const Eigen::VectorXd& x) {
  Eigen::VectorXd y = Eigen::VectorXd::Zero(x.size());
  for (int i = 0; i < x.size(); i++) {
    y(i) = (i > 0) ? y(i - 1) + x(i) : x(i);
  }
  return y;
}

// Function to estimate thresholds from an Eigen VectorXd object
Eigen::VectorXd estimate_thresholds(
    const Eigen::VectorXd& X,
    const double& correct = 1e-08
) {
  // Filter out NA values
  Eigen::VectorXd x = filter(X, !is_na(X));
  int N = x.size();
  
  int m = x.maxCoeff();
  Eigen::VectorXd counts = Eigen::VectorXd::Zero(m);
  for (int i = 0; i < N; i++) {
    counts(x(i)-1) += 1.0;
  }
  for (int i = 0; i < m; i++) {
    double cnt = counts(i);
    counts(i) = (cnt == 0.0) ? correct : cnt;
  }
  for (int i = 1; i < m; i++) {
    double cnt = counts(i);
    counts(i) += counts(i-1);
  }
  Eigen::VectorXd p(m-1);
  p = counts.head(m-1).array() / N;
  Eigen::VectorXd q(m-1);
  q = Phi_inv_vec(p);
  return q;
}

// [[Rcpp::export]]
double polyserial_(const Eigen::VectorXd& x, const Eigen::VectorXi& y) {
  // …
}