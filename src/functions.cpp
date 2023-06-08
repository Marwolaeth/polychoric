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

/*
 TOMS462 is a C++ library which evaluates the upper right tail of the bivariate
 normal distribution; that is, the probability that normal variables
 X and Y with correlation R will satisfy H <= X and K <= Y.
 
 Library provided by the University of South Carolina under the GNU LGPL license.
 See https://people.math.sc.edu/Burkardt/cpp_src/toms462/toms462.html
 
 Citations:
 1. Thomas Donnelly,
 Algorithm 462: Bivariate Normal Distribution,
 Communications of the ACM,
 October 1973, Volume 16, Number 10, page 638.
 2. Donald Owen,
 Tables for Computing Bivariate Normal Probabilities,
 Annals of Mathematical Statistics,
 Volume 27, Number 4, pages 1075-1090, December 1956.
 */

/*
 LBFGS++ is a header-only C++ library by Yixuan Qiu that implements the
 Limited-memory BFGS algorithm (L-BFGS) for unconstrained minimization problems,
 and a modified version of the L-BFGS-B algorithm for box-constrained ones.
 
 The code for the L-BFGS solver is derived and modified from the libLBFGS
 library developed by Naoaki Okazaki.
 
 LBFGS++ is implemented as a header-only C++ library, whose only dependency,
 Eigen, is also header-only.
 
 The development page of LBFGS++ is at https://github.com/yixuan/LBFGSpp.
 
 Library provided under the MIT license.
 See https://lbfgspp.statr.me/
 
 Citations:
 1. Thomas Donnelly,
 Algorithm 462: Bivariate Normal Distribution,
 Communications of the ACM,
 October 1973, Volume 16, Number 10, page 638.
 2. Donald Owen,
 Tables for Computing Bivariate Normal Probabilities,
 Annals of Mathematical Statistics,
 Volume 27, Number 4, pages 1075-1090, December 1956.
 */

const double pi = 3.1415926535; // set pi

// ## UTILS
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
  // Create empty vector to store selected values
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

Eigen::MatrixXd list_to_matrix(Rcpp::List lst) {
  // Determine the number of rows and columns
  int n_rows = Rcpp::as<Rcpp::NumericVector>(lst[0]).size();
  int n_cols = lst.size();
  
  // Create the Eigen::Matrix object
  Eigen::MatrixXd mat(n_rows, n_cols);
  
  // Loop over the columns and copy the data
  for (int j = 0; j < n_cols; j++) {
    Rcpp::NumericVector col = Rcpp::as<Rcpp::NumericVector>(lst[j]);
    for (int i = 0; i < n_rows; i++) {
      mat(i, j) = col[i];
    }
  }
  
  return mat;
}

// [[Rcpp::export]]
Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> shadow_matrix(
    const Eigen::MatrixXd& M
) {
  Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> shadow = !(M.array().isNaN());
  
  return shadow.matrix();
}

// Function to count number of unique values in a vector
int n_unique(
    const Eigen::VectorXd& X,
    Eigen::Array<bool, Eigen::Dynamic, 1> keep = {}
) {
  // Filter out NA values
  if (keep.size() == 0) {
    keep = !(X.array().isNaN());
  }
  Eigen::VectorXd x = filter(X, keep);
  
  int n = x.size();
  std::unordered_set<double> seen;
  for (int i = 0; i < n; i++) {
    seen.insert(x(i));
  }
  return seen.size();
}

// Function to calculate rank of vector elements with ties resolved by mean rank
Eigen::VectorXd rank_(
    const Eigen::VectorXd& v,
    Eigen::Array<bool, Eigen::Dynamic, 1> keep = {}
) {
  // Filter out NA values
  if (keep.size() == 0) {
    keep = !(v.array().isNaN());
  }
  Eigen::VectorXd x = filter(v, keep);
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

// Check ordinal data for zero value to avoid d(-1) reference
Eigen::VectorXd correct_data(Eigen::VectorXd x) {
  double minx = x.minCoeff();
  if (minx <= 0) {
    Rcpp::warning("Variable is scaled to have a minimum of 1\n");
    double correct = 1.0 - minx;
    x = x.array() + correct;
  }
  return x;
}

// ## DISTRIBUTIONS
// Normal probability density function of x given e and s
double normal_pdf(const double& x, const double& e = 0.0, const double& s = 1.0) {
  return std::exp(-0.5 * std::pow((x - e) / s, 2)) / (s * std::sqrt(2.0 * pi));
}

// Standard normal probability density function
double phi(const double& z) {
  return std::exp(-0.5 * std::pow(z, 2)) / (std::sqrt(2.0 * pi));
}

// Normal probability density function of a numerical vector X given e and s
Eigen::VectorXd normal_pdf_(const Eigen::VectorXd& x, const double& e = 0.0, const double& s = 1.0) {
  const int n = x.size();
  const double denom = s * std::sqrt(2.0 * pi);
  
  Eigen::VectorXd pdf(n);
  for (int i = 0; i < n; i++) {
    pdf(i) = std::exp(-0.5 * std::pow((x(i) - e) / s, 2)) / denom;
  }
  return pdf;
}

Eigen::VectorXd phi_(const Eigen::VectorXd& z) {
  const int n = z.size();
  const double denom = std::sqrt(2.0 * pi);
  
  Eigen::VectorXd pdf(n);
  for (int i = 0; i < n; i++) {
    pdf(i) = std::exp(-0.5 * std::pow(z(i), 2)) / denom;
  }
  return pdf;
}

// Normal distribution function
double normal_cdf(const double& x, const double& e = 0.0, const double& s = 1.0) {
  double z = (x - e) / (s * std::sqrt(2)); // calculate z-score
  double cdf = 0.5 * (1 + std::erf(z)); // calculate CDF using erf function
  return cdf; // return CDF value
}

double Phi(const double& x) {
  double z = x / std::sqrt(2);
  double cdf = 0.5 * (1 + std::erf(z)); // calculate CDF using erf function
  return cdf; // return CDF value
}

// Inverse normal distribution
double normal_quantile(const double& p, const double& e = 0.0, const double& s = 1.0) {
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

double Phi_inv(const double& p) {
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
  return z;
}

Eigen::VectorXd normal_quantile_(
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
    x(i) = normal_quantile(p(i), e, s);
  }
  return x;
}

Eigen::VectorXd Phi_inv_(const Eigen::VectorXd& p) {
  // check for valid input standard deviations
  // if (s <= 0) {
  //   Rcpp::stop("Standard deviation must be positive");
  // }
  int n = p.size();
  Eigen::VectorXd z(n);
  for (int i = 0; i < n; i++) {
    z(i) = Phi_inv(p(i));
  }
  return z;
}

double bivariate_normal_pdf(
    const double& x,
    const double& y,
    const double& r = 0.0,
    const double& e1 = 0.0,
    const double& s1 = 1.0,
    const double& e2 = 0.0,
    const double& s2 = 1.0
) {
  // if (s1 <= 0 || s2 <= 0) {
  //   Rcpp::stop("Standard deviation must be positive");
  // }
  double z1 = (x - e1) / s1; // calculate standardized value of x
  double z2 = (y - e2) / s2; // calculate standardized value of y
  // calculate PDF value
  double pdf = 1.0 / (2 * pi * s1 * s2 * std::sqrt(1 - r*r)) * 
    std::exp(-(z1*z1 - 2*r*z1*z2 + z2*z2 + 1e-8) / (2 * (1 - r*r) + 1e-8));
  return pdf; // return PDF value
}

double phi2(const double& z1, const double& z2, const double& r = 0.0) {
  // calculate PDF value
  double pdf = 1.0 / (2 * pi * std::sqrt(1 - r*r + + 1e-8)) * 
    std::exp(-(z1*z1 - 2*r*z1*z2 + z2*z2 + 1e-8) / (2 * (1 - r*r) + 1e-8));
  return pdf; // return PDF value
}

double bivariate_normal_cdf(
    const double& x,
    const double& y,
    const double& r = 0.0,
    const double& e1 = 0.0,
    const double& s1 = 1.0,
    const double& e2 = 0.0,
    const double& s2 = 1.0
) {
  // if (s1 <= 0 || s2 <= 0) {
  //   Rcpp::stop("Standard deviation must be positive");
  // }
  double z1 = (x - e1) / s1; // calculate standardized value of x
  double z2 = (y - e2) / s2; // calculate standardized value of y
  double cdf = bivnor(-z1, -z2, r);
  return cdf; // return CDF value
}

double Phi2(
    const double& z1,
    const double& z2,
    const double& r = 0.0
) {
  double cdf = bivnor(-z1, -z2, r);
  return cdf; // return CDF value
}

// ## MOMENTS
// #### Univariate
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

// #### Bivariate
double Cov(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
  Eigen::Index n = x.size();
  Eigen::VectorXd x_centered = x.array() - x.mean();
  Eigen::VectorXd y_centered = y.array() - y.mean();
  double covar = (x_centered.dot(y_centered)) / (n - 1);
  return covar;
}

double cor_pearson(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
  return Cov(x, y) / (SD(x) * SD(y));
}

double cor_spearman(const Eigen::VectorXd& x, const Eigen::VectorXd& y) {
  // Filter for pairwise complete observations
  // Eigen::Array<bool, Eigen::Dynamic, 1> pco = pairwise_complete(X, Y);
  // Eigen::VectorXd x = filter(X, pco);
  // Eigen::VectorXd y = filter(Y, pco);
  
  long long int n = x.size();
  Eigen::VectorXd d(n);
  double rho;
  // check for equal length
  // if (y.size() != n) {
  //   Rcpp::stop("Vectors are of different dimensionality");
  // }
  d = rank_(x).array() - rank_(y).array();
  rho = 1.0 - (6.0 * (d.dot(d))) / (n*(n*n - 1));
  rho = (rho < -1.0) ? -1.0 : rho;
  rho = (rho >  1.0) ?  1.0 : rho;
  return rho;
}

// Polychoric correlation class to be optimised
class Polychoric {
private:
  // double rho;                  // Polychoric correlation coefficient estimate
  // Eigen::VectorXd gamma;       // Thresholds for latent variable X
  // Eigen::VectorXd tau;         // Thresholds for latent variable Y
  Eigen::VectorXd g;              // Augmented thresholds: –∞ · gamma · ∞
  Eigen::VectorXd t;              // Augmented thresholds: –∞ · tau · ∞
  // Eigen::MatrixXd pmf;         // Matrix of hypothesized probabilities
  const Eigen::MatrixXd G;        // (adjusted) Observed contingency table
  int r;                          // # of rows
  int s;                          // # of columns
  
  // Probabilities of latent X and Y to fall into intervals formed by gamma and tau
  // Analogous to psych::polyBinBvn()
  Eigen::MatrixXd grid_discretised_normal_pmf(
      const Eigen::VectorXd& X,
      const Eigen::VectorXd& Y,
      const double& r = 0.0
  ) {
    int n = X.size() - 2;
    int m = Y.size() - 2;
    // if (s1 <= 0 || s2 <= 0) {
    //   Rcpp::stop("Standard deviation must be positive");
    // }
    
    // create matrix to store CDF values
    Eigen::MatrixXd P = Eigen::MatrixXd::Zero(n+2, m+2);
    for (int i = 1; i < n+2; i++) {
      double z1 = X(i);
      for (int j = 1; j < m+2; j++) {
        double z2 = Y(j);
        double cdf = bivnor(-z1, -z2, r);
        P(i, j) = cdf; // store CDF value in matrix
      }
    }
    P(n+1, m+1) = 1.0;
    // return P;
    Eigen::MatrixXd pmf_matrix(n+1, m+1);
    for (int i = 1; i < n+2; i++) {
      for (int j = 1; j < m+2; j++) {
        // Calculate P((x_i-1 < X < x_i) & (y_j-1 < Y < y_j))
        pmf_matrix(i-1, j-1) = (P(i, j) - P(i-1, j)) - (P(i, j-1) - P(i-1, j-1));
      }
    }
    return pmf_matrix; // return matrix of probabilities
  }
  
  // Log likelihood of the observed contingency table given hypothesised probabilities
  double contingency_table_loglik(
      const Eigen::MatrixXd& G,
      const Eigen::MatrixXd& P
  ) {
    // Compute log likelihood of contingency table G given matrix of probabilities P
    double loglik = 0.0;
    for (int i = 0; i < r; i++) {
      for (int j = 0; j < s; j++) {
        if (P(i, j) > 0.0) {
          // adjust value for the logarithm to be finite
          loglik += G(i, j) * std::log(P(i, j));
        }
      }
    }
    return -loglik;
  }
  
  // Partial derivative of log likelihood w/ respect to rho
  // Olsson 1979
  double dl_drho(
      const double& rho,
      const Eigen::VectorXd& gamma,
      const Eigen::VectorXd& tau,
      const Eigen::MatrixXd& G,
      const Eigen::MatrixXd& P
  ) {
    // Compute the partial derivative of the log likelihood with respect to rho
    double dl = 0.0;
    for (int i = 1; i < r+1; i++) {
      for (int j = 1; j < s+1; j++) {
        if (P(i-1,j-1) > 0.0) {
          double dphi1 = phi2(g(i), t(j), rho) - phi2(g(i-1), t(j), rho);
          double dphi2 = phi2(g(i), t(j-1), rho) - phi2(g(i-1), t(j-1), rho);
          double m = G(i-1,j-1)/(P(i-1, j-1)) * (dphi1 - dphi2);
          dl += m;
        }
      }
    }
    return -dl;
  }
  
public:
  // Define constructor
  Polychoric(
    const Eigen::MatrixXd& G_,    // Observed contingency table
    const Eigen::VectorXd& gamma, // Thresholds of latent variable X
    const Eigen::VectorXd& tau    // Thresholds of latent variable Y
  ) : G(G_) {
    r = G.rows();
    s = G.cols();
    
    // append large values to gamma and tau
    g = Eigen::VectorXd(r+1);
    g(0) = -1e6;
    for (int i = 1; i < r; i++) {
      g(i) = gamma(i-1);
    }
    g(r) = 1e6;
    
    t = Eigen::VectorXd(s+1);
    t(0) = -1e6;
    for (int j = 1; j < s; j++) {
      t(j) = tau(j-1);
    }
    t(s) = 1e6;
  }
  
  // Define the objective function for the optimizer
  double operator ()(const Eigen::VectorXd& x, Eigen::VectorXd& grad) {
    double rho = x(0);
    double ll = 0.0;
    Eigen::MatrixXd pmf = grid_discretised_normal_pmf(g, t, rho);
    ll = contingency_table_loglik(G, pmf);
    grad[0] = dl_drho(rho, g, t, G, pmf);
    return ll;
  }
};

// Function to estimate thresholds from an Eigen VectorXd object
Eigen::VectorXd estimate_thresholds(
    const Eigen::VectorXd& X,
    Eigen::Array<bool, Eigen::Dynamic, 1> keep = {},
    const double& correct = 1e-08
) {
  // Filter out NA values
  if (keep.size() == 0) {
    keep = !(X.array().isNaN());
  }
  Eigen::VectorXd x = filter(X, keep);
  
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
    counts(i) += counts(i-1);
  }
  Eigen::VectorXd p(m-1);
  p = counts.head(m-1).array() / N;
  Eigen::VectorXd q(m-1);
  q = Phi_inv_(p);
  return q;
}

// Function to estimate thresholds from a contingency table
std::vector<Eigen::VectorXd> estimate_thresholds(
    const Eigen::MatrixXd& G,
    const int& r, // # of rows
    const int& s, // # of columns
    const int& N  // # of observations
) {
  std::vector<Eigen::VectorXd> thresholds;
  
  // cumulative probabilities of latent X variable
  Eigen::VectorXd row_cumsum = G.rowwise().sum();
  for (int i = 1; i < r; i++) {
    row_cumsum(i) = row_cumsum(i) + row_cumsum(i-1);
  }
  // cumulative probabilities of latent Y variable
  Eigen::VectorXd col_cumsum = G.colwise().sum();
  for (int j = 1; j < s; j++) {
    col_cumsum(j) = col_cumsum(j) + col_cumsum(j-1);
  }
  
  // Compute the thresholds
  Eigen::VectorXd gamma = Phi_inv_(row_cumsum.head(r-1) / N);
  Eigen::VectorXd tau   = Phi_inv_(col_cumsum.head(s-1) / N);
  
  thresholds.push_back(gamma);
  thresholds.push_back(tau);
  
  return thresholds;
}

// Setting up an L-BFGS-B optimizer
double optim_polychoric(
    const Eigen::MatrixXd& G,
    const Eigen::VectorXd& gamma,
    const Eigen::VectorXd& tau
) {
  // Define the initial point for the optimization
  int n = 1;
  Eigen::VectorXd x0 = Eigen::VectorXd::Constant(n, 0.0);
  
  // Set the optimization parameters
  LBFGSpp::LBFGSBParam<double> params;
  params.max_iterations = 20;
  // params.max_linesearch = 99;
  params.epsilon = 1e-9;
  params.delta   = 1e-9;
  params.max_step = 0.9;
  params.min_step = 1e-10;
  Eigen::VectorXd lb = Eigen::VectorXd::Constant(n, -.999);
  Eigen::VectorXd ub = Eigen::VectorXd::Constant(n,  .999);
  
  // Create solver and function object
  LBFGSpp::LBFGSBSolver<double, LBFGSpp::LineSearchMoreThuente> solver(params);
  Polychoric poly(G, gamma, tau);
  
  // Run the optimization
  double f = 0.0;
  Eigen::VectorXd x = x0;
  int niter = solver.minimize(poly, x, f, lb, ub);
  
  // Return the optimal value of rho
  return x(0);
}

// Polychoric correlation given a contingency table
// [[Rcpp::export(.poly_tab)]]
double poly_tab(
    Eigen::MatrixXd G,
    const double& correct = 1e-08
) {
  int r = G.rows();
  int s = G.cols();
  int N = G.sum();
  // Check for perfect correlation first
  double mds = G.diagonal().sum();        // Main diagonal sum
  double sds = G.rowwise().reverse().diagonal().sum(); // Secondary diagonal sum
  if (mds == N) {return  1.0;}
  if (sds == N) {return -1.0;}
  
  // Correct for continuity
  for (int i = 0; i < r; i++) {
    for (int j = 0; j < s; j++) {
      double x = G(i,j);
      double v = (x == 0.0) ? (x + correct) : x;
      G(i,j) = v;
    }
  }
  
  std::vector<Eigen::VectorXd> thresholds = estimate_thresholds(G, r, s, N);
  Eigen::VectorXd gamma = thresholds[0];
  Eigen::VectorXd tau   = thresholds[1];
  
  double rho = optim_polychoric(G, gamma, tau);
  return(rho);
}

// Polychoric correlation given a contingency table to be called inside poly_df()
double poly_tab_inside(
    Eigen::MatrixXd G,
    const Eigen::VectorXd& gamma,
    const Eigen::VectorXd& tau,
    const double& correct = 1e-08
) {
  int r = G.rows();
  int s = G.cols();
  int N = G.sum();
  // Check for perfect correlation first
  double mds = G.diagonal().sum();                      // Main diagonal sum
  double sds = G.rowwise().reverse().diagonal().sum(); // Secondary diagonal sum
  if (mds == N) {return  1.0;}
  if (sds == N) {return -1.0;}
  
  // Correct for continuity
  for (int i = 0; i < r; i++) {
    for (int j = 0; j < s; j++) {
      double x = G(i,j);
      double v = (x == 0.0) ? (x + correct) : x;
      G(i,j) = v;
    }
  }
  
  double rho = optim_polychoric(G, gamma, tau);
  return(rho);
}

// Faster Eigen contingency table from a pair of Eigen::VectorXd
Eigen::MatrixXd contingency_table(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& y
) {
  int n = x.size();
  int k1 = x.maxCoeff();
  int k2 = y.maxCoeff();
  Eigen::MatrixXd table = Eigen::MatrixXd::Zero(k1, k2);
  for (int i = 0; i < n; i++) {
    table(x(i)-1, y(i)-1) += 1;
  }
  return table;
}

// Polychoric correlation given two ordinal (integer) vectors
// [[Rcpp::export(.poly_xy)]]
double poly_xy(
    const Eigen::VectorXd& X,
    const Eigen::VectorXd& Y,
    const double& correct = 1e-08
) {
  // Filter for pairwise complete observations
  // auto start = std::chrono::system_clock::now(); // TIME IT
  Eigen::Array<bool, Eigen::Dynamic, 1> pco = pairwise_complete(X, Y);
  Eigen::VectorXd x = filter(X, pco);
  Eigen::VectorXd y = filter(Y, pco);
  // auto end = std::chrono::system_clock::now(); // TIME IT
  // std::chrono::duration<double> elapsed_seconds = end-start;
  // Rcpp::Rcout << "Filtering:              " << elapsed_seconds.count() << "s" << std::endl;
  
  // Check if variables are discrete with reasonable number of levels
  // start = std::chrono::system_clock::now(); // TIME IT
  if ((n_unique(x) > 10) || (n_unique(y) > 10)) {
    Rcpp::warning(
      "Too many levels or continuous input: returning Spearman's rho\n"
    );
    return(cor_spearman(x, y));
  }
  // end = std::chrono::system_clock::now(); // TIME IT
  // elapsed_seconds = end-start;
  // Rcpp::Rcout << "Checking descreteness:  " << elapsed_seconds.count() << "s" << std::endl;
  
  // Check if data types are polychoric-friendly
  // start = std::chrono::system_clock::now(); // TIME IT
  x = correct_data(x);
  y = correct_data(y);
  // end = std::chrono::system_clock::now(); // TIME IT
  // elapsed_seconds = end-start;
  // Rcpp::Rcout << "Correcting values:      " << elapsed_seconds.count() << "s" << std::endl;
  
  // Construct a contingency table
  // start = std::chrono::system_clock::now(); // TIME IT
  Eigen::MatrixXd G = contingency_table(x, y);
  // end = std::chrono::system_clock::now(); // TIME IT
  // elapsed_seconds = end-start;
  // Rcpp::Rcout << "Constructing table:     " << elapsed_seconds.count() << "s" << std::endl;
  // Use the function for a table
  // start = std::chrono::system_clock::now(); // TIME IT
  // double rho = poly_tab(G, correct);
  // end = std::chrono::system_clock::now(); // TIME IT
  // elapsed_seconds = end-start;
  // Rcpp::Rcout << "Estimating coefficient: " << elapsed_seconds.count() << "s" << std::endl;
  return poly_tab(G, correct);
}

// Polychoric correlation given two ordinal vectors to be called inside poly_df()
double poly_xy_inside(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    const Eigen::VectorXd& gamma,
    const Eigen::VectorXd& tau,
    const double& correct = 1e-08
) {
  // Construct a contingency table
  Eigen::MatrixXd G = contingency_table(x, y);
  // Use the function for a table
  return poly_tab_inside(G, gamma, tau, correct);
}

// Polychoric correlation given a numeric matrix
Eigen::MatrixXd poly_mat(
    const Eigen::MatrixXd& X,  // A numeric (integer) matrix
    Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> shadow,
    Rcpp::List thresholds,     // Pre-estimated thresholds
    double correct = 1e-08     // Continuity correction value
) {
  const int m = X.cols(); // number of items
  const int N = X.rows(); // number of respondents
  
  Eigen::MatrixXd Corr = Eigen::MatrixXd::Ones(m, m);
  // Create item vectors and threshold templates
  Eigen::VectorXd x(N), y(N), gamma, tau;
  // Boolean vectors of whether an observation is not missing
  Eigen::Array<bool, Eigen::Dynamic, 1> x_complete, y_complete;
  
  // Check for discreteness
  bool cont_x, cont_y;
  double rho;
  
  // Populate the matrices and vectors with pairwise correlations
  for (int i = 1; i < m; i++) {
    x = X.col(i);
    gamma = thresholds[i];
    // I will search for reasons \ Turning pain to gold©
    x_complete = shadow.col(i);
    // Check if x is discrete
    cont_x = (gamma.size() == 0);
    for (int j = 0; j < i; j++) {
      y = X.col(j);
      tau = thresholds[j];
      // Feel the shadow call. And shimmer as the ice.©
      y_complete = shadow.col(j);
      // Filter for pairwise complete observations
      Eigen::VectorXd xc = filter(x, x_complete && y_complete);
      Eigen::VectorXd yc = filter(y, x_complete && y_complete);
      // If at least one variable suspected to be continuous
      cont_y = (tau.size() == 0);
      if (cont_x || cont_y) {
        Rcpp::warning(
          "Too many levels or continuous input: returning Spearman's rho\n"
        );
        rho = cor_spearman(xc, yc);
      } else {
        rho = poly_xy_inside(xc, yc, gamma, tau, correct);
      }
      Corr(i,j) = rho;
      Corr(j,i) = rho;
    }
  }
  // Return
  return Corr;
}

// Polychoric correlation matrix for a data frame of ordinal items
// [[Rcpp::export(.poly_df)]]
Eigen::MatrixXd poly_df(Rcpp::List X, double correct = 1e-08) {
  const int m = X.size();    // number of items
  
  // Convert to matrix
  Eigen::MatrixXd M = list_to_matrix(X);
  // Create matrix of non-missing values
  Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> shadow = shadow_matrix(M);
  // Create empty lists to store thresholds
  Rcpp::List thresholds(m); // Estimate threshold values early
  // Create threshold templates
  Eigen::VectorXd gamma, tau;
  
  // Populate the items and thresholds lists
  for (int i = 0; i < m; i++) {
    M.col(i) = correct_data(M.col(i));
    if (n_unique(M.col(i), shadow.col(i)) > 10) {
      gamma = {};
    } else {
      gamma = estimate_thresholds(M.col(i), shadow.col(i), correct);
    }
    thresholds[i] = gamma;
  }
  
  // Populate the matrices and vectors with pairwise correlations
  Eigen::MatrixXd Corr = poly_mat(M, shadow, thresholds, correct);
  // Return
  return Corr;
}

// p-value of the polychoric correlation coefficient
double poly_pval(const double& rho, const int& n) {
  if (rho == 1.0 || rho == -1.0) return 2e-16;
  
  double z = 0.5 * log((1 + rho) / (1 - rho));
  double sez = 1 / std::sqrt(n - 3);
  double p;
  
  if (z > 0) {
    p = (1.0 - normal_cdf(z, 0, sez)) * 2.0;
  } else {
    p = normal_cdf(z, 0, sez) * 2.0;
  }
  
  return p;
}

// [[Rcpp::export(.poly_tab_full)]]
Rcpp::List poly_tab_full(
    Eigen::MatrixXd G,
    const double& correct
) {
  // Return values
  double rho, pval;
  Rcpp::List res(4);
  res.attr("class") = "polychoric";
  
  int r = G.rows();
  int s = G.cols();
  int N = G.sum();
  
  // Correct for continuity
  for (int i = 0; i < r; i++) {
    for (int j = 0; j < s; j++) {
      double x = G(i,j);
      double v = (x == 0.0) ? (x + correct) : x;
      G(i,j) = v;
    }
  }
  
  std::vector<Eigen::VectorXd> thresholds = estimate_thresholds(G, r, s, N);
  Eigen::VectorXd gamma = thresholds[0];
  Eigen::VectorXd tau   = thresholds[1];
  
  // Check for perfect correlation first
  double mds = G.diagonal().sum();                     // Main diagonal sum
  double sds = G.rowwise().reverse().diagonal().sum(); // Secondary diagonal sum
  if (mds == N) {
    res = Rcpp::List::create(
      Rcpp::Named("rho")   = 1.0,
      Rcpp::Named("pval")  = 2e-16,
      Rcpp::Named("gamma") = gamma,
      Rcpp::Named("tau")   = tau
    );
    return res;
  }
  if (sds == N) {
    res = Rcpp::List::create(
      Rcpp::Named("rho")   = -1.0,
      Rcpp::Named("pval")  = 2e-16,
      Rcpp::Named("gamma") = gamma,
      Rcpp::Named("tau")   = tau
    );
    return res;
  }
  
  rho  = optim_polychoric(G, gamma, tau);
  pval = poly_pval(rho, N);
  
  res = Rcpp::List::create(
    Rcpp::Named("rho")   = rho,
    Rcpp::Named("pval")  = pval,
    Rcpp::Named("gamma") = gamma,
    Rcpp::Named("tau")   = tau
  );
  return res;
}

// [[Rcpp::export(.poly_xy_full)]]
Rcpp::List poly_xy_full(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    const double& correct = 1e-08
) {
  // Construct a contingency table
  Eigen::MatrixXd G = contingency_table(x, y);
  // Use the function for a table
  return poly_tab_full(G, correct);
}

// Polychoric correlation + supplementary information
// [[Rcpp::export(.poly_df_full)]]
Rcpp::List poly_df_full(Rcpp::List X, double correct = 1e-08) {
  const int m = X.size();    // number of items
  int nc = (m * (m-1)) / 2;  // number of correlations
  
  // Convert to matrix
  Eigen::MatrixXd M = list_to_matrix(X);
  int N = M.rows();          // number of respondents
  // Create matrix of non-missing values
  Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> shadow = shadow_matrix(M);
  // Create empty lists to store thresholds
  Rcpp::List thresholds(m); // Estimate threshold values early
  // Create threshold templates
  Eigen::VectorXd gamma, tau;
  
  // Populate the items and thresholds lists
  for (int i = 0; i < m; i++) {
    M.col(i) = correct_data(M.col(i));
    if (n_unique(M.col(i), shadow.col(i)) > 10) {
      gamma = {};
    } else {
      gamma = estimate_thresholds(M.col(i), shadow.col(i), correct);
    }
    thresholds[i] = gamma;
  }
  
  // Calculate matrix of coefficients
  Eigen::MatrixXd Corr = poly_mat(M, shadow, thresholds, correct);
  
  // Create p-value matrix template
  Eigen::MatrixXd Pval = Eigen::MatrixXd::Constant(m, m, 2e-16*nc);
  
  // Create matrix of adjusted p-values
  for (int i = 1; i < m; i++) {
    for (int j = 0; j < i; j++) {
      // Calculate p-value, not allowing zeros
      double pval = poly_pval(Corr(i,j), N);
      pval = std::max(2e-16, pval);
      // Correct p-value (Bonferroni)
      pval = std::min(1.0, pval*nc);
      Pval(i,j) = pval;
      Pval(j,i) = pval;
    }
  }
  // Return
  Rcpp::List res = Rcpp::List::create(
    Rcpp::Named("rho")  = Corr,
    Rcpp::Named("pval") = Pval,
    Rcpp::Named("tau")  = thresholds
  );
  res.attr("class") = "polychoric";
  return res;
}

// Polyserial correlation class to be optimised
class Polyserial {
private:
  Eigen::VectorXd x;    // Continuous variable X
  Eigen::VectorXd d;    // Ordinal variable D
  int n;                // # of observations
  int s;                // # of category levels in D
  Eigen::VectorXd z;    // Standardized values (z-scores) of x
  double correct;       // Correct for continuity
  
  // Pr(D = d_k | x_k)
  Eigen::VectorXd polyserial_pd(const double& rho, const Eigen::VectorXd& rz) {
    // Declare vector to return
    Eigen::VectorXd P_d_x(n);
    // Constant components
    double denom = std::sqrt(1.0 - rho*rho + 1e-8);
    
    // Iterate over observations
    for (int k = 0; k < n; k++) {
      int j = d(k)-1;
      double tau_star_j   = (d(k)<s) ? (tau(j)   - rz(k)) / denom :  1e6;
      double tau_star_jm1 = (j>1)    ? (tau(j-1) - rz(k)) / denom : -1e6;
      P_d_x(k) = Phi(tau_star_j) - Phi(tau_star_jm1);
    }
    return P_d_x;
  }
  
  // Log likelihood of observations given hypothesised probabilities
  double polyserial_loglik(const Eigen::VectorXd& Pd) {
    
    Eigen::VectorXd logP = (Pd.array() + 1e-8).log();
    Eigen::VectorXd logfx = (normal_pdf_(x, mu, sigma).array() + 1e-8).log();
    double loglik = (logfx + logP).sum();
    
    return -loglik;
  }
  
  // Partial derivative of log likelihood w/ respect to rho
  // Drasgow 1986
  double dl_drho(
      const double& rho,
      const Eigen::VectorXd& Pd,
      const Eigen::VectorXd& rz
  ) {
    // Compute the partial derivative of the log likelihood with respect to rho
    double dl = 0.0;
    // Constants depending on rho only
    double denom = std::sqrt(1.0 - rho*rho + 1e-6);
    double dtau_denom = 1 / std::pow(1.0 - rho*rho + 1e-6, 3/2);
    
    for (int k = 0; k < n; k++) {
      int j = d(k)-1;
      // Rename tau_star for t for compactness
      double phi_t_j   = (d(k)<s) ? phi((tau(j)   - rz(k)) / denom) : 0.0;
      double phi_t_jm1 = (j>1)    ? phi((tau(j-1) - rz(k)) / denom) : 0.0;
      
      if (Pd(k) > 0.0) {
        double dlk = (1/Pd(k) * dtau_denom) *
          ((phi_t_j*(tau(j)*rho - z(k))) - (phi_t_jm1*(tau(j-1)*rho - z(k))));
        dl += dlk;
      }
    }
    return -dl;
  }
  
public:
  // Variable parameters made available for the wrapper output
  double rho;           // Optimal value of rho
  Eigen::VectorXd tau;  // Thresholds of Y underlying D
  double mu;            // Mean of X
  double sigma;         // Standard deviation of X
  
  // Define constructor
  Polyserial(
    const Eigen::VectorXd& X, // Continuous covariate
    const Eigen::VectorXd& D, // Ordinal covariate
    const double& correct     // Continuity correction value
  ) : x(X), d(D), correct(correct) {
    n = d.size();
    s = d.maxCoeff();
    z = (x.array() - x.mean()) / SD(x);
    // Assert missing values were already filtered out
    Eigen::Array<bool, Eigen::Dynamic, 1> mask = Eigen::Array<bool, Eigen::Dynamic, 1>::Ones(n).cast<bool>();
    tau = estimate_thresholds(d, mask, correct);
    mu = x.mean();
    sigma = SD(x);
  }
  
  // Define the objective function for the optimizer
  double operator ()(const Eigen::VectorXd& x, Eigen::VectorXd& grad) {
    double rho = x(0);
    double ll = 0.0;
    Eigen::VectorXd rz = z.array() * rho;    // \rho × z
    Eigen::VectorXd Pd = polyserial_pd(rho, rz);
    ll = polyserial_loglik(Pd);
    grad[0] = dl_drho(rho, Pd, rz);
    return ll;
  }
};

// Setting up an L-BFGS-B optimizer
Polyserial optim_polyserial(
    const Eigen::VectorXd& X,
    const Eigen::VectorXd& D,
    double correct = 1e-08
) {
  // Filter and correct
  Eigen::Array<bool, Eigen::Dynamic, 1> pco = pairwise_complete(X, D);
  Eigen::VectorXd x = filter(X, pco);
  Eigen::VectorXd d = correct_data(filter(D, pco));
  
  // Define the initial point for the optimization: Pearson's coefficient
  int n = 1;
  Eigen::VectorXd x0 = Eigen::VectorXd::Constant(n, cor_pearson(x, d));
  
  // Set the optimization parameters
  LBFGSpp::LBFGSBParam<double> params;
  params.max_iterations = 20;
  params.epsilon = 1e-8;
  params.delta   = 1e-8;
  params.max_step = 1e20;
  params.min_step = 1e-20;
  Eigen::VectorXd lb = Eigen::VectorXd::Constant(n, -1.);
  Eigen::VectorXd ub = Eigen::VectorXd::Constant(n,  1.);
  
  // Create solver and function object
  LBFGSpp::LBFGSBSolver<double, LBFGSpp::LineSearchMoreThuente> solver(params);
  Polyserial poly(x, d, correct);
  
  // Run the optimization
  double f = 0.0;
  Eigen::VectorXd rho = x0;
  int niter = 0;
  
  // Optimizer may fail when Pearson's r is too good an estimate
  try {
    niter = solver.minimize(poly, rho, f, lb, ub);
    // Append the optimal value of rho
    poly.rho = rho(0);
  } catch (const std::domain_error& e) {
    Rcpp::Rcout << "Optimizer problem: “" << e.what() << "”" << std::endl; 
    Rcpp::warning("The optimizer failed to improve the initial guess: \
                    returning Pearson's coefficient\n");
    poly.rho = x0(0);
  }
  
  return poly;
}

// For `coef.only = TRUE`
// [[Rcpp::export(.cor_polyserial)]]
double cor_polyserial(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& d,
    const double& correct = 1e-08
) {
  Polyserial correlation = optim_polyserial(x, d, correct);
  return correlation.rho;
}

// For `coef.only = FALSE`
// [[Rcpp::export(.cor_polyserial_full)]]
Rcpp::List cor_polyserial_full(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& d,
    const double& correct = 1e-08
) {
  Polyserial correlation = optim_polyserial(x, d, correct);
  
  double pval = std::max(2e-16, poly_pval(correlation.rho, d.size()));
  // Parameters
  Rcpp::List params = Rcpp::List::create(
    Rcpp::Named("tau")   = correlation.tau,
    Rcpp::Named("mu")    = correlation.mu,
    Rcpp::Named("sigma") = correlation.sigma
  );
  
  Rcpp::List res = Rcpp::List::create(
    Rcpp::Named("rho")  = correlation.rho,
    Rcpp::Named("pval") = pval,
    Rcpp::Named("parameters") = params
  );
  res.attr("class") = "polyserial";
  
  return res;
}