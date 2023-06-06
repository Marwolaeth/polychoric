// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#include <Rcpp.h>
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

double phi2(
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
    std::exp(-(z1*z1 - 2*r*z1*z2 + z2*z2 + 1e-12) / (2 * (1 - r*r) + 1e-12));
  return pdf; // return PDF value
}

double Phi2(
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
  return cdf; // return PDF value
}

// Polychoric correlation and thresholds as class to be returned
class Polychoric {
private:
  // double rho;            // Polychoric correlation coefficient estimate
  Eigen::VectorXd gamma;    // Thresholds for latent variable X
  Eigen::VectorXd tau;      // Thresholds for latent variable Y
  Eigen::VectorXd g;        // Augmented thresholds: –∞ · gamma · ∞
  Eigen::VectorXd t;        // Augmented thresholds: –∞ · tau · ∞
  // Eigen::MatrixXd pmf;   // Matrix of hypothesized probabilities
  Eigen::MatrixXd G;        // (adjusted) Observed contingency table
  
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
  
  // Log likelihood given true contingency table and hypothesised probabilities
  double contingency_table_loglik(
      const Eigen::MatrixXd& G,
      const Eigen::MatrixXd& P
  ) {
    // Compute log likelihood of contingency table G given matrix of probabilities P
    int r = G.rows();
    int s = G.cols();
    // if (r != P.rows() || s != P.cols()) {
    //   Rcpp::stop("Input matrices must have same dimensions");
    // }
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
  double dl_drho(
      const double& rho,
      const Eigen::VectorXd& gamma,
      const Eigen::VectorXd& tau,
      const Eigen::MatrixXd& G,
      const Eigen::MatrixXd& P
  ) {
    int r = G.rows();
    int s = G.cols();
    
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
    const Eigen::MatrixXd& G_,
    const Eigen::VectorXd& gamma,
    const Eigen::VectorXd& tau
  ) : G(G_) {
    int r = G.rows();
    int s = G.cols();
    
    // append large values to gamma and tau
    g = Eigen::VectorXd(r+1);
    g(0) = -1e9;
    for (int i = 1; i < r; i++) {
      g(i) = gamma(i-1);
    }
    g(r) = 1e9;
    
    t = Eigen::VectorXd(s+1);
    t(0) = -1e9;
    for (int j = 1; j < s; j++) {
      t(j) = tau(j-1);
    }
    t(s) = 1e9;
  }
  
  // Define the objective function
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
    int N,
    const double& correct = 0.1
) {
  // Filter out NA values
  Eigen::VectorXd x = filter(X, !is_na(X));
  N = x.size();
  
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
  Eigen::VectorXd gamma = Phi_inv_vec(row_cumsum.head(r-1) / N);
  Eigen::VectorXd tau   = Phi_inv_vec(col_cumsum.head(s-1) / N);
  
  thresholds.push_back(gamma);
  thresholds.push_back(tau);
  
  return thresholds;
}

// Setting up an L-BFGS-B optimizer
double poly_optim(
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
  params.delta = 1e-9;
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
    const double& correct = 0.1
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
  
  double rho = poly_optim(G, gamma, tau);
  return(rho);
}

// Polychoric correlation given a contingency table
double poly_tab_inside(
    Eigen::MatrixXd G,
    const Eigen::VectorXd& gamma,
    const Eigen::VectorXd& tau,
    const double& correct = 0.1
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
  
  double rho = poly_optim(G, gamma, tau);
  return(rho);
}

// Faster Eigen contingency table from Eigen::VectorXd
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

// Polychoric correlation of two ordinal (integer) vectors
// [[Rcpp::export(.poly_xy)]]
double poly_xy(
    const Eigen::VectorXd& X,
    const Eigen::VectorXd& Y,
    double correct = 0.1
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
      "Too many levels or continuous input: returning Spearman's rho"
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

// Polychoric correlation to be used poly_df()
double poly_xy_inside(
    const Eigen::VectorXd& x,
    const Eigen::VectorXd& y,
    const Eigen::VectorXd& gamma,
    const Eigen::VectorXd& tau,
    double correct = 0.1
) {
  // Construct a contingency table
  Eigen::MatrixXd G = contingency_table(x, y);
  // Use the function for a table
  return poly_tab_inside(G, gamma, tau, correct);
}

// Polychoric correlation matrix for a data frame of ordinal items
// [[Rcpp::export(.poly_df)]]
Eigen::MatrixXd poly_df(Rcpp::List X, double correct = 0.1) {
  int m = X.size();                           // number of items
  int N = ((Rcpp::IntegerVector)X[0]).size(); // number of respondents
  
  // Create correlation matrix template
  Eigen::MatrixXd Corr = Eigen::MatrixXd::Ones(m, m);
  Eigen::VectorXd x(N), y(N), gamma, tau;
  // Create empty lists to store items and their thresholds
  Rcpp::List items(m);
  Rcpp::List thresholds(m);
  bool cont_x, cont_y;
  double rho;
  
  // Populate the items and thresholds lists
  for (int i = 0; i < m; i++) {
    x = correct_data(X[i]);
    items[i] = x;
    if (n_unique(x) > 10) {
      gamma = {};
    } else {
      gamma = estimate_thresholds(x, N, correct);
    }
    thresholds[i] = gamma;
  }
  
  // Populate the matrices and vectors with pairwise correlations
  for (int i = 1; i < m; i++) {
    x = items[i];
    gamma = thresholds[i];
    cont_x = (gamma.size() == 0);
    Eigen::Array<bool, Eigen::Dynamic, 1> x_complete = !is_na(x);
    for (int j = 0; j < i; j++) {
      y = items[j];
      tau = thresholds[j];
      // Filter for pairwise complete observations
      Eigen::Array<bool, Eigen::Dynamic, 1> y_complete = !is_na(y);
      Eigen::VectorXd xc = filter(x, x_complete && y_complete);
      Eigen::VectorXd yc = filter(y, x_complete && y_complete);
      // If at least one variable suspected to be continuous
      cont_y = (tau.size() == 0);
      if (cont_x || cont_y) {
        Rcpp::warning(
          "Too many levels or continuous input: returning Spearman's rho"
        );
        rho = cor_spearman(x, y);
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

double poly_pval(const double& rho, const int& n) {
  if (rho == 1.0 || rho == -1.0) return 2e-16;
  
  double z = 0.5 * log((1 + rho) / (1 - rho));
  double sez = 1 / std::sqrt(n - 3);
  double p;
  
  if (z > 0) {
    p = (1.0 - Phi(z, 0, sez)) * 2.0;
  } else {
    p = Phi(z, 0, sez) * 2.0;
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
  double mds = G.diagonal().sum();        // Main diagonal sum
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
  
  rho = poly_optim(G, gamma, tau);
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
    double correct = 0.1
) {
  // Construct a contingency table
  Eigen::MatrixXd G = contingency_table(x, y);
  // Use the function for a table
  return poly_tab_full(G, correct);
}

// Polychoric correlation + supplementary information
// [[Rcpp::export(.poly_df_full)]]
Rcpp::List poly_df_full(Rcpp::List X, double correct = 0.1) {
  int m = X.size();                           // number of items
  int N = ((Rcpp::IntegerVector)X[0]).size(); // number of respondents
  int nc = (m * (m-1)) / 2;                   // number of correlations
  
  // Create correlation matrix template
  Eigen::MatrixXd Corr = Eigen::MatrixXd::Ones(m, m);
  Eigen::VectorXd x(N), y(N), gamma, tau;
  // Create p-value matrix template
  Eigen::MatrixXd Pval = Eigen::MatrixXd::Constant(m, m, 2e-16*nc);
  // Create empty lists to store items and their thresholds
  Rcpp::List items(m);
  Rcpp::List thresholds(m);
  bool cont_x, cont_y;
  double rho;
  
  // Populate the items and thresholds lists
  for (int i = 0; i < m; i++) {
    x = correct_data(X[i]);
    items[i] = x;
    if (n_unique(x) > 10) {
      gamma = {};
    } else {
      gamma = estimate_thresholds(x, N, correct);
    }
    thresholds[i] = gamma;
  }
  
  // Populate the matrices and vectors with pairwise correlations
  for (int i = 1; i < m; i++) {
    x = items[i];
    gamma = thresholds[i];
    cont_x = (gamma.size() == 0);
    Eigen::Array<bool, Eigen::Dynamic, 1> x_complete = !is_na(x);
    for (int j = 0; j < i; j++) {
      y = items[j];
      tau = thresholds[j];
      // Filter for pairwise complete observations
      Eigen::Array<bool, Eigen::Dynamic, 1> y_complete = !is_na(y);
      Eigen::VectorXd xc = filter(x, x_complete && y_complete);
      Eigen::VectorXd yc = filter(y, x_complete && y_complete);
      // If at least one variable suspected to be continuous
      cont_y = (tau.size() == 0);
      if (cont_x || cont_y) {
        Rcpp::warning(
          "Too many levels or continuous input: returning Spearman's rho"
        );
        rho = cor_spearman(x, y);
      } else {
        rho = poly_xy_inside(xc, yc, gamma, tau,correct);
      }
      Corr(i,j) = rho;
      Corr(j,i) = rho;
      // Calculate p-value
      double pval = poly_pval(rho, xc.size());
      // Correct p-value (Bonferroni)
      pval = std::min(1.0, pval*nc);
      Pval(i,j) = pval;
      Pval(j,i) = pval;
    }
  }
  // Return
  Rcpp::List res = Rcpp::List::create(
    Rcpp::Named("rho") = Corr,
    Rcpp::Named("pval") = Pval,
    Rcpp::Named("tau") = thresholds
  );
  res.attr("class") = "polychoric";
  return res;
}
