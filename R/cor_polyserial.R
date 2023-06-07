#' Estimate polyserial correlation coefficients between a continuous and an ordinal variables
#'
#' @param x A numeric vector, generally a sample from a normally distributed variable.
#' @param d A vector of discrete scores: ordinal or integer.
#' @param correct Correction value to use to correct for continuity in the case of zero entry cell of a contingency table. Can be used with any `x` input type.
#' @param coef.only If `TRUE`, returns only correlation coefficients (see Value).
#'
#' @return Polyserial correlation coefficients: numeric of length one (if `coef.only=TRUE`) or a list containing the estimated coefficient `rho`, the p-value `pval` and a list with estimated parameters of the variables used: `tau` for estimated thresholds of ordinal variable `d`, `mu` and `sigma` are the mean and the (estimated population) standard deviation of continuous variable `x`.
#' @export
#'
#' @details
#' Polyserial correlation coefficients are a type of correlation coefficient used to measure the relationship between a continuous variable and an ordinal variable. They are computed by estimating the correlation between the observed continuous variable and a latent continuous variable that is derived from the observed ordinal variable. The `cor_polyserial()` function estimates a latent Pearson correlation coefficient under the assumption that the latent traits of interest are standard normal random variables (Drasgow 1986).
#' 
#' The computation of polyserial correlation coefficients involves estimating the threshold values `tau` that separate the ordinal categories for the variable. The `cor_polyserial()` function currently estimates the thresholds using a two-step maximum likelihood estimation, where first the thresholds are deduced from the univariate distribution of the ordinal variable and then the `L-BFGS-B` optimization algorithm (Yixuan 2023) is used to find the value of the correlation coefficient `rho` that maximizes the likelihood of the observed data.
#'
#' @note
#' The `cor_polyserial()` function always uses pairwise complete observations. Therefore, the user need not worry about missing data. However, depending on the analysis design and the ratio of missing data, it may be essential to check for patterns of missingness and consider imputation.
#'
#' @references Drasgow, Fritz (1986). Polychoric and polyserial correlations. The Encyclopedia of Statistics. 7. 68-74.
#' @references Qiu, Yixuan (2023). LBFGS++: A Header-only C++ Library for L-BFGS and L-BFGS-B Algorithms. Available at: https://lbfgspp.statr.me/
#'
#' @examples
#' n <- 200
#' x <- rnorm(n, 100, 11)
#' y <- (x - 100) / 11 + rnorm(n)
#' C <- cut(x, breaks = 5, ordered_result = TRUE)
#' D <- cut(y, breaks = 5, ordered_result = TRUE)
#' cor_polyserial(x, D)
#' cor_polyserial(x, D, coef.only = FALSE)
#' cor_polyserial(x, C) # currently unable to handle perfect correlation
#' cor_polyserial(y, C)
#' cor_polyserial(y, D) # currently unable to handle perfect correlation
#' 
#' @keywords category
#' @keywords multivariate
#' @keywords nonparametric
#' @concept correlation
#' @concept polyserial
#' @concept polyserial correlation
#' @concept C++
#' @concept social science
#' @concept ordinal
#' @concept ordinal data
#' @concept ordinal factor
#' @concept likert
#' @concept likert scale
cor_polyserial <- function(
    x,
    d,
    correct = 1e-08,
    coef.only = TRUE
) {
  if (coef.only) {
    return(.cor_polyserial(x, d, correct))
  } else {
    return(.cor_polyserial_full(x, d, correct))
  }
}

