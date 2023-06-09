% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cor_polychoric.R
\name{cor_polychoric}
\alias{cor_polychoric}
\title{Estimate polychoric correlation coefficients}
\usage{
cor_polychoric(x, d = NULL, correct = 1e-08, coef.only = TRUE)
}
\arguments{
\item{x}{A vector of discrete scores: ordinal or integer, a contingency table or a data.frame (or a tibble) of discrete scores.}

\item{d}{A vector of discrete scores: ordinal or integer; only used when `x` is a vector.}

\item{correct}{Correction value to use to correct for continuity in the case of zero entry cell of a contingency table. Can be used with any `x` input type.}

\item{coef.only}{If `TRUE`, returns only correlation coefficients (see Value).}
}
\value{
Polychoric correlation coefficients: numeric of length one (for a pair of vectors or a contingency table) or an m×m correlation matrix (if `x` is a data.frame), where m is the number of items in the dataset. `coef.only=FALSE`, returns a list with (a matrix of) coefficients, a list of threshold estimates for every item used (length two for a pair of vectors or a table, length m for a data.frame) and (a matrix of Bonferroni-adjusted) p-values of correlation coefficients.
}
\description{
Estimate polychoric correlation coefficients
}
\details{
Polychoric correlation coefficients are a type of correlation coefficient used to measure the relationship between two ordinal variables. They are computed by estimating the correlation between two underlying continuous variables that are assumed to give rise to the observed ordinal data. The `cor_polychoric()` function estimates latent Pearson correlation coefficients under the assumption that the latent traits of interest are standard normal random variables.

The computation of polychoric correlation coefficients involves estimating the thresholds (here called `gamma` and `tau`, like in (Drasgow 1986), or just `tau` as in `psych` package (Revelle 2023)) that separate the ordinal categories for each variable. These thresholds are used to transform the ordinal data into a set of continuous variables, which can then be used to estimate the correlation coefficient using standard methods. The `cor_polychoric()` function currently estimates the coefficient using a two-step maximum likelihood estimation, where first the thresholds are deduced from univariate distributions of ordinal variables and then the `L-BFGS-B` optimization algorithm (Yixuan 2023) is used to find the value of the correlation coefficient `rho` that maximizes the likelihood of the observed contingency table. The `toms462` (Donnelly 1973, Owen 1956) algorithm is used to approximate the bivariate normal distribution (quadrant probabilities) of threshold values.
}
\note{
The `cor_polychoric()` function always uses pairwise complete observations. Therefore, the user need not worry about missing data. However, depending on the analysis design and the ratio of missing data, it may be essential to check for patterns of missingness and consider imputation.
}
\examples{
## GSS 2012 SCHWARTZ VALUES MODULE ----
data("gss12_values", package = 'polychoric')
# A pair of discrete vectors
cor_polychoric(gss12_values$valorig, gss12_values$valeql)
cor_polychoric(gss12_values$valorig, gss12_values$valeql, coef.only = FALSE)

# A contingency table
(G <- table(gss12_values$valorig, gss12_values$valeql))
cor_polychoric(G)

# A data.frame
cor_polychoric(gss12_values)

# For safety, returns Spearman's rho if at least one of the vectors
# is presumably continuous (n_distinct(x) > 10 | n_distinct(d) > 10)
# (with a warning)
y <- rnorm(nrow(gss12_values))
cor_polychoric(gss12_values$valspl, y)

}
\references{
Drasgow, Fritz (1986). Polychoric and polyserial correlations. The Encyclopedia of Statistics. 7. 68-74.

Revelle, William (2023). _psych: Procedures for Psychological, Psychometric, and Personality Research_. Northwestern University, Evanston, Illinois. R package version 2.3.3, <https://CRAN.R-project.org/package=psych>.

Olsson, Ulf (1979). Maximum Likelihood Estimation of the Polychoric Correlation Coefficient, Psychometrika, 44:443-460.

Donnelly, Thomas (1973). Algorithm 462: Bivariate Normal Distribution, Communications of the ACM, October 1973, Volume 16, Number 10, page 638.

Owen, Donald (1956). Tables for Computing Bivariate Normal Probabilities, Annals of Mathematical Statistics, December 1956, Volume 27, Number 4, pages 1075-1090.

Qiu, Yixuan (2023). LBFGS++: A Header-only C++ Library for L-BFGS and L-BFGS-B Algorithms. Available at: https://lbfgspp.statr.me/
}
\concept{C++}
\concept{correlation}
\concept{correlation matrix}
\concept{likert}
\concept{likert scale}
\concept{ordinal}
\concept{ordinal data}
\concept{ordinal factor}
\concept{polychoric}
\concept{polychoric correlation}
\concept{social science}
\concept{survey}
\keyword{category}
\keyword{multivariate}
\keyword{nonparametric}
