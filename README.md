# polychoric
Fast Polychoric Correlation for Vectors and Data Frames

## About
Polychoric is a package that provides a wrapper for C++ routines used to calculate polychoric correlation coefficients, which are often used in social science or marketing research. The single `polycorr()` function can take in ordinal factor vectors, a contingency table, or a data frame and returns corresponding polychoric correlation estimates in the form of a single numeric value or a correlation matrix.

<details>
  <summary>Spoiler</summary>
  
  *Pathetic amateurish craftsmanship*
 </details>
  
## The Purpose
Polychoric correlation coefficients are a type of correlation coefficient used to measure the relationship between two ordinal variables. They are computed by estimating the correlation between two underlying continuous variables that are assumed to give rise to the observed ordinal data. The `polycorr()` function estimates latent Pearson correlation coefficients under the assumption that the latent traits of interest are standard normal random variables.

The computation of polychoric correlation coefficients involves estimating the thresholds (here called `gamma` and `tau`, like in (*Drasgow 1986*), or just `tau` as in `psych` package (*Revelle 2023*)) that separate the ordinal categories for each variable. These thresholds are used to transform the ordinal data into a set of continuous variables, which can then be used to estimate the correlation coefficient using standard methods. The `polycorr()` function currently estimates the thresholds using a two-step maximum likelihood estimation, where first the thresholds are deduced from univariate distributions of ordinal variables and then the `L-BFGS-B` optimization algorithm (*Yixuan 2023*) is used to find the value of the correlation coefficient $\rho$ that maximizes the likelihood of the observed contingency table. The `toms462` (*Donnelly 1973*, *Owen 1956*) [algorithm](https://people.sc.fsu.edu/~jburkardt/cpp_src/toms462/toms462.html) is used to approximate the bivariate normal distribution (quadrant probabilities) of threshold values.

Polychoric correlation coefficients are useful for analyzing data that involve ordinal variables, such as Likert scales or survey responses. They provide a measure of the strength and direction of the relationship between two ordinal variables, which can be useful for understanding patterns in the data.

## Examples

```
## GSS 2012 Cultural Module
data("gss12_values", package = 'polychoric')
# A pair of discrete vectors
polycorr(gss12_values$valorig, gss12_values$valeql)
polycorr(gss12_values$valorig, gss12_values$valeql, coef.only = FALSE)

# A contingency table
(G <- table(gss12_values$valorig, gss12_values$valeql))
polycorr(G)

# A data.frame
polycorr(gss12_values)

# For safety, returns Spearman's rho if at least one of the vectors
# is presumably continuous (n_distinct(x) > 10 | n_distinct(y) > 10)
# (with a warning)
x <- rnorm(nrow(gss12_values))
polycorr(gss12_values$valspl, x)
```

## To do
This is an alpha version of the package. It is under development.

<details>
  <summary>The upcoming steps</summary>
  
  1. Optimise `Polychoric` class in the source code.
  2. Provide a polyserial correlation estimation function.
  3. Implement (optional) more robust distributional assumptions, e.g. a skew normal distribution (*Jin and Yang-Wallentin 2017*).
 </details>

## References
1. Drasgow, Fritz (1986). Polychoric and polyserial correlations. The Encyclopedia of Statistics. 7. 68-74.

2. Revelle, William (2023). _psych: Procedures for Psychological, Psychometric, and Personality Research_. Northwestern University, Evanston, Illinois. R package version 2.3.3, <https://CRAN.R-project.org/package=psych>.

3. Olsson, Ulf (1979). Maximum Likelihood Estimation of the Polychoric Correlation Coefficient, Psychometrika, 44:443-460.

4. Donnelly, Thomas (1973). Algorithm 462: Bivariate Normal Distribution, Communications of the ACM, October 1973, Volume 16, Number 10, page 638.

5. Owen, Donald (1956). Tables for Computing Bivariate Normal Probabilities, Annals of Mathematical Statistics, December 1956, Volume 27, Number 4, pages 1075-1090.

6. Qiu, Yixuan (2023). LBFGS++: A Header-only C++ Library for L-BFGS and L-BFGS-B Algorithms. Available at: https://lbfgspp.statr.me/

7. Jin S, Yang-Wallentin F (2017). Asymptotic Robustness Study of the Polychoric Correlation Estimation. Psychometrika. 2017 Mar;82(1):67-85.
