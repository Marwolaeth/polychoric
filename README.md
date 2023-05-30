# polychoric
Fast Polychoric Correlation for Vectors and Data Frames

## About
A wrapper for C++ routines for calculating polychoric correlation coefficients often used in social science or marketing research. Single function `polycorr()` can take in ordinal factor (possibly integer) vectors, a contingency table or a data frame. It returns corresponding polychoric correlation estimates in a form of single numeric value or correlation matrix.

<details>
  <summary>Spoiler</summary>
  
  *Pathetic amateurish craftsmanship*
 </details>
  
## The Purpose
Polychoric correlation coefficients are a type of correlation coefficient used to measure the relationship between two ordinal variables. They are computed by estimating the correlation between two underlying continuous variables that are assumed to give rise to the observed ordinal data. This function currently only estimates latent Pearson correlation coefficients under the assumption that the latent threats of interest are standard normal random variables.

The computation of polychoric correlation coefficients involves estimating the thresholds (here called `gamma` and `tau`, like in (Drasgow 1986), or just `tau` as in `psych` package (Revelle 2023)) that separate the ordinal categories for each variable. These thresholds are used to transform the ordinal data into a set of continuous variables, which can then be used to estimate the correlation coefficient using standard methods.

There are several methods for computing polychoric correlation coefficients, including maximum likelihood estimation and Bayesian estimation. These methods involve different assumptions about the distribution of the underlying continuous variables and the relationship between them. This package currently only uses a two-step maximum likelihood estimation, where first the thresholds are deduced from univariate distributions of ordinal variables and then L-BFGS-B optimisation algorithm (Yixuan 2023) is used to find the value of the correlation coefficient $\rho$ that maximises the likelihood of the observed contingency table.

Polychoric correlation coefficients are useful for analyzing data that involve ordinal variables, such as Likert scales or survey responses. They provide a measure of the strength and direction of the relationship between two ordinal variables, which can be useful for understanding patterns in the data and making predictions about future observations.

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

## References
1. Drasgow, Fritz (1986). Polychoric and polyserial correlations. The Encyclopedia of Statistics. 7. 68-74.

2. Revelle, William (2023). _psych: Procedures for Psychological, Psychometric, and Personality Research_. Northwestern University, Evanston, Illinois. R package version 2.3.3, <https://CRAN.R-project.org/package=psych>.

3. Olsson, U (1979). Maximum Likelihood Estimation of the Polychoric Correlation Coefficient, Psychometrika, 44:443-460.

4. Qiu, Yixuan (2023). LBFGS++: A Header-only C++ Library for L-BFGS and L-BFGS-B Algorithms. Available at: https://lbfgspp.statr.me/
