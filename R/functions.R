#' Estimate polychoric correlation coefficients
#'
#' @param x A vector of discrete scores: ordinal or integer, a contingency table or a data.frame (or a tibble) of discrete scores.
#' @param y A vector of discrete scores: ordinal or integer; only used when x is a vector.
#' @param correct Correction value to use to correct for continuity in the case of zero entry cell of a contingency table. Can be used with any x input type.
#' @param coef.only If TRUE, returns only correlation coefficients (see Value).
#'
#' @return Polychoric correlation coefficients: numeric of length one (for a pair of vectors or a contingency table) or an m×m correlation matrix (if x is a data.frame), where m is the number of items in the dataset. coef.only=FALSE, returns a list with (a matrix of) coefficients, a list of threshold estimates for every item used (length two for a pair of vectors or a table, length m for a data.frame) and (a matrix of) p-values of correlation coefficients.
#' @export
#'
#' @details
#' 
#'
#' @examples
#' ## GSS 2012 Cultural Module
#' data("gss12_values", package = 'polychoric')
#' # A pair of discrete vectors
#' polycorr(gss12_values$valorig, gss12_values$valeql)
#' polycorr(gss12_values$valorig, gss12_values$valeql, coef.only = FALSE)
#' 
#' # A contingency table
#' (G <- table(gss12_values$valorig, gss12_values$valeql))
#' polycorr(G)
#' 
#' # A data.frame
#' polycorr(gss12_values)
#' 
#' # For safety, returns Spearman's rho if at least one of the vectors
#' # is presumably continuous (n_distinct(x) > 10 | n_distinct(y) > 10)
#' # (with a warning)
#' x <- rnorm(nrow(gss12_values))
#' polycorr(gss12_values$valspl, x)
polycorr <- function(x, y = NULL, correct = 0.1, coef.only = TRUE) {
  x <- na.omit(x)
  if (!is.null(y)) {
    y <- na.omit(y)
    if (coef.only) {
      return(.poly_xy(x, y, correct = correct))
    } else {
      return(.poly_xy_full(x, y, correct = correct))
    }
  }
  if (is.matrix(x)) {
    if (coef.only) {
      return(.poly_tab(x, correct = correct))
    } else {
      return(.poly_tab_full(x, correct = correct))
    }
  }
  items <- names(x)
  if (coef.only) {
    rho <- .poly_df(x, correct = correct)
    colnames(rho) <- items
    rownames(rho) <- items
    return(rho)
  } else {
    res <- .poly_df_full(x, correct = correct)
    dimnames(res[['rho']])  <- list(items, items)
    dimnames(res[['pval']]) <- list(items, items)
    names(res[['tau']]) <- items
    return(res)
  }
}
