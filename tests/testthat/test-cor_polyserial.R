## DATA ----
#### Generate Data ----
n <- 200
x <- rnorm(n, 100, 11)
y <- (x - 100) / 11 + rnorm(n)
C <- cut(x, breaks = 5, ordered_result = TRUE)
D <- cut(y, breaks = 5, ordered_result = TRUE)
l <- sample(c(TRUE, FALSE), n, replace = TRUE)

#### Missing values ----
xm <- x
ym <- y
cm <- C
dm <- D
xm[sample(1:n, 13)] <- NaN
ym[sample(1:n, 13)] <- NaN
cm[sample(1:n, 13)] <- NA
dm[sample(1:n, 13)] <- NA

## TESTS ----
test_that("cor_polyserial() warns when finds zero values in an ordinal variable", {
  expect_warning(cor_polyserial(x, l))
})

test_that("cor_polyserial() handles missing values", {
  expect_true(!is.na(cor_polyserial(xm, dm)))
})

test_that("cor_polyserial(coef.only=FALSE) returns an object of class `polyserial`", {
  expect_s3_class(cor_polyserial(ym, cm, coef.only = FALSE), class = 'polyserial')
})
