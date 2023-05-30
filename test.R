library(microbenchmark)
library(mvtnorm)

# UTILITIES ----
(x <- c(1, 3, 2, 3, 5, 9, 6, 6, 1))
rank(x)
rank_vec(x)
# polychoric::rcppeigen_bothproducts(x)

Phi(x[1]/10)
pnorm(x[1]/10)
qnorm(x/10)
Phi_inv_vec(x/10)

length(unique(x))
n_unique(x)

(y <- sample(c(TRUE, FALSE), 13, replace = TRUE))
as.integer(y)
correct_data(x)
correct_data(y)

n <- 2e3
x <- rnorm(n)
bench_rank <- microbenchmark(
  base = rank(x),
  rcpp = rank_vec(x),
  times = 20000,
  control = list(warmup = 1000)
)
bench_rank

bench_unique <- microbenchmark(
  base = sapply(gss, \(x) length(unique(x))),
  rcpp = sapply(gss, n_unique),
  times = 40000,
  control = list(warmup = 1000)
)
bench_unique

(x <- sample(c(TRUE, FALSE), 1500L, replace = TRUE))
bench_correct <- microbenchmark(
  base = sapply(gss_num, \(x) if (min(x) <= 0) x + (1 - min(x)) else x),
  rcpp = sapply(gss, correct_data),
  times = 40000,
  control = list(warmup = 1000)
)
bench_correct

# MOMENTS ----

## Univariate ----
n <- 2e3
x <- rnorm(n)

bench_sd <- microbenchmark(
  base = sd(x),
  rcpp = sd_(x),
  eign = SD(x),
  times = 20000,
  control = list(warmup = 1000)
)
bench_sd

## Bivariate ----
(z1 <- rnorm(13))
(z2 <- rnorm(13, 1, 1.2))
mvtnorm::dmvnorm(c(z1[1], z2[1]))
phi2(z1[1], z2[1])
mvtnorm::pmvnorm(upper = c(z1[1], z2[1]))
Phi2(z1[1], z2[1])

y <- rnorm(n)

bench_cov <- microbenchmark(
  base = cov(x, y),
  eign = Cov(x, y),
  times = 20000,
  control = list(warmup = 1000)
)
bench_cov

bench_pearson <- microbenchmark(
  base = cor(x, y),
  eign = corr_pearson(x, y),
  times = 20000,
  control = list(warmup = 1000)
)
bench_pearson

(x <- rnorm(30000))
(y <- rchisq(30000, 4))
cor(x, y, method = 'spearman')
rank(x)
rank(y)
rank(x) - rank(y)
cor_spearman(x, y)

bench_spearman <- microbenchmark(
  base = cor(x, y, method = 'spearman'),
  eign = cor_spearman(x, y),
  times = 20000,
  control = list(warmup = 1000)
)
bench_spearman

## Multivariate ----
X <- rmvnorm(n, mean = c(100, 10, 13, -10))

bench_cov_matrix <- microbenchmark(
  base = cov(X),
  eign = cov_matrix(X),
  times = 20000,
  control = list(warmup = 1000)
)
bench_cov_matrix

bench_cor_matrix_pearson <- microbenchmark(
  base = cor(X),
  eign = corr_matrix_pearson(X),
  times = 20000,
  control = list(warmup = 1000)
)
bench_cor_matrix_pearson

rm(x, y)
rm(list = grep('^bench', ls(), value = TRUE))
