## MIXED TYPE DATA CORRELATION ----

#### Create Dataset ----
likert5_lvls <- c(
  'Disagree',
  'Moderately Disagree',
  'Neutral',
  'Moderately Agree',
  'Agree'
)
(likert5 <- factor(likert5_lvls, ordered = TRUE, levels = likert5_lvls))

likert6_lvls <- c(
  'Not like me at all',
  'Not like me',
  'Probably not like me',
  'Sligtly like me',
  'Like me',
  'Definitely like me'
)
(likert6 <- factor(likert6_lvls, ordered = TRUE, levels = likert6_lvls))

likert2_lvls <- c(
  'No',
  'Yes'
)
(likert2 <- factor(likert2_lvls, ordered = TRUE, levels = likert2_lvls))

(n <- 200)
df <- data.frame(
  f1 = sample(likert5, n, replace = TRUE),
  s1 = sample(likert6, n, replace = TRUE),
  t1 = sample(likert2, n, replace = TRUE),
  f2 = sample(likert5, n, replace = TRUE),
  s2 = sample(likert6, n, replace = TRUE),
  t2 = sample(likert2, n, replace = TRUE),
  f3 = sample(likert5, n, replace = TRUE),
  s3 = sample(likert6, n, replace = TRUE),
  t3 = sample(likert2, n, replace = TRUE),
  l1 = sample(c(TRUE, FALSE), n, replace = TRUE)
)
head(df, 13)

#### Correlate ----
cor_polychoric(df)
cor_polychoric(df, coef.only = FALSE)

#### Missing values ----
df_miss <- df
mask <- matrix(
  sample(c(TRUE, FALSE), n*ncol(df), replace = TRUE, prob = c(.9, .1)),
  nrow = n
)
df_miss[!mask] <- NA
summary(df_miss)
cor_polychoric(df_miss)
cor_polychoric(df_miss, coef.only = FALSE)

#### Correlate + continuous ----
df$c1 <- rlnorm(n, 4, .4)
head(df, 13)
cor_polychoric(df)

## POLYSERIAL CORRELATION ----

#### Create Dataset ----
library(mvtnorm)
(rho <- .77)    # True rho
(s1 <- 11)      # True X standard deviation
(sigma <- matrix(c(s1^2, s1*rho, s1*rho, 1), ncol = 2))
df <- rmvnorm(200, mean = c(100, 0), sigma = sigma) |>
  as.data.frame() |>
  setNames(c('x', 'y'))
head(df, 13)
psych::corr.test(df)
df$c <- cut(df$x, breaks = 5, labels = likert5_lvls, ordered_result = TRUE)
df$d <- cut(df$y, breaks = 5, labels = likert5_lvls, ordered_result = TRUE)
summary(df)

#### Test correlations ----
cor_polyserial(df$x, df$d)
cor_polyserial(df$x, df$c)
cor_polyserial(df$y, df$d)
psych::polyserial(df$x, df[,3:4] |> lapply(as.integer) |> as.data.frame())

bm_ps <- microbenchmark(
  standard = polyserial(df$x, df[,3:4] |> lapply(as.integer) |> as.data.frame()),
  polyseri = cor_polyserial_full(df$x, df$d),
  times = 26L,
  control = list(warmup = 3L)
)
bm_ps

x <- df$x
d <- df$d
x[sample(1:200, 6)] <- NaN
d[sample(1:200, 6)] <- NA

cor_polyserial(x, d)
cor_polyserial(x, d, coef.only = FALSE)

## GSS 2012 SCHWARTZ VALUES MODULE ----
data("gss12_values", package = 'polychoric')

#### A pair of discrete vectors ----
cor_polychoric(gss12_values$valorig, gss12_values$valeql)
cor_polychoric(gss12_values$valorig, gss12_values$valeql, coef.only = FALSE)

#### A contingency table ----
(G <- table(gss12_values$valorig, gss12_values$valeql))
cor_polychoric(G)

## IN LIBRARY ----
library(polychoric)
?cor_polychoric
?cor_polyserial

## BENCHMARK ----
library(microbenchmark)
library(psych)
gss_num <- gss12_values |> lapply(as.integer) |> as.data.frame()
head(gss_num, 13)
gc()

bm <- microbenchmark(
  standard = polychoric(gss_num),
  cor_polychoric = cor_polychoric(gss12_values, coef.only = FALSE),
  times = 32L,
  control = list(warmup = 6)
)
bm

bm_xy <- microbenchmark(
  standard = polychoric(table(gss12_values$valorig, gss12_values$valeql)),
  cor_polychoric = cor_polychoric(gss12_values$valorig, gss12_values$valeql),
  times = 32L,
  control = list(warmup = 6)
)
bm_xy
