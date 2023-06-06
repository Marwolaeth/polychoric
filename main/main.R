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
polycorr(df)
polycorr(df, coef.only = FALSE)

#### Missing values ----
df_miss <- df
mask <- matrix(
  sample(c(TRUE, FALSE), n*ncol(df), replace = TRUE, prob = c(.9, .1)),
  nrow = n
)
df_miss[!mask] <- NA
summary(df_miss)
polycorr(df_miss)
polycorr(df_miss, coef.only = FALSE)

#### Correlate + continuous ----
df$c1 <- rlnorm(n, 4, .4)
head(df, 13)
polycorr(df)

## GSS 2012 SCHWARTZ VALUES MODULE ----
data("gss12_values", package = 'polychoric')

#### A pair of discrete vectors ----
polycorr(gss12_values$valorig, gss12_values$valeql)
polycorr(gss12_values$valorig, gss12_values$valeql, coef.only = FALSE)

#### A contingency table ----
(G <- table(gss12_values$valorig, gss12_values$valeql))
polycorr(G)

#### A data.frame ----
polycorr(gss12_values)
.poly_df_full(gss12_values)

#### Mixed variable types ----
# For safety, returns Spearman's rho if at least one of the vectors
# is presumably continuous (n_distinct(x) > 10 | n_distinct(y) > 10)
x <- rnorm(nrow(gss12_values))
polycorr(gss12_values$valspl, x)

## IN LIBRARY ----
library(polychoric)
?polycorr

## BENCHMARK ----
library(microbenchmark)
library(psych)
gss_num <- gss12_values |> lapply(as.integer) |> as.data.frame()
head(gss_num, 13)
gc()

bm <- microbenchmark(
  standard = polychoric(gss_num),
  polycorr = polycorr(gss12_values, coef.only = FALSE),
  times = 32L,
  control = list(warmup = 6)
)
bm

bm <- microbenchmark(
  standard = polychoric(table(gss12_values$valorig, gss12_values$valeql)),
  polycorr = polycorr(gss12_values$valorig, gss12_values$valeql),
  times = 32L,
  control = list(warmup = 6)
)
bm
