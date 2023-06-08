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

#### Test correlations ----
cor_polychoric(df)
cor_polychoric(df, coef.only = FALSE)
cor_polychoric(df$t1, df$t2, coef.only = FALSE)

## POLYSERIAL CORRELATION ----

#### Create Dataset ----
n <- 200
x <- rnorm(n, 100, 11)
y <- (x - 100) / 11 + rnorm(n)
C <- cut(x, breaks = 5, ordered_result = TRUE)
D <- cut(y, breaks = 5, ordered_result = TRUE)

#### Test correlations ----
cor_polyserial(x, D)
cor_polyserial(x, D, coef.only = FALSE)
cor_polyserial(x, C) # currently unable to handle perfect correlation
cor_polyserial(y, C)
cor_polyserial(y, D) # currently unable to handle perfect correlation

x[sample(1:n, 6)] <- NaN
D[sample(1:n, 6)] <- NA
cor_polyserial(x, D)
