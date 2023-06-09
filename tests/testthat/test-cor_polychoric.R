## DATA ----
#### Create Labels ----
likert5_lvls <- c(
  'Disagree',
  'Moderately Disagree',
  'Neutral',
  'Moderately Agree',
  'Agree'
)
likert5 <- factor(likert5_lvls, ordered = TRUE, levels = likert5_lvls)

likert6_lvls <- c(
  'Not like me at all',
  'Not like me',
  'Probably not like me',
  'Sligtly like me',
  'Like me',
  'Definitely like me'
)
likert6 <- factor(likert6_lvls, ordered = TRUE, levels = likert6_lvls)

likert2_lvls <- c(
  'No',
  'Yes'
)
likert2 <- factor(likert2_lvls, ordered = TRUE, levels = likert2_lvls)

#### Create Dataset ----
n <- 200
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

#### Missing values ----
df_miss <- df
mask <- matrix(
  sample(c(TRUE, FALSE), n*ncol(df), replace = TRUE, prob = c(.9, .1)),
  nrow = n
)
df_miss[!mask] <- NA

## TESTS ----
test_that("cor_polychoric() warns when finds zero values", {
  expect_warning(cor_polychoric(df$l1, df$f1))
})

test_that("cor_polychoric() warns when suspeсts continuous input", {
  expect_warning(cor_polychoric(df$f1, rbinom(n, 20, .5)))
})

test_that("cor_polychoric() handles missing values", {
  expect_length(cor_polychoric(df_miss[,1:9]), 9*9)
})

test_that("cor_polychoric(coef.only=FALSE) returns an object of class `polychoric`", {
  expect_s3_class(cor_polychoric(df[,1:9], coef.only = FALSE), class = 'polychoric')
})

test_that("cor_polychoric() is fast", {
  expect_lt(
    system.time(cor_polychoric(df[,1:9], coef.only = FALSE))[['user.self']],
    .2
  )
})
