library(microbenchmark)
library(mvtnorm)
library(latex2exp)
library(psych)

## POLYSERIAL CORRELATION ----

likert5_lvls <- c(
  'Disagree',
  'Moderately Disagree',
  'Neutral',
  'Moderately Agree',
  'Agree'
)
(likert5 <- factor(likert5_lvls, ordered = TRUE, levels = likert5_lvls))

#### Create Dataset ----
(rho <- .77)    # True rho
(s1 <- 11)      # True X standard deviation
(sigma <- matrix(c(s1^2, s1*rho, s1*rho, 1), ncol = 2))
set.seed(111)
df <- rmvnorm(200, mean = c(100, 0), sigma = sigma) |>
  as.data.frame() |>
  setNames(c('x', 'y'))
head(df, 13)
psych::corr.test(df)
df$c <- cut(df$x, breaks = 5, labels = likert5_lvls, ordered_result = TRUE)
df$d <- cut(df$y, breaks = 5, labels = likert5_lvls, ordered_result = TRUE)
summary(df)

#### Test objective and the derivative ----
(tau <- ((table(df$d) |> cumsum()) / 200) |> qnorm())
(pd <- polyserial_pd(rho, df$x, df$d))
polyserial_loglik(df$x, df$d, pd)

loglik_test <- function(rho) {
  pd <- polyserial_pd(rho, df$x, df$d)
  polyserial_loglik_(df$x, df$d, pd)
}

dl_test <- function(rho) {
  pd <- polyserial_pd(rho, df$x, df$d)
  dl_drho(rho, pd, df$x, df$d)
}

ll <- sapply(rhos, loglik_test)
rhos[which.min(ll)]
psych::polyserial(df$x, df[,3:4] |> lapply(as.integer) |> as.data.frame())
plot(
  rhos, ll,
  type = 'b',
  xaxt = 'n',
  xlim = c(-1, 1),
  xlab = TeX(r'(Candidate $\rho$)'),
  ylab = TeX(r'(\ln L)'),
  main = TeX(r'(Log likelihood of observed data given $\rho$)')
)
axis(1, seq(-1, 1, .25), labels = TRUE, cex.axis = .8)
axis(1, seq(-1, 1, .05), labels = FALSE)
grid(nx = 9)
abline(v = rhos[which.min(ll)], col = 'red3')

dl <- sapply(rhos, dl_test)
plot(rhos, dl, type = 'l', lty = 'dashed', col = 'royalblue')
abline(v = rhos[which.min(abs(dl))], col = 'red3')
rhos[which.min(ll)]
rhos[which.min(abs(dl))]
dl[which.min(abs(dl))]

## RANK ----
mat <- df |> sapply(as.integer)
bm_rank <- microbenchmark(
  base = apply(df, 2, rank),
  eign = apply(mat, 2, rank_),
  times = 1e3,
  control = list(warmup = 400L)
)
bm_rank
