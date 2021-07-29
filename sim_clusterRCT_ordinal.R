
pacman::p_load(plyr, ordinal, pbapply, assertthat)


# for a logistic or ordinal regression, the level 1 variance = p1^2/3
# See Agresti A. Categorical data analysis. 2. Hoboken, NJ: Wiley; 2002.



# The ordinal model is parameterised by cumulative response probabilities:



sim_cRCT <- function(k1, # number of clusters in ctr
                     k2, # numbner of clusters in the int
                     m, # average cluster size
                     p, # vector of control probabilities of each outcome level
                     or, # intervention effect as an odds ratio
                     rho) { # intra-class correlation

  assertthat::assert_that(sum(p) == 1)
  # random effects
  sigma2 <- (pi^2) / 3
  theta <- sqrt((rho * sigma2) / (1 - rho)) # TODO check this

  x <- rep(c(0, 1), times = c(k1 * m, k2 * m))
  site <- c(rep(1:(k1 + k2), each = m))
  a <- rnorm((k1 + k2), mean = 0, theta) # site level random effect
  arep <- rep(a, each = m)

  # for 3 categories we have 2 intercepts
  beta0 <- log((p[2] + p[3]) / p[1])
  beta2 <- log(p[3] / (p[1] + p[2]))
  beta1 <- log(or)
  lp1 <- beta0 + beta1 * x + arep
  lp2 <- beta2 + beta1 * x + arep
  p1 <- 1 - (exp(lp1) / (1 + exp(lp1)))
  p3 <- exp(lp2) / (1 + exp(lp2))
  p2 <- (exp(lp1) / (1 + exp(lp1))) - p3
  probs <- cbind(p1, p2, p3)
  y <- apply(probs, 1, function(x) sample(1:3, 1, prob = x))
  fit <- clmm(factor(y) ~ x + (1 | site)) # nolint
  p <- summary(fit)$coefficients[3, 4] # nolint
  # return(data.frame(y, x, site, probs))
  return(p)
}


# now the power calculations
nsims <- 1000
fit_all <- pbreplicate(
  nsims,
  sim_cRCT(
    k1 = 6, k2 = 6, m = 38,
    p = c(.3, .4, .3), or = 2.5, rho = 0.05
  )
)

sum(is.na(fit_all)) / nsims
mean(fit_all < 0.05, na.rm = TRUE)


## validate agains the popower command in Hmisc
n <- 648
# clusters
k <- 12
# intra class correlation
icc <- 0.05

# Design effect
deff <- 1 + icc * ((n / k) - 1)
# effective sample size
ess <- n / deff
# power with ess
p1 <- c(0.525, 0.3, 0.175)
p2 <- pomodm(p = p, odds.ratio = 2.5)
pavg <- (p + p2) / 2
# check wth average powers
popower(pavg, 2.5, n = ceiling(ess), alpha = 0.05)
