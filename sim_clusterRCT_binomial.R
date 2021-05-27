library(plyr)
library(lme4)
library(MASS) ## for glmmPQL
library(pbapply)





###########################################################################################
## Using glmPQL (better for lower numbers of clusters)
###########################################################################################

designcRCT <- function(m,
                       k1,
                       k2) {
  int <- expand.grid(iid = c(1:m), site = c(1:k1), ttt = 0)
  ctr <- expand.grid(iid = c(1:m), site = c((k1 + 1):(k1 + k2)), ttt = 1)
  df <- rbind(int, ctr)
  return(df)
}



simcRCT <- function(df,
                    nsim,
                    p1, # prev of outcome at pre intervetion
                    p2, # prev of outcome at post intervetion
                    rho) { # intra-class correlation

  o2 <- p2 / (1 - p2)
  o1 <- p1 / (1 - p1)
  beta <- c(log(o2), log(o1 / o2))
  names(beta) <- c("(Intercept)", "ttt")
  sigma2 <- (pi^2) / 3
  theta <- sqrt((rho * sigma2) / (1 - rho))
  names(theta) <- c("site.(Intercept)")
  ss <- simulate(~ ttt + (1 | site),
    nsim = nsim, family = binomial,
    newdata = df, newparams = list(theta = theta, beta = beta)
  )
  return(ss)
}

# fits a linear mixed model to the simulated data

fitsimPQL <- function(i, dat = simdat, out = ss) {
  # dat is the data structure
  # out is the list of simulated outcomes

  return <- tryCatch(
    {
      dat$resp <- out[[i]]
      res <- coef(summary(MASS::glmmPQL(resp ~ ttt,
        random = ~ 1 | site,
        family = binomial,
        data = dat,
        verbose = FALSE
      )))["ttt", ]
      names(res) <- c("est", "se", "df", "t", "p")
      return(res)
    },
    warning = function(e) {
      # message(e)  # print error message
      return(c(est = NA, se = NA, df = NA, t = NA, p = NA))
    },
    error = function(e) {
      # message(e)  # print error message
      return(c(est = NA, se = NA, df = NA, t = NA, p = NA))
    }
  )

  return(return)
}


# now the power calculations


simdat <- designcRCT(
  m = 38,
  k1 = 6,
  k2 = 6
)

ss <- simcRCT(simdat, p1 = .3, p2 = .15, rho = .05, nsim = 10000)

fitAll <- t(pbsapply(seq(1000), function(i) fitsimPQL(i)))

pow <- mean(fitAll[, "p"] < 0.05, na.rm = TRUE)
pow
