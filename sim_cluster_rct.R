# Program: sim_cluster_RCT.R
# Author: Chirs Oldmeadow
# Date: 15-12-2020
# Summary: Some functions for estimating power in a cluster RCT with a
#          continuous outcome through simulations. This version allows for
#          baseline adjustment ANCOVA style


library(plyr)
library(lme4)
library(lmerTest)
library(pbapply)
library(here)

set.seed(101)


# function to create the design matrix
 create_design <- function(m,k1,k2,m1,s){
     # m = number of student per school
     # k = number of schools per arm
     # s = standard deviation of outcome (error sd)  
  
    int <- expand.grid(iid = c(1:m),  site = c(1:k1), ttt = 0)
    ctr <- expand.grid(iid = c(1:m),  site = c((k1+1):(k1+k2)), ttt = 1)
    
    expdat <- rbind(int,ctr)
    expdat$site <- as.factor(expdat$site)
    expdat$baseline <- rnorm(n = ((k1+k2)*m), mean = m1, sd = s) 
    return(expdat)  
 }


# function to simulate the outcome for a given design matrix, effect size, and
# other charactertistcs of the data (baseline followup correlation and ICC)  
 sim_outcomes_cRCT<-function(design,m1,m2,s,r01,rho,nreps){
    #m1 = mean in the control group
    #m2 = mean in the intervention group
    # s = error sd
    # r01 = pre post correlation 
    b1 <- m2-m1
    b2 <- r01  # assume sd followup = sd pre
    beta <- c(m1, b1, b2)
    names(beta)<-c("(Intercept)","ttt", "baseline")
    sigma2 <- s*s 
    theta <-sqrt((rho*sigma2)/(1-rho))/s
    
    names(theta)<-c("site.(Intercept)")
    ss <- simulate(~ ttt + baseline + (1 | site), nsim = nreps, 
                   family = gaussian,
                   newdata = design, 
                   newparams = list(theta = theta,sigma = s,   beta = beta))
    return(ss)
}


# fits a linear mixed model to the simulated data

 fitsim <- function(i,dat = expdat,out = ss) {
# dat is the data structure
# out is the list of simulated outcomes

   return <- tryCatch({
     dat$y <- out[[i]]
     res<-coef(summary(lmerTest::lmer(y ~ ttt + baseline + (1|site), data = dat)))["ttt", ]
     names(res)<-c("est","se","df","t","p")
     return(res)
   },
   warning =function(e) {
     #message(e)  # print error message
     return(c(est=NA, se=NA, df = NA, t = NA, p=NA))
   },
   
   error=function(e) {
     #message(e)  # print error message
     return(c(est=NA, se=NA, df = NA, t=NA, p=NA))
   })
   
   return(return)
 }



 # now the power calculations
expdat <- create_design(m=10,k1 = 10, k2 = 10, m1 = 80, s = 90)
ss <- sim_outcomes_cRCT(expdat, m1 = 80 ,m2 = 160 ,s = 90,r01 = 0.6, rho = .2,nreps = 1000)
fitAll <- t(pbsapply(seq(1000), function(i) fitsim(i)))

pow <- mean( fitAll[,"p"] < 0.05, na.rm = TRUE)
pow    

save.image(file=here('power.RData'))

# load(file = here('power.RData'))






