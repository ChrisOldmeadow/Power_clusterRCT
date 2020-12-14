library(plyr)
library(lme4)
library(MASS) ## for glmmPQL





pow_simclusterRCT<-function(m,k1,k2,p1,p2,rho,alpha,nreps){

    int <- expand.grid(iid = c(1:m), site = c(1:k1), ttt = 0)
    ctr <- expand.grid(iid = c(1:m), site = c((k1+1):(k1+k2)), ttt = 1)
    
    expdat <- rbind(int,ctr)
    
    
    #expdat$site <- factor(ifelse(expdat$ttt == 0,expdat$site, expdat$site + k))
    
    
    set.seed(101)
    nsim <- 200
    
    o2<-p2/(1-p2)
    o1<- p1/(1-p1)
    
    beta <- c(log(o2), log(o1/o2))
    names(beta)<-c("(Intercept)","ttt")
    
    
    sigma2 <-(pi ^ 2) / 3
    theta <-sqrt((rho*sigma2)/(1-rho))
    
    
    names(theta)<-c("site.(Intercept)")
    
    
    ss <- simulate(~ttt  + (1 | site), nsim = nsim, family = binomial, 
                   newdata = expdat, newparams = list(theta = theta,   beta = beta))
    
    
    expdat$resp <- ss[, 1]
    fit1 <- lme4::glmer(resp ~ ttt + (1 | site), family = binomial, data=expdat, control=glmerControl(optimizer="bobyqa", check.conv.singular = .makeCC(action = "ignore",  tol = 1e-4)))
    
   # fit1 <- MASS::glmmPQL(resp ~ ttt , random = 1 | site, family = binomial, data=expdat)
    
    fitsim <- function(i) {
      
      return <- tryCatch({
        res<-coef(summary(lme4::refit(fit1, ss[[i]])))["ttt", ]
        names(res)<-c("est","se","z","p")
        return(res)
      },
      warning =function(e) {
        #message(e)  # print error message
        return(c(est=NA, se=NA, z=NA, p=NA))
      },
      
      error=function(e) {
        #message(e)  # print error message
        return(c(est=NA, se=NA, z=NA, p=NA))
      })
      
      return(return)
    }
    
    fitAll <- ldply(seq(nsim), function(i) fitsim(i))
    
    pow<-with(fitAll, mean( p < 0.05, na.rm = TRUE))
    
    
    return(pow)

}

pow_simclusterRCT(m=15,k1=9,k2=11,p1=.08,p2=.15,rho=.02,alpha=.05,nreps = 500)
pow_simclusterRCT(m=15,k1=9,k2=11,p1=.1,p2=.2,rho=.01,alpha=.05,nreps = 500)



#70% response rate of teachers = 120*.7
#126 schools = 42 per arm, 30% of schools droppped out = 32*.7 = 29 schools per arm


pow_simclusterRCT(m=84,k1=29,k2=29,p1=.125,p2=.25,rho=.2,alpha=.025,nreps = 500)




library(lmerTest)

simcRCT<-function(k,m,mu1,mu2,sigma,rho,alpha){
  
  treat <- c(rep(0,(k*m)),rep(1,(k*m)))
  site <- c(rep(1:(k*2),each=m))
  
  
  var_b <-(rho*sigma^2)/(1-rho)
  
  err <- rnorm(k*2*m,mean=0,sigma)
  a <- rnorm(k*2,mean=0,sqrt(var_b)) # site level random effect
  
  arep <- rep(a,each = m)
  
  b1 <- mu2-mu1
  y <- mu1 + b1*treat + arep + err
  
  
  
  res <- lmer(y ~ treat + (1|site))
  summary(res)
  tests<-anova(res,type = 2)
  p<-tests$`Pr(>F)`
}

sims<-replicate(500,simcRCT(k=29,
                              m=85,
                              rho=.2,
                              sigma=1,
                              mu1=1,
                              mu2=1.4))




power<-mean(ifelse(sims<0.025,1,0))
power
