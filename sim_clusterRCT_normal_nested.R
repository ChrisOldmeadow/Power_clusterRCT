library(plyr)
library(lme4)
library(MASS) ## for glmmPQL




    set.seed(101)

pow_simclusterRCT<-function(m,t,k1,k2,m1,m2,s,rho,alpha,nreps){

  # m = number of student per class
  # t = number of classes per school
  # k = number of schools per arm
  
  
    int <- expand.grid(iid = c(1:m), class = c(1:t), site = c(1:k1), ttt = 0)
    ctr <- expand.grid(iid = c(1:m), class = c(1:t),  site = c((k1+1):(k1+k2)), ttt = 1)
    
    expdat <- rbind(int,ctr)
    
    
    #expdat$site <- factor(ifelse(expdat$ttt == 0,expdat$site, expdat$site + k))
    
    
    
    b1 <- m2-m1
  

    beta <- c(m1, b1)
    names(beta)<-c("(Intercept)","ttt")
    
   sigma2 <- s*s 
    theta <-sqrt((rho*sigma2)/(1-rho))
    
    
    names(theta)<-c("site.(Intercept)")
    
    
    ss <- simulate(~ttt  + (1 | site), nsim = nreps, 
                   family = gaussian,
                   newdata = expdat, 
                   newparams = list(theta = theta,sigma = sigma2,   beta = beta))
    
    
    expdat$resp <- ss[, 1]
    fit1 <- lme4::lmer(resp ~ ttt + (1 | site),  data=expdat )
    
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
    
    fitAll <- ldply(seq(nreps), function(i) fitsim(i))
    
    pow<-with(fitAll, mean( p < 0.05, na.rm = TRUE))
    
    
    return(pow)

}

pow_simclusterRCT(m=50, t = 4,k1=12,k2=12,m1=1,m2=1.3,s = 1,rho=.05,alpha=.025,nreps = 1500)








