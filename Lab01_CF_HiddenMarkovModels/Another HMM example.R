library(mhsmm)
K <- 2 # number of clusters
# Initialize the HMM parameters
start.val <- hmmspec(init = c(1,rep(0,K-1)),
              trans = matrix(c(.8,.2,.2,.8),
                             ncol=K,nrow = K,byrow=TRUE),
                dens.emission = dnorm.hsmm,
              parms.emission = list(mu=c(0,0),sigma=c(1,3)))
mod.k2 <- hmmfit(x = ftse$V1,
                start.val = start.val,mstep = mstep.norm)
plot(mod.k2$loglik, type = "b", 
     ylab = "Log-likelihood", xlab = "Iteration")
clusters <- mod.k2$yhat
plot(ftse$V1,col=clusters)
table(clusters)
mod.k2$model

K <- 3 # number of clusters
# Initialize the HMM parameters
start.val <- hmmspec(init = c(1,rep(0,K-1)),
                     trans = matrix(c(.8,.1,.1,.1,.7,.2,
                                      .1,.1,.8),
                                    ncol=K,nrow = K,byrow=TRUE),
                     dens.emission = dnorm.hsmm,
                     parms.emission = list(mu=c(0,0,0),
                                           sigma=c(1,3,5)))
mod.k3 <- hmmfit(x = ftse$V1,
                 start.val = start.val,mstep = mstep.norm)
plot(mod.k3$loglik, type = "b", 
     ylab = "Log-likelihood", xlab = "Iteration")
clusters <- mod.k3$yhat
plot(ftse$V1,col=clusters)
table(clusters)
mod.k3$model

npar <- NULL
for(k in 2:3)
{
  npar <- c(npar,(k-1) + k*(k-1) + 2*k) # inital + tpm + emission
}
npar
AIC <- c(-2*max(mod.k2$loglik) + 2*npar[1], -2*max(mod.k3$loglik) + 2*npar[2])
BIC <- c(-2*max(mod.k2$loglik) + log(2112)*npar[1], -2*max(mod.k3$loglik) + log(2112)*npar[2])

library(mclust)
fm2 <- Mclust(ftse$V1,G=2,modelNames = "V")
K <- 2 # number of clusters
# Initialize the HMM parameters
start.val.mclust <- hmmspec(init = c(1,rep(0,K-1)),
                     trans = transit(fm2$classification,2),
                     dens.emission = dnorm.hsmm,
                     parms.emission = list(mu=fm2$parameters$mean,
                                           sigma=fm2$parameters$variance$sigmasq))
mod.k2.mclust <- hmmfit(x = ftse$V1,
                 start.val = start.val.mclust,mstep = mstep.norm)
plot(mod.k2.mclust$loglik, type = "b", 
     ylab = "Log-likelihood", xlab = "Iteration")
clusters <- mod.k2.mclust$yhat
plot(ftse$V1,col=clusters)
table(clusters)
mod.k2.mclust$model

#############################################
## Compute a Transition probability matrix ##
## given a vector of states                ##
#############################################

transit <- function(clas, J){
  
  # clas is a vector of integers
  
  if(J<max(clas))
    stop('J must be greater or equal to max(clas)')
  
  T <- length(clas)
  
  if(J==1)
    PI <- matrix(1,1,1)
  
  if(J>1){
    
    PI <- matrix(0,nrow=J,ncol=J)
    
    for(t in 2:T){
      PI[clas[t-1],clas[t]] = PI[clas[t-1],clas[t]] + 1;
    }
    PI <- diag(1/rowSums(PI)) %*% PI
    
  }  
  
  return(PI)
  
}

#
# Use the t distribution
#
library(gamlss)
dtf.hsmm <- function(x, j, model)
{
  ret = dTF(x = x, mu = model$parms.emission$mu[j],
            sigma = model$parms.emission$sigma[j],
            nu = sqrt(model$parms.emission$nu[j]))
  ret[is.na(ret)] = 1
  ret
}
mstep.tf <- function (x, wt) 
{
  k = ncol(wt)
  mu = numeric(k)
  sigma = numeric(k)
  nu = numeric(k)
  for (i in 1:k) {
    tmp <- gamlssML(x,weights = wt[,i],family = "TF")
    mu[i] <- tmp$mu
    sigma[i] <- tmp$sigma
    nu[i] <- tmp$nu
  }
  list(mu = mu, sigma = sigma, nu = nu)
}

#
# Multivariate Gaussian
K <- 2
start.val <- hmmspec(init = c(1,0),
                     trans = matrix(c(.9,.1,.1,.9),byrow=T, nrow = K, ncol = K),
                     parms.emis = list(mu = list(c(0, 0.1, 0.3, -0.1),c(.1,.2, .3,.1)), 
                                       sigma=list(diag(4),diag(4))),
                     dens.emis = dmvnorm.hsmm)
mod.hmm.k2 <- hmmfit(matrix(unlist(dax4[,1:4]),ncol=4), 
                     start.val, mstep = mstep.mvnorm)
