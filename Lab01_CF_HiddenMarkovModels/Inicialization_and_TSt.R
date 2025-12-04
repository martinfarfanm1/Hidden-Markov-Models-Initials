library(mhsmm)
require(tidyverse)
require(gamlss)

ftse <-
  read.csv(url("http://www.hmms-for-time-series.de/second/data/ftse.txt"),
           col.names = "V1")

plot(ftse$V1, type = "l")

# COMO ES COUNTS, VEO 2-3 PROMEDIOS ERGO 2-3 POISSONS
histDist(ftse$V1, "NO")

# Financial almost never follows Gaussian.

K <- 2 # Number of clusters
# Initalialize the HMM parameters

# Seccion 4.2.4

start.val <- hmmspec(
  init = c(1, rep(0, K - 1)),
  #Initial probabilities by 4.2.4
  trans = matrix(
    c(.8, .2, .2, .8)
    ,
    nrow = K,
    ncol = K,
    byrow = T
  ),
  #initial values for the transition matrix
  parms.emission  = list(mu = c(0, 0), sigma = c(1, 3)),
  #depends on type of density // needs to be a list // num of variables and each one repeated because we have two clusters
  dens.emission  = dnorm.hsmm
)
start.val

mod.k2 <- hmmfit(ftse$V1, start.val, mstep = mstep.norm)

# mstep.norm depende del x y de los weighted probabilites posterior probability
# wt is a matrix de las probabilidades para pasar w_tk = Pz(Ct = k|X)

plot(
  mod.k2$loglik,
  type = "b",
  ylab = "Log-likelihood",
  xlab = "Iteration",
  main = "Convergence K=2"
)

clusters <- mod.k2$yhat #Secuencia de 1 y 2, como se clusterizo
plot(ftse$V1, col = clusters)
# captures volatility
clusters %>% table

mod.k2$model
# 1st cluster is very persistent, hay 99 de prob que me quede ahí
# mu de 1st es muy cercano a 0, pero en el otro es alto, y negativo (pierdo plata)
# sigma es la volatilidad, segunda cluster is super volatile

# 3 clusters

K <- 3 # Number of clusters
# Initalialize the HMM parameters

# Seccion 4.2.4

start.val3 <- hmmspec(
  init = c(1, rep(0, K - 1)),
  #Initial probabilities by 4.2.4
  trans = matrix(
    c(.8, .1, .1, .1, .8, .1, .1, .1, .8)
    ,
    nrow = K,
    ncol = K,
    byrow = T
  ),
  #initial values for the transition matrix
  parms.emission  = list(mu = c(0, 0, 0), sigma = c(1, 3, 5)),
  #depends on type of density // needs to be a list // num of variables and each one repeated because we have 3 clusters
  dens.emission  = dnorm.hsmm
)
start.val3

mod.k3 <- hmmfit(ftse$V1, start.val3, mstep = mstep.norm)

# mstep.norm depende del x y de los weighted probabilites posterior probability
# wt is a matrix de las probabilidades para pasar w_tk = Pz(Ct = k|X)

plot(
  mod.k3$loglik,
  type = "b",
  ylab = "Log-likelihood",
  xlab = "Iteration",
  main = "Convergence K=3"
)

clusters3 <- mod.k3$yhat #Secuencia de 1, 3 y 2, como se clusterizo
plot(ftse$V1, col = clusters3)
# captures volatility
clusters3 %>% table

mod.k3$model

# How to do model selection?

npar <- NULL

for (k in 2:3) {
  npar <-
    c(npar, (k - 1) + k * (k - 1) + 2 * k) #inital + tpm + emission (observations)
}
npar

mod.k2$loglik

AIC <-
  c(-2 * max(mod.k2$loglik) + 2 * npar[1],
    -2 * max(mod.k3$loglik) + 2 * npar[2])
BIC <-
  c(-2 * max(mod.k2$loglik) + log(2112) * npar[1],
    -2 * max(mod.k3$loglik) + log(2112) * npar[2])
# -2 and + is minimize
# la solucion con más parámetros "siempre es mejor", Smallest absolute value


# initialization


library(mclust)
fm2 <-
  Mclust(ftse$V1, G = 2, modelNames = "V") #same variance for all classes
fm2$parameters

K <- 2 # Number of clusters
# Initalialize the HMM parameters con McLust


start.val.mc <- hmmspec(
  init = c(1, rep(0, K - 1)),
  #Initial probabilities by 4.2.4
  trans = transit(fm2$classification,2),
  #initial values for the transition matrix
  parms.emission  = list(
    mu = fm2$parameters$mean,
    sigma = fm2$parameters$variance$sigmasq
  ),
  #depends on type of density // needs to be a list // num of variables and each one repeated because we have two clusters
  dens.emission  = dnorm.hsmm
)
start.val.mc

mod.k2.mc <- hmmfit(ftse$V1, start.val.mc, mstep = mstep.norm)

# mstep.norm depende del x y de los weighted probabilites posterior probability
# wt is a matrix de las probabilidades para pasar w_tk = Pz(Ct = k|X)

plot(
  mod.k2.mc$loglik,
  type = "b",
  ylab = "Log-likelihood",
  xlab = "Iteration",
  main = "Convergence K=2 / McLust"
)

clusters.mc <-
  mod.k2.mc$yhat #Secuencia de 1 y 2, como se clusterizo
plot(ftse$V1, col = clusters.mc)
# captures volatility
clusters.mc %>% table

mod.k2.mc$model



# Use the t distribution
dnorm.hsmm
# t
dtf.hsmm <- function(x, j, model) {
  ret <- dTF(
    x = x,
    mu = model$parms.emission$mu[j],
    sigma = sqrt(model$parms.emission$sigma[j]),
    nu = sqrt(sqrt(model$parms.emission$nu[j]))
  )
  ret[is.na(ret)] = 1
  ret
}

mstep.tf <- function(x, wt) {
  k = ncol(wt)
  mu = numeric(k)
  sigma = numeric(k)
  nu = numeric(k)
  for (i in 1:k) {
    tmp <- gamlssML(x, weights = wt[, i], family = "TF")
    mu[i] <- tmp$mu
    sigma[i] <- tmp$sigma
    nu[i] <- tmp$nu
  }
  list(mu = mu, sigma = sigma, nu = nu)
  
}

#For doing the transition matrix
#############################################
## Compute a Transition probability matrix ##
## given a vector of states                ##
#############################################

transit <- function(clas, J) {
  # clas is a vector of integers
  
  if (J < max(clas))
    stop('J must be greater or equal to max(clas)')
  
  T <- length(clas)
  
  if (J == 1)
    PI <- matrix(1, 1, 1)
  
  if (J > 1) {
    PI <- matrix(0, nrow = J, ncol = J)
    
    for (t in 2:T) {
      PI[clas[t - 1], clas[t]] = PI[clas[t - 1], clas[t]] + 1
      
    }
    PI <- diag(1 / rowSums(PI)) %*% PI
    
  }
  
  return(PI)
  
}
