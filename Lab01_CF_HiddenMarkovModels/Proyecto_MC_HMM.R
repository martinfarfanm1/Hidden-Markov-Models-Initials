library(mhsmm)
require(tidyverse)
require(gamlss)
require(mclust)
require(pROC)
library(dplyr)

cont <- Rieti[, 1:9]

# load(AirPassengers)
# AirPassengers %>% plot
# #always use this to see cyclical data
# plot(decompose(AirPassengers)) # obs is the sum of these three things.

plot(cont)


# MClust ------------------------------------------------------------------

# Analisis preliminar
mod1.mc <- Mclust(cont)
mod1.mc #VVE, 6 Clusters

mod1.mc$BIC
bic.mod.mc <- mod1.mc$bic

mod1.mc$parameters
diag(mod1.mc$parameters$variance$sigma[, , 1])
mod1.mc$G
mod1.mc %>% summary()

mod1.mc %>% plot
mod1.mc$classification
mod1.mc$classification %>% table

# Levels of contamination on time
plot(x = cont$NO2,
     y = cont$Oxylene,
     col = mod1.mc$classification)

# Demonstrate seasonality
plot(cont$Toluene, col = mod1.mc$classification)
plot(cont$O3, col = mod1.mc$classification)
plot(cont$SO2, col = mod1.mc$classification)

seq.values <- list()

for (j in 1:ncol(cont)) {
  denf <- 0
  den <- list()
  pi <- list()
  seq.values[[j]] <-
    seq(min(cont[[j]]), max(cont[[j]]), length = 2000)
  for (i in 1:mod1.mc$G) {
    # recorrer las 6 clusters
    den[[i]] <- dnorm(
      seq.values[[j]],
      mean = mod1.mc$parameters$mean[j, i],
      sd = mod1.mc$parameters$variance$sigma[j, j, i]
    )
    
    pi[[i]] <- mod1.mc$parameters$pro[i]
    denf <- denf + (den[[i]] * pi[[i]])
  }
  hist(
    cont[[j]],
    breaks = 20,
    freq = FALSE,
    main = paste("Histogram of", names(cont[j])),
    xlab = names(cont[j]),
    ylim = c(0, max(denf))
  )
  lines(seq.values[[j]], denf, lwd = 3)
  for (i in 1:mod1.mc$G) {
    lines(seq.values[[j]], den[[i]] * pi[[i]], col = i + 1, lwd = 2)
  }
}

plot(cont, col = mod1.mc$classification)

# VARIACIONES DEL MODELO ORIGINAL

# cHANGINNG the # clusters // Sensitivity analysis
G <- 10
mod1.mc.G <- list()
mod1.mc.G.BIC <- list()

for (i in 1:G) {
  mod1.mc.G[[i]] <- Mclust(cont, G = i)
  mod1.mc.G.BIC[[i]] <- mod1.mc.G[[i]]$bic
}

# Maxima accuracy
mod1.mc.G.BIC %>% unlist() %>% max()

# How many K (neighbours) we need to achive it
mod1.mc.G.BIC %>% unlist() %>% which.max() %>% as.numeric

# this shows that the algorithm was right in the first attempt

# HMM ---------------------------------------------------------------------

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

# since we determined the number of clusters is 6
K <- mod1.mc$G

init.parm.mu <- function(mats, k) {
  mus <- list()
  for (i in 1:k) {
    mus[[i]] <- mats[, i]
  }
  return(mus)
}

init.parm.sigma <- function(mats, k) {
  sig <- list()
  for (i in 1:k) {
    sig[[i]] <- mats[, , i]
  }
  return(sig)
}

start.val.6 <- hmmspec(
  init = c(1, rep(0, K - 1)),
  trans = transit(mod1.mc$classification, K),
  parms.emission = list(
    mu = init.parm.mu(mod1.mc$parameters$mean, K),
    sigma = init.parm.sigma(mod1.mc$parameters$variance$sigma, K)
  ),
  dens.emission = dmvnorm.hsmm
)

mod.hmm.mc6 <-
  hmmfit(matrix(unlist(cont), ncol = ncol(cont)), start.val.6, mstep = mstep.mvnorm)

plot(
  mod.hmm.mc6$loglik,
  type = "b",
  ylab = "Log-likelihood",
  xlab = "Iteration",
  main = "Convergence K=6 // MClust"
)

clust.hmm.mc6 <-
  mod.hmm.mc6$yhat #Secuencia de 1:k como se clusterizo

# Demonstrate seasonality
plot(cont$Toluene, col = clust.hmm.mc6)

# Levels of contamination on time
plot(x = cont$Etilbenzene,
     y = cont$Oxylene,
     col = clust.hmm.mc6)


plot(cont, col = clust.hmm.mc6)

clust.hmm.mc6 %>% table

mod.hmm.mc6$model

conf.mat <-
  table(mod1.mc$classification, clust.hmm.mc6) # truth table // confussion matrix

# Metrics of the categorization

accuracy <- sum(diag(conf.mat)) / nrow(cont)

# for doing a confussion matrix and a lot of other things, package caret
library(caret)

conf.mat.caret <- confusionMatrix(data = as.factor(clust.hmm.mc6),
                                  reference = as.factor(mod1.mc$classification))

conf.mat.caret
conf.mat.caret$overall["Accuracy"]
conf.mat.caret$table
conf.mat.caret$byClass

# ROC
# create data partition para el otro

# split = createDataPartition(data$default, p=0.81, list=FALSE)
# train_data = data[split, ]
# test_data = data[-split, ]

roc.mod <- roc(
  mod1.mc$classification,
  # Valores reales
  clust.hmm.mc6,
  #valores reales
  plot = TRUE,
  legacy.axes = TRUE,
  col = "midnightblue",
  lwd = 3,
  auc.polygon = T,
  auc.polygon.col = "lightblue",
  print.auc = T
)

roc.mod$auc

# How to do model selection? According to BIC

num_term <- function(clust) {
  kn <-
    (clust - 1) + clust * (clust - 1) + 2 * clust #inital + tpm + emission (observations)
  return(kn)
}

n.term <- num_term(K)

# (AIC <- -2 * max(mod.hmm.mc6$loglik) + 2 * n.term)
(BIC <- -2 * max(mod.hmm.mc6$loglik) + log(365) * n.term)
# -2 and + is minimize
# the solution with more parameters "is always/usually better", Smallest absolute value

(prop <- 100 - (BIC / bic.mod.mc) * 100)

# ENVIROMENTAL  -----------------------------------------------------------

mod.mc.tot <- Mclust(Rieti)
mod.mc.tot$BIC
mod.mc.tot$bic

mod.mc.tot$parameters
mod.mc.tot$G
mod.mc.tot %>% summary()
mod.mc.tot %>% plot
mod.mc.tot$classification
mod.mc.tot$classification %>% table

plot(Rieti, col = mod.mc.tot$classification)


K.tot <- mod.mc.tot$G

start.val.tot <- hmmspec(
  init = c(1, rep(0, K.tot - 1)),
  trans = transit(mod.mc.tot$classification, K.tot),
  parms.emission = list(
    mu = init.parm.mu(mod.mc.tot$parameters$mean, K.tot),
    sigma = init.parm.sigma(mod.mc.tot$parameters$variance$sigma, K.tot)
  ),
  dens.emission = dmvnorm.hsmm
)

mod.hmm.mc.tot <-
  hmmfit(matrix(unlist(Rieti), ncol = ncol(Rieti)), start.val.tot,
         mstep = mstep.mvnorm)

plot(
  mod.hmm.mc.tot$loglik,
  type = "b",
  ylab = "Log-likelihood",
  xlab = "Iteration",
  main = "Convergence K=6 // MClust"
)

clust.hmm.mc.tot <-
  mod.hmm.mc.tot$yhat #Secuencia de 1:k como se clusterizo
plot(Rieti, col = clust.hmm.mc.tot)

clust.hmm.mc.tot %>% table

mod.hmm.mc.tot$model

conf.mat.caret.tot <-
  confusionMatrix(
    data = as.factor(clust.hmm.mc.tot),
    reference = as.factor(mod.mc.tot$classification)
  )

conf.mat.caret.tot
conf.mat.caret.tot$overall["Accuracy"]
conf.mat.caret.tot$table
conf.mat.caret.tot$byClass

roc.mod.tot <- roc(
  mod.mc.tot$classification,
  # Valores reales
  clust.hmm.mc.tot,
  #valores reales
  plot = TRUE,
  legacy.axes = TRUE,
  col = "midnightblue",
  lwd = 3,
  auc.polygon = T,
  auc.polygon.col = "lightblue",
  print.auc = T
)

roc.mod.tot$auc

n.term.tot <- num_term(K.tot)

# (AIC <- -2 * max(mod.hmm.mc6$loglik) + 2 * n.term)
(BIC.tot <- -2 * max(mod.hmm.mc.tot$loglik) + log(365) * n.term.tot)

(prop.tot <- 100 - (abs(BIC.tot) / abs(mod.mc.tot$bic) * 100))

#without enviromental data
bic.mod.mc #mc
BIC #hmm

# with
mod.mc.tot$bic
BIC.tot

#According to BIC, HMM always outperforms the MClust
