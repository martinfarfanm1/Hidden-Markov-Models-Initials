require(tidyverse)
require(gamlss)
library(mhsmm)

x <-
  read.csv(
    url(
      "http://www.hmms-for-time-series.de/second/data/earthquakes.txt"
    ),
    sep = "",
    col.names = c("Year", "Counts")
  )
x %>% view()
x %>% plot(type = "l")
plot(x$Counts, type = "l")
hist(x$Counts)
# COMO ES COUNTS, USAMOS POISSON
# POISSON ES LAMBDA (MEAN O VARIANZA)
# COMO ES COUNTS, VEO 2-3 PROMEDIOS ERGO 2-3 POISSONS
histDist(x$Counts, "PO")

# 2 Clusters
K <- 2
dpois.hsmm
start_val <- hmmspec(
  init = rep(1 / K, K),
  trans = matrix(
    c(.8, .2, .1, .9),
    nrow = K,
    ncol = K,
    byrow = T
  ),
  parms.emission = list(lambda = c(15, 25)),
  dens.emission = dpois.hsmm
)
# INICIALIZAMOS CON VALORES RANDOMS POR OBSERVACIONES
mod.hmm.k2 <- hmmfit(x[, 2], start_val, mstep = mstep.pois)
plot(mod.hmm.k2$loglik,
     type = "b",
     ylab = "Log-likelihood",
     xlab = "Iteration")
states <- mod.hmm.k2$yhat
plot(x[, 2], col = states, main = "# of earthquakes")
abline(h = mod.hmm.k2$model$parms.emission$lambda[1])
abline(h = mod.hmm.k2$model$parms.emission$lambda[2], col = 2)

mod.hmm.k2 %>% summary
mod.hmm.k2$loglik %>% max()

# 3 Clusters
K3 <- 3
start_val3 <- hmmspec(
  init = rep(1 / K3, K3),
  trans = matrix(
    1 / K3,
    nrow = K3,
    ncol = K3,
    byrow = T
  ),
  parms.emission = list(lambda = c(10, 20, 30)),
  dens.emission = dpois.hsmm
)
# INICIALIZAMOS CON VALORES RANDOMS POR OBSERVACIONES
mod.hmm.k3 <- hmmfit(x[, 2], start_val3, mstep = mstep.pois)
plot(
  mod.hmm.k3$loglik,
  type = "b",
  ylab = "Log-likelihood",
  xlab = "Iteration",
  main = "Convergence K=3"
)
states <- mod.hmm.k3$yhat
plot(x[, 2], col = states, main = "# of earthquakes")
abline(h = mod.hmm.k3$model$parms.emission$lambda[1])
abline(h = mod.hmm.k3$model$parms.emission$lambda[2], col = 2)
abline(h = mod.hmm.k3$model$parms.emission$lambda[3], col = 3)

mod.hmm.k3 %>% summary
mod.hmm.k3$model

# 4 Clusters
K4 <- 4
start_val4 <- hmmspec(
  init = rep(1 / K4, K4),
  trans = matrix(
    1 / K4,
    nrow = K4,
    ncol = K4,
    byrow = T
  ),
  parms.emission = list(lambda = c(12, 15, 22, 30)),
  dens.emission = dpois.hsmm
)
# INICIALIZAMOS CON VALORES RANDOMS POR OBSERVACIONES
mod.hmm.k4 <- hmmfit(x[, 2], start_val4, mstep = mstep.pois)
plot(
  mod.hmm.k4$loglik,
  type = "b",
  ylab = "Log-likelihood",
  xlab = "Iteration",
  main = "Convergence K=4"
)
states <- mod.hmm.k4$yhat # IS JUST CLUSTERING
plot(x[, 2], col = states, main = "# of earthquakes")
abline(h = mod.hmm.k4$model$parms.emission$lambda[1])
abline(h = mod.hmm.k4$model$parms.emission$lambda[2], col = 2)
abline(h = mod.hmm.k4$model$parms.emission$lambda[3], col = 3)
abline(h = mod.hmm.k4$model$parms.emission$lambda[4], col = 4)

mod.hmm.k4 %>% summary
mod.hmm.k4$model
round(mod.hmm.k4$p, 3)

# Transient Cluster --> CUANDO TENGO 0 EN LA DIAGONAL, ESTOY OBLIGADO A IRME DE
# ESA POSICION EN EL T CUANDO ANDO EN ESA EN T-1
# SI HAY 0 ES PROBABLEMENTE QUE ESTEMOS EN OVERFITTING


# MODEL SELECTION
# Cuántos parámetros hay en un HMM


