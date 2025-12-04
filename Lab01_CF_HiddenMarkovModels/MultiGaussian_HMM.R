library(mhsmm)
require(tidyverse)
require(gamlss)

dax4 <-
  read.csv(
    url("http://www.hmms-for-time-series.de/second/data/dax4.txt"),
    sep = "",
    col.names = c("Allianz", "Daimlerchrisler", "Deustche Bank", "Siemens")
  )

K <- 2
start.val.m <- hmmspec(
  init = c(1, 0),
  trans = matrix(
    c(.8, .2, .2, .8),
    nrow = K,
    ncol = K,
    byrow = T
  ),
  parms.emission = list(mu = list(c(0, 0.1, 0.3, -0.1), # deben haber promedios de cada variable, y c() por cada cluster
                                  c(.1, .2, .3, .1)),
                        sigma = list(diag(4), diag(4))),
  dens.emission = dmvnorm.hsmm
)

mod.k2.m <- hmmfit(matrix(unlist(dax4[, 1:4]), ncol = 4), start.val.m,
                   mstep = mstep.mvnorm)

# mstep.norm depende del x y de los weighted probabilites posterior probability
# wt is a matrix de las probabilidades para pasar w_tk = Pz(Ct = k|X)

plot(
  mod.k2.m$loglik,
  type = "b",
  ylab = "Log-likelihood",
  xlab = "Iteration",
  main = "Convergence K=2 / Multivariate"
)

clusters.m <- mod.k2.m$yhat #Secuencia de 1 y 2, como se clusterizo
plot(dax4, col = clusters.m)
# captures volatility
clusters.m %>% table

mod.k2$model