library(mhsmm)
plot(earthquake[, 2], type = "l")
# Hidden Markov models
K <- 2 #number of clusters
start.val <- hmmspec(
  init = rep(1 / K, K),
  trans = matrix(1 / K, nrow = K, ncol = K),
  parms.emission  = list(lambda = c(15, 25)),
  dens.emission  = dpois.hsmm
)
mod.hmm.k2.bis <-
  hmmfit(earthquake[, 2], start.val, mstep = mstep.pois)
plot(
  mod.hmm.k2.bis$loglik,
  type = "b",
  ylab = "Log-likelihood",
  xlab = "Iteration"
)
states <- mod.hmm.k2$yhat
plot(earthquake[, 2], col = states, main = "# of earthquakes")
abline(h = mod.hmm.k2$model$parms.emission$lambda[1])
abline(h = mod.hmm.k2$model$parms.emission$lambda[2], col = 2)
K <- 3
start.val <- hmmspec(
  init = rep(1 / K, K),
  trans = matrix(1 / K, nrow = K, ncol = K),
  parms.emis = list(lambda = c(10, 20, 30)),
  dens.emis = dpois.hsmm
)
mod.hmm.k3 <- hmmfit(earthquake[, 2], start.val, mstep = mstep.pois)
plot(mod.hmm.k3$loglik,
     type = "b",
     ylab = "Log-likelihood",
     xlab = "Iteration")
states <- mod.hmm.k3$yhat
plot(earthquake[, 2], col = states, main = "# of earthquakes")
abline(h = mod.hmm.k3$model$parms.emission$lambda[1])
abline(h = mod.hmm.k3$model$parms.emission$lambda[2], col = 2)
abline(h = mod.hmm.k3$model$parms.emission$lambda[3], col = 3)
K <- 4
start.val <- hmmspec(
  init = rep(1 / K, K),
  trans = matrix(1 / K, nrow = K, ncol = K),
  parms.emis = list(lambda = c(10, 20, 30, 40)),
  dens.emis = dpois.hsmm
)
mod.hmm.k4 <- hmmfit(earthquake[, 2], start.val, mstep = mstep.pois)
states <- mod.hmm.k4$yhat
plot(earthquake[, 2], col = states, main = "# of earthquakes")
abline(h = mod.hmm.k4$model$parms.emission$lambda[1])
abline(h = mod.hmm.k4$model$parms.emission$lambda[2], col = 2)
abline(h = mod.hmm.k4$model$parms.emission$lambda[3], col = 3)
abline(h = mod.hmm.k4$model$parms.emission$lambda[4], col = 4)
#
#
# Bootstrap/Simulation: an example
#
B <- 1000 # replicates
lambda.boot <- matrix(NA, 2, B)
for (b in 1:B)
{
  true.par <- hmmspec(
    init = mod.hmm.k2$model$init,
    trans = mod.hmm.k2$model$transition,
    parms.emis = list(lambda = mod.hmm.k2$model$parms.emission$lambda),
    dens.emis = dpois.hsmm
  )
  train <-
    simulate(true.par,
             nsim = 107,
             seed = b,
             rand.emis = rpois.hsmm)
  mod.boot <- hmmfit(train, true.par, mstep = mstep.pois)
  lambda.boot[, b] <- mod.boot$model$parms.emission$lambda
}
apply(lambda.boot, 1, mean)
apply(lambda.boot, 1, sd)
#

dmvnorm.hsmm <- function(x, j, model)
{
  dmvnorm(x,
          mean = model$parms.emission$mu[[j]],
          sigma = model$parms.emission$sigma[[j]])
}

mstep.mvnorm <- function(x, wt) {
  emission <- list(mu = list(), sigma = list())
  for (i in 1:ncol(wt)) {
    tmp <- cov.wt(x, wt[, i])
    emission$mu[[i]] <- tmp$center
    emission$sigma[[i]] <- tmp$cov
  }
  emission
}

K <- 2
start.val <- hmmspec(
  init = c(1, 0),
  trans = matrix(
    c(.9, .1, .1, .9),
    byrow = T,
    nrow = K,
    ncol = K
  ),
  parms.emis = list(mu = list(c(0, 0.1, 0.3, -0.1), c(.1, .2, .3, .1)),
                    sigma = list(diag(4), diag(4))),
  dens.emis = dmvnorm.hsmm
)
mod.hmm.k2 <- hmmfit(matrix(unlist(dax4[, 1:4]), ncol = 4),
                     start.val, mstep = mstep.mvnorm)
#
# Faithful data
data("faithful")
K <- 2
start.val <- hmmspec(
  init = rep(1 / K, K),
  trans = matrix(1 / K, nrow = K, ncol = K),
  parms.emis = list(mu = c(50, 80), sigma = c(1, 1)),
  dens.emis = dnorm.hsmm
)
mod.hmm.k2 <- hmmfit(faithful[, 2], start.val, mstep = mstep.norm)
plot(mod.hmm.k2$loglik,
     type = "b",
     ylab = "Log-likelihood",
     xlab = "Iteration")
states <- mod.hmm.k2$yhat
plot(faithful[, 2], col = states, main = "Waiting times")
abline(h = mod.hmm.k2$model$parms.emission$mu[1])
abline(h = mod.hmm.k2$model$parms.emission$mu[2], col = 2)


K <- 3
start.val <- hmmspec(
  init = rep(1 / K, K),
  trans = matrix(
    c(.8, .1, .1, .1, .8, .1, .1, .1, .8),
    byrow = T,
    nrow = K,
    ncol = K
  ),
  parms.emis = list(mu = c(50, 80, 60), sigma = c(1, 1, 1)),
  dens.emis = dnorm.hsmm
)
mod.hmm.k3 <- hmmfit(faithful[, 2], start.val, mstep = mstep.norm)
plot(mod.hmm.k3$loglik,
     type = "b",
     ylab = "Log-likelihood",
     xlab = "Iteration")
states <- mod.hmm.k3$yhat
plot(faithful[, 2], col = states, main = "Waiting times")
abline(h = mod.hmm.k3$model$parms.emission$mu[1])
abline(h = mod.hmm.k3$model$parms.emission$mu[2], col = 2)
abline(h = mod.hmm.k3$model$parms.emission$mu[3], col = 3)

K <- 2
start.val <- hmmspec(
  init = c(1, 0),
  trans = matrix(
    c(.9, .1, .1, .9),
    byrow = T,
    nrow = K,
    ncol = K
  ),
  parms.emis = list(mu = list(c(2, 80), c(4, 50)),
                    sigma = list(diag(2), diag(2))),
  dens.emis = dmvnorm.hsmm
)
mod.hmm.k2 <- hmmfit(matrix(unlist(faithful), ncol = 2),
                     start.val, mstep = mstep.mvnorm)
plot(faithful, col = mod.hmm.k2$yhat)

initialization <- kmeans(faithful, 2)
kmeans.means <- initialization$centers
sigma1 <- cov(faithful[initialization$cluster == 1,])
sigma2 <- cov(faithful[initialization$cluster == 2,])

K <- 2
start.val <- hmmspec(
  init = c(1, 0),
  trans = matrix(
    c(.9, .1, .1, .9),
    byrow = T,
    nrow = K,
    ncol = K
  ),
  parms.emis = list(
    mu = list(kmeans.means[1,], kmeans.means[2,]),
    sigma = list(sigma1, sigma2)
  ),
  dens.emis = dmvnorm.hsmm
)
mod.hmm.k2.v1 <- hmmfit(matrix(unlist(faithful), ncol = 2),
                        start.val, mstep = mstep.mvnorm)
plot(faithful, col = mod.hmm.k2.v1$yhat)


initialization <- kmeans(dax4[, 1:4], 2)
kmeans.means <- initialization$centers
sigma1 <- cov(dax4[initialization$cluster == 1, 1:4])
sigma2 <- cov(dax4[initialization$cluster == 2, 1:4])

K <- 2
start.val <- hmmspec(
  init = c(1, 0),
  trans = matrix(
    c(.9, .1, .1, .9),
    byrow = T,
    nrow = K,
    ncol = K
  ),
  parms.emis = list(
    mu = list(kmeans.means[1,], kmeans.means[2,]),
    sigma = list(sigma1, sigma2)
  ),
  dens.emis = dmvnorm.hsmm
)
mod.hmm.k2 <- hmmfit(matrix(unlist(dax4[, 1:4]), ncol = 4),
                     start.val, mstep = mstep.mvnorm)
plot(dax4[, 1:4], col = mod.hmm.k2$yhat)

K <- 2
start.val <- hsmmspec(
  init = c(1, 0),
  trans = matrix(
    c(0, 1, 1, 0),
    byrow = T,
    nrow = K,
    ncol = K
  ),
  parms.emis = list(
    mu = list(kmeans.means[1,], kmeans.means[2,]),
    sigma = list(sigma1, sigma2)
  ),
  sojourn = list(
    lambda = c(10, 20),
    shift = c(1, 1),
    type = "poisson"
  ),
  dens.emis = dmvnorm.hsmm
)
mod.hsmm.k2 <- hsmmfit(matrix(unlist(dax4[, 1:4]), ncol = 4),
                       start.val, mstep = mstep.mvnorm)
plot(dax4[, 1:4], col = mod.hmm.k2$yhat)
