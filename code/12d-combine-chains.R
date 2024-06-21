# 11d-combine-chains.R: this script combines the three chains that are run in
#                       parallel into a single spAbundance model object, such that
#                       the object can subsequently be used with all spAbundance
#                       model functions.
rm(list = ls())
library(spAbundance)
library(coda)
library(abind)

# Read in the model objects from the three chains -------------------------
load('results/hbef-sfMsNMix-results-chain-1.rda')
out.1 <- out
load('results/hbef-sfMsNMix-results-chain-2.rda')
out.2 <- out
load('results/hbef-sfMsNMix-results-chain-3.rda')
out.3 <- out

out.full <- list()
# Combine all MCMC samples together ---------------------------------------
out.full$n.chains <- 3
# Community-level abundance effects
out.full$beta.comm.samples <- mcmc(rbind(out.1$beta.comm.samples,
                                         out.2$beta.comm.samples,
                                         out.3$beta.comm.samples))
# Community-level abundance variances
out.full$tau.sq.beta.samples <- mcmc(rbind(out.1$tau.sq.beta.samples,
                                           out.2$tau.sq.beta.samples,
                                           out.3$tau.sq.beta.samples))
# Community-level detection effects
out.full$alpha.comm.samples <- mcmc(rbind(out.1$alpha.comm.samples,
                                         out.2$alpha.comm.samples,
                                         out.3$alpha.comm.samples))
# Community-level detection variances
out.full$tau.sq.alpha.samples <- mcmc(rbind(out.1$tau.sq.alpha.samples,
                                           out.2$tau.sq.alpha.samples,
                                           out.3$tau.sq.alpha.samples))
# Species-level abundance regression coefficients
out.full$beta.samples <- mcmc(rbind(out.1$beta.samples,
                                    out.2$beta.samples,
                                    out.3$beta.samples))
# Species-level detection regression coefficients
out.full$alpha.samples <- mcmc(rbind(out.1$alpha.samples,
                                    out.2$alpha.samples,
                                    out.3$alpha.samples))
# Fitted values
out.full$N.samples <- abind(out.1$N.samples, out.2$N.samples, out.3$N.samples, along = 1)
# Expected abundance values
out.full$mu.samples <- abind(out.1$mu.samples, out.2$mu.samples, out.3$mu.samples, along = 1)
# Spatial decay parameters
out.full$theta.samples <- mcmc(rbind(out.1$theta.samples,
                                    out.2$theta.samples,
                                    out.3$theta.samples))
# Spatial factor loadings
out.full$lambda.samples <- mcmc(rbind(out.1$lambda.samples,
                                    out.2$lambda.samples,
                                    out.3$lambda.samples))
# Spatial factors
out.full$w.samples <- abind(out.1$w.samples, out.2$w.samples, out.3$w.samples, along = 1)
# Get RHat values for the main parameters ---------------------------------
out.full$rhat <- out.1$rhat
# beta.comm
tmp <- mcmc.list(out.1$beta.comm.samples, out.2$beta.comm.samples, out.3$beta.comm.samples)
out.full$rhat$beta.comm <- as.vector(gelman.diag(tmp, autoburnin = FALSE)$psrf[, 2])
# tau.sq.beta
tmp <- mcmc.list(out.1$tau.sq.beta.samples, out.2$tau.sq.beta.samples, out.3$tau.sq.beta.samples)
out.full$rhat$tau.sq.beta <- as.vector(gelman.diag(tmp, autoburnin = FALSE)$psrf[, 2])
# beta
tmp <- mcmc.list(out.1$beta.samples, out.2$beta.samples, out.3$beta.samples)
out.full$rhat$beta <- as.vector(gelman.diag(tmp, autoburnin = FALSE)$psrf[, 2])
# alpha.comm
tmp <- mcmc.list(out.1$alpha.comm.samples, out.2$alpha.comm.samples, out.3$alpha.comm.samples)
out.full$rhat$alpha.comm <- as.vector(gelman.diag(tmp, autoburnin = FALSE)$psrf[, 2])
# tau.sq.alpha
tmp <- mcmc.list(out.1$tau.sq.alpha.samples, out.2$tau.sq.alpha.samples, out.3$tau.sq.alpha.samples)
out.full$rhat$tau.sq.alpha <- as.vector(gelman.diag(tmp, autoburnin = FALSE)$psrf[, 2])
# alpha
tmp <- mcmc.list(out.1$alpha.samples, out.2$alpha.samples, out.3$alpha.samples)
out.full$rhat$alpha <- as.vector(gelman.diag(tmp, autoburnin = FALSE)$psrf[, 2])
# theta
tmp <- mcmc.list(out.1$theta.samples, out.2$theta.samples, out.3$theta.samples)
out.full$rhat$theta <- as.vector(gelman.diag(tmp, autoburnin = FALSE)$psrf[, 2])
# lambda
tmp <- mcmc.list(out.1$lambda.samples, out.2$lambda.samples, out.3$lambda.samples)
out.full$rhat$lambda <- as.vector(gelman.diag(tmp, autoburnin = FALSE, multivariate = FALSE)$psrf[, 2])

# Get ESS values for the main parameters ----------------------------------
out.full$ESS <- list()
out.full$ESS$beta.comm <- effectiveSize(out.full$beta.comm.samples)
out.full$ESS$tau.sq.beta <- effectiveSize(out.full$tau.sq.beta.samples)
out.full$ESS$beta <- effectiveSize(out.full$beta.samples)
out.full$ESS$alpha.comm <- effectiveSize(out.full$alpha.comm.samples)
out.full$ESS$tau.sq.alpha <- effectiveSize(out.full$tau.sq.alpha.samples)
out.full$ESS$alpha <- effectiveSize(out.full$alpha.samples)
out.full$ESS$lambda <- effectiveSize(out.full$lambda.samples)
out.full$ESS$theta <- effectiveSize(out.full$theta.samples)

# Other stuff -------------------------------------------------------------
# This stuff is used under the hood for use with summaries, plotting,
# prediction, WAIC,  PPCs.
# Design matrix for fixed abundance effects
out.full$X <- out.1$X
# Design matrix for fixed design effects
out.full$X.p <- out.1$X.p
# Design matrix for abundance random effects
out.full$X.re <- out.1$X.re
out.full$X.random <- out.1$X.random
# Design matrix for detection random effects
out.full$X.p.re <- out.1$X.p.re
out.full$X.p.random <- out.1$X.p.random
# The count data
out.full$y <- out.1$y
# Offset if there is one
out.full$offset <- out.1$offset
# Information on the call to the function
out.full$call <- out.1$call
# Overall number of samples run per chain
out.full$n.samples <- out.1$n.samples
# Names of columns in X
out.full$x.names <- out.1$x.names
# Names of columns in X.p
out.full$x.p.names <- out.1$x.p.names
# Species names
out.full$sp.names <- out.1$sp.names
# Number of posterior samples saved for each chain
out.full$n.post <- out.1$n.post
# Thinning rate
out.full$n.thin <- out.1$n.thin
# Amount of burn-in
out.full$n.burn <- out.1$n.burn
# Number of chains
out.full$n.chains <- 3
# Distribution used for abundance (NB or Poisson)
out.full$dist <- out.1$dist
# Names of columns with random effects
out.full$re.cols <- out.1$re.cols
out.full$re.det.cols <- out.1$re.det.cols
# Names of spatial factors
out.full$theta.names <- out.1$theta.names
# NNGP or GP
out.full$type <- out.1$type
# Index for covariance function
out.full$cov.model.indx <- out.1$cov.model.indx
# Coordinates
out.full$coords <- out.1$coords
# Number of neighbors
out.full$n.neighbors <- out.1$n.neighbors
# Number of factors
out.full$q <- out.1$q
# Logical indicating if there were any random effects
out.full$muRE <- out.1$muRE
out.full$pRE <- out.1$pRE
# Run time
# Setting the "overall" run time to be the longest run time across the three chains
tmp <- which.max(c(out.1$run.time[3], out.2$run.time[3], out.3$run.time[3]))
if (tmp == 1) out.full$run.time <- out.1$run.time
if (tmp == 2) out.full$run.time <- out.2$run.time
if (tmp == 3) out.full$run.time <- out.3$run.time

# Make sure the new object has class sfMsNMix ------------------------------
class(out.full) <- 'sfMsNMix'

# Check to make sure it's working -----------------------------------------
# These should all just return without an error.
summary(out.full)
plot(out.full, 'beta.comm', density = FALSE)

# Save full object --------------------------------------------------------
save(out.full, file = 'results/hbef-sfMsNmix-three-chains.rda')
