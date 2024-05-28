# 3b-nngp-sim-example.R: this script fits a spatial occupancy model using
#                        a Nearest Neighbor Gaussian Process with a
#                        simulated data set for direct comparison with a
#                        spatial occupancy model that uses a full GP.
rm(list = ls())
library(spOccupancy)

# Simulate Data -----------------------------------------------------------
# Set the random seed. Do not change this, as otherwise you will fit the model
# for a different data set compared to the 3a-gp-sim-example.R script.
set.seed(400)
J.x <- 40
J.y <- 40
# Total number of spatial locations.
J <- J.x * J.y
# Number of replicate surveys.
n.rep <- rep(3, J)
# Occupancy intercept and regression coefficient for a single covariate.
beta <- c(0, 0.2)
# Detection intercept and regression coefficient for a single covariate.
alpha <- c(0.3, 0.5)
# Spatial decay parameter (the data are simulated across a 1 x 1 unit square
# and so the effective spatial range of 0.7 is quite large relative to the
# study area).
phi <- 3 / .7
# Spatial variance
sigma.sq <- 1.5
# These can be used to specify random effects, but here we don't do that.
p.RE <- list()
psi.RE <- list()
# Fit a spatial model
sp <- TRUE
# Use an exponential covariance function.
cov.model = 'exponential'

# Generate the data with simOcc
dat <- simOcc(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
              p.RE = p.RE, psi.RE = psi.RE, sigma.sq = sigma.sq, phi = phi, sp = TRUE,
              cov.model = cov.model)

# Get everything in spOccupancy format ------------------------------------
y <- dat$y
X.p <- dat$X.p
X <- dat$X
coords <- as.matrix(dat$coords)

# Occupancy covariates
occ.covs <- X
colnames(occ.covs) <- c('int', 'occ.cov.1')
# Detection covariates
det.covs <- list(int = X.p[, , 1],
                 det.cov.1 = X.p[, , 2])
# List in the format necessary for spOccupancy
data.list <- list(y = y, # Detection-nondetection data
                  occ.covs = occ.covs, # Occupancy covariates
                  det.covs = det.covs, # Detection covariates
                  coords = coords.fit) # Coordinates

# Priors ------------------------------------------------------------------
# Manually specifying priors for the spatial parameters, then setting all
# the other ones using the defaults.
prior.list <- list(sigma.sq.ig = c(2, 1),
                   phi.unif = c(3 / 1, 3 / .1))
# Initial values and tuning values ----------------------------------------
inits.list <- list(alpha = 0, beta = 0, phi = 3 / .5, sigma.sq = 1,
                   z = apply(y, 1, max, na.rm = TRUE))
# Tuning
tuning.list <- list(phi = 0.5)

# MCMC criteria -----------------------------------------------------------
batch.length <- 25
n.batch <- 2000
# Total number of MCMC samples
n.batch * batch.length
n.burn <- 30000
n.thin <- 5
n.chains <- 1

# Run the model -----------------------------------------------------------
# Approx run time: 7 min
out <- spPGOcc(occ.formula = ~ occ.cov.1,
               det.formula = ~ det.cov.1,
               data = data.list,
               n.batch = n.batch,
               batch.length = batch.length,
               inits = inits.list,
               priors = prior.list,
               accept.rate = 0.43,
               cov.model = "exponential",
               tuning = tuning.list,
               n.omp.threads = 1,
               verbose = TRUE,
               NNGP = FALSE,
               n.report = 10,
               n.burn = n.burn,
               n.thin = n.thin,
               n.chains = n.chains)

# Save spatial effects and simulated data to hard drive -------------------
w.means.nngp <- apply(out$w.samples, 2, mean) + mean(out$beta.samples[, 1])
save(dat, w.means.nngp, file = '~/Dropbox/talks/switzerland24/workshop/results/spPGOcc-NNGP-sims.rda')
