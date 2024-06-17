# 10b-sim-spatial-nmix.R: script to fit a single-species spatial N-mixture model
#                         to predict the abundance of a simulated species across a
#                         simulated landscape.
rm(list = ls())
library(spAbundance)
library(ggplot2)
library(sf)
library(stars)
library(patchwork)
# For a color palette that I like.
library(pals)
# Make sure to set seed to get same simulated data set.
set.seed(392)

# Simulate Data -----------------------------------------------------------
# Simulate abundance across a 50 x 50 grid (2500 locations)
J.x <- 50
J.y <- 50
J <- J.x * J.y
# A maximum of 5 replicates at each location
n.rep <- sample(5, J, replace = TRUE)
# Abundance is modeled according to an intercept and a single covariate
beta <- c(0, 0.5)
p.abund <- length(beta)
# Detection is modeled according to an intercept and two covariates
alpha <- c(0.5, -0.3, -0.5)
p.det <- length(alpha)
# Spatial decay parameter
phi <- 3 / .8
# Spatial variance parameter
sigma.sq <- 0.5
# Simulate data with a spatial effect
sp <- TRUE
# Use an exponential covariance model.
cov.model <- 'exponential'
# Simulate data from a Poisson distribution.
family <- 'Poisson'
# Run simNMix to simulate the data.
dat <- simNMix(J.x = J.x, J.y = J.y, n.rep = n.rep, beta = beta, alpha = alpha,
               sp = sp, phi = phi, sigma.sq = sigma.sq, cov.model = cov.model,
               family = family)
# Plot abundance across the landscape
plot.df <- data.frame(x = dat$coords[, 1],
                      y = dat$coords[, 2],
                      N = dat$N,
                      w = dat$w)
# Visualize the simulated abundance and spatial effect across the landscape
pred.stars <- st_as_stars(plot.df, dims = c('x', 'y'))
ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = N),interpolate = FALSE) +
  scale_fill_gradientn("", colors = ocean.tempo(1000),
                       na.value = NA) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "Easting", y = "Northing", title = 'Latent abundance')
ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = w),interpolate = FALSE) +
  scale_fill_gradientn("", colors = ocean.tempo(1000),
                       na.value = NA) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "Easting", y = "Northing", title = 'Spatial random effect')

# Randomly select 300 cells for model fitting and format for spAbundance
sampled.indx <- sample(1:J, 300, replace = FALSE)
y <- dat$y[sampled.indx, ]
abund.covs <- data.frame(habitat = dat$X[sampled.indx, 2])
det.covs <- list(det.cov.1 = dat$X.p[sampled.indx, , 2],
                 det.cov.2 = dat$X.p[sampled.indx, , 3])
coords <- dat$coords[sampled.indx, ]
data.list <- list(y = y, abund.covs = abund.covs, det.covs = det.covs,
                  coords = coords)
str(data.list)
# We will fit the model with these 300 locations, and then predict across the
# landscape to see how well our model approximates the true simulated patterns.

# Full abundance covariates and coordinates for prediction across the area
X.0 <- as.matrix(data.frame(intercept = 1,
                            habitat = dat$X[, 2]))
coords.0 <- dat$coords

# Fit a spatial N-mixture model -------------------------------------------
# Specify prior distributions
prior.list <- list(beta.normal = list(0, 10),
                   alpha.normal = list(0, 2.72),
                   # Data are simulated along a unit square study region
                   phi.unif = c(3 / 1, 3 / .1),
                   sigma.sq.ig = c(2, 1))
# Specify initial values
inits.list <- list(beta = 0, alpha = 0, phi = 3 / .5, sigma.sq = 1,
                   N = apply(data.list$y, 1, max, na.rm = TRUE))
# Specify initial tuning variances
tuning <- list(beta = 0.5, alpha = 0.1, phi = 1, w = 1)

# Specify MCMC criteria
n.batch <- 1600
batch.length <- 25
n.burn <- 20000
n.thin <- 10
n.chains <- 3

# Approx run time: 3 min
out <- spNMix(abund.formula = ~ habitat,
              det.formula = ~ det.cov.1 + det.cov.2,
              data = data.list, priors = prior.list, inits = inits.list,
              tuning = tuning, n.batch = n.batch, batch.length = batch.length,
              n.burn = n.burn, n.thin = n.thin, n.chains = n.chains, n.report = 400,
              NNGP = TRUE, n.neighbors = 8, cov.model = 'exponential')

summary(out)

# Compare model estimates to simulated data -------------------------------
# Compare estimated expected abundance (mu) with the true expected abundance
# Posterior means of the expected abundance values
mu.means <- apply(out$mu.samples, 2, mean)
# 95% credible interval
mu.quants <- apply(out$mu.samples, 2, quantile, c(0.025, 0.975))
mu.true <- dat$mu[sampled.indx]
plot(mu.true, mu.means, pch = 19, xlab = 'True expected abundance',
     ylab = 'Estimated expected abundance', ylim = range(mu.quants))
segments(mu.true, mu.quants[1, ], mu.true, mu.quants[2, ])
abline(0, 1)

# Compare the estimated spatial random effect (w) with the true effect
# Posterior means of the spatial random effect
w.means <- apply(out$w.samples, 2, mean)
# 95% credible interval
w.quants <- apply(out$w.samples, 2, quantile, c(0.025, 0.975))
w.true <- dat$w[sampled.indx]
plot(w.true, w.means, pch = 19, xlab = 'True spatial effect (w)',
     ylab = 'Estimated spatial effect (w)', ylim = range(w.quants))
segments(w.true, w.quants[1, ], w.true, w.quants[2, ])
abline(0, 1)

# Predict across the entire square region ---------------------------------
# The prediction design matrix
str(X.0)
# The prediction coordinates
str(coords.0)
# Just for fun let's compare run time when using 4 threads vs. 1 thread. Note
# that this may not work on all of your machines.
# Four threads
curr.time <- Sys.time()
out.pred <- predict(out, X.0, coords.0, n.report = 100, n.omp.threads = 4)
run.time.4 <- Sys.time() - curr.time
# One thread
curr.time <- Sys.time()
out.pred.1 <- predict(out, X.0, coords.0, n.report = 100, n.omp.threads = 1)
run.time.1 <- Sys.time() - curr.time
# Slight improvements in run time, for even just this small example.
run.time.4
run.time.1

# Generate map of predicted abundance and the predicted spatial effect, and
# compare with the truth
# Expected abundance means
mu.0.means <- apply(out.pred$mu.0.samples, 2, mean)
w.0.means <- apply(out.pred$w.0.samples, 2, mean)

plot.pred.df <- data.frame(x = coords.0[, 1], y = coords.0[, 2],
                           mu.0 = mu.0.means, w.0 = w.0.means,
                           mu.true = dat$mu, w.true = dat$w)
pred.stars <- st_as_stars(plot.pred.df, dims = c('x', 'y'))
mu.plot.est <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = mu.0),interpolate = FALSE) +
  scale_fill_gradientn("", colors = ocean.tempo(1000),
                       na.value = NA) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "Easting", y = "Northing", title = '(a) Predicted expected abundance')
mu.plot.true <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = mu.true),interpolate = FALSE) +
  scale_fill_gradientn("", colors = ocean.tempo(1000),
                       na.value = NA) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "Easting", y = "Northing", title = '(b) True expected abundance')
w.plot.est <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = w.0),interpolate = FALSE) +
  scale_fill_gradientn("", colors = ocean.tempo(1000),
                       na.value = NA) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "Easting", y = "Northing", title = '(c) Predicted spatial effect')
w.plot.true <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = w.true),interpolate = FALSE) +
  scale_fill_gradientn("", colors = ocean.tempo(1000),
                       na.value = NA) +
  theme_bw() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "Easting", y = "Northing", title = '(d) True spatial effect')
# Looks pretty good!
mu.plot.est + mu.plot.true + w.plot.est + w.plot.true
