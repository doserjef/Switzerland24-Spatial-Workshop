# 11c-european-goldfinch-spatial-nmix.R: in the previous exercise, we used simulated
#                                        data to show how to fit spatial N-mixture
#                                        models, and saw that everything went very
#                                        smoothly. Here we will work with a real
#                                        data set (the same European Goldfinch
#                                        data we worked with in our occupancy modeling
#                                        section but now using the raw counts). We will
#                                        see how spatial (and non-spatial) N-mixture models
#                                        can be difficult to estimate in practice and can
#                                        show very bad convergence properties.
# Data source citations:
#   Kéry, M. & Royle, J.A. (2016) _Applied Hierarchical Modeling in Ecology_ AHM1 - 11.3.
#   Swiss Federal Statistical Office (http://www.bfs.admin.ch)
#   Data were derived from objects included in the AHMBook and unmarked R packages.
rm(list = ls())
library(spAbundance)
library(ggplot2)
library(sf)
# For some color palettes I like.
library(pals)
set.seed(32020)

# Load the data set and explore the data ----------------------------------
# Note the data are currently formatted for multi-species N-mixture models.
# Loads an object called data.swiss.mhb
load('data/swiss-mhb-count-data.rda')
str(data.swiss.mhb)
# Plot the locations to refamiliarize yourself with the MHB data. Remember that
# we've seen these data a couple of times before.
plot(data.swiss.mhb$coords, pch = 19)
# Subset the data for just working with European goldfinch.
sp.names <- dimnames(data.swiss.mhb$y)[[1]]
indx <- which(sp.names == 'European Goldfinch')
data.goldfinch <- data.swiss.mhb
data.goldfinch$y <- data.swiss.mhb$y[indx, , ]
str(data.goldfinch)

# Plot the maximum observed counts for the European goldfinch at each site
coords.sf <- st_as_sf(as.data.frame(data.goldfinch$coords),
                      coords = c('x', 'y'))
coords.sf$max.count <- apply(data.goldfinch$y, 1, max, na.rm = TRUE)
ggplot() +
  geom_sf(data = coords.sf, aes(fill = max.count), pch = 21, size = 2) +
  scale_fill_gradientn("", colors = ocean.tempo(1000), na.value = NA) +
  theme_bw()
table(data.goldfinch$y)
# Note the seemingly massive overdispersion in the observed count data, where most
# counts are 0, but every once in a while we have very large count values. Such
# situations can make fitting spatial abundance models pretty tricky. Spatial
# N-mixture models in particular are difficult to fit under such situations as
# there are specific statistical identifiability limitations.

# Poisson non-spatial -----------------------------------------------------
# Let's start off with a modest number of MCMC samples. Again, remember from lecture
# that models in spAbundance use a less efficient algorithm than those in spOccupancy,
# and so you will generally find yourself needing to run models for more MCMC samples
# in spAbundance compared to spOccupancy.
n.batch <- 1800
batch.length <- 25
n.batch * batch.length
n.burn <- 20000
n.thin <- 25
n.chains <- 3
# Fit the model -----------------------
# Approx. run time: ~1.5 min
out.non.sp <- NMix(abund.formula = ~ scale(elevation) + I(scale(elevation)^2) +
                                     scale(forest),
                   det.formula = ~ scale(date) + I(scale(date^2)) + scale(dur) +
                                   (1 | obs),
                   data = data.goldfinch, family = 'Poisson',
                   n.omp.threads = 1, n.report = 200, n.batch = n.batch,
                   batch.length = batch.length, n.burn = n.burn, n.thin = n.thin,
                   n.chains = n.chains)
# Take a look at the model summary and some traceplots to assess convergence
summary(out.non.sp)
# Abundance coefficients
plot(out.non.sp, 'beta', density = FALSE)
# Detection coefficients
plot(out.non.sp, 'alpha', density = FALSE)
# Haven't quite reached convergence, so we would in reality like to run this model
# longer, but let's proceed for now.
# Assess model fit --------------------
# Let's generate a plot of the fitted values vs. the actual data values
out.non.sp.fitted <- fitted(out.non.sp, type = 'marginal')
str(out.non.sp.fitted)
# Mean fitted values
y.rep.means <- apply(out.non.sp.fitted$y.rep.samples, c(2, 3), mean)
plot(data.goldfinch$y, y.rep.means, pch = 19,
     xlim = range(data.goldfinch$y, na.rm = TRUE),
     ylim = range(data.goldfinch$y, na.rm = TRUE))
abline(0, 1)
# We see pretty decent correspondence between the true data values and the
# model generated data values.
# Plot estimated abundance at the site locations
coords.sf$mu.non.sp.means <- apply(out.non.sp$mu.samples, 2, mean)
ggplot() +
  geom_sf(data = coords.sf, aes(fill = mu.non.sp.means), pch = 21, size = 2) +
  scale_fill_gradientn("", colors = ocean.tempo(1000), na.value = NA) +
  theme_bw()
table(data.goldfinch$y)

# Fit a spatial N-mixture model -------------------------------------------
# Let's try and fit a spatial N-mixture model with the default priors
out.sp <- spNMix(abund.formula = ~ scale(elevation) + I(scale(elevation)^2) +
                                   scale(forest),
                 det.formula = ~ scale(date) + I(scale(date^2)) + scale(dur) +
                                 (1 | obs),
                 data = data.goldfinch, family = 'Poisson',
                 NNGP = TRUE, n.neighbors = 7, cov.model = 'exponential',
                 n.omp.threads = 1, n.report = 200,
                 n.batch = n.batch, batch.length = batch.length,
                 n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)

# Let's see if the model converged...
summary(out.sp)
# Definitely have not reached convergence. Let's look at some trace plots to dig
# more into this.
# Abundance coefficients
plot(out.sp, 'beta', density = FALSE)
# Detection coefficients
plot(out.sp, 'alpha', density = FALSE)

# Look at the traceplots for both the abundance intercept and detection intercept.
# Notice that they show seemingly opposing patterns, indicating that they may
# not be identifiable from each other. This is the concept we discussed in lecture
# related to it being very difficult to separately distinguish between overdispersion
# in latent abundance versus overdispersion in detection probability.
# We can further look at this by looking at the correlation between the posterior
# samples of the abundance intercept and the detection intercept
cor(out.sp$beta.samples[, 1], out.sp$alpha.samples[, 1])
plot(c(out.sp$beta.samples[, 1]), c(out.sp$alpha.samples[, 1]), pch = 19,
     xlab = 'Abundance Intercept MCMC Samples', ylab = 'Detection Intercept MCMC Samples')
# Such a pattern is a clear indication that there is some sort of "identifiability"
# problem in the model, and that we may need to place more strict requirements on
# the model parameters.

# Fit a spatial N-mixture model with more restrictions --------------------
# Let's manually specify the prior distribution for the spatial decay parameter
# based on the inntersite distance matrix
dist.mat <- dist(data.goldfinch$coords)
quantile(dist.mat)
quantile(dist.mat, .1)
# Let's set the lower bound on the effective spatial range to 40,000m.
priors <- list(phi.unif = c(3 / max(dist.mat), 3 / 40000))
# We will also set more informative priors for the abundance regression
# coefficients by setting the variance to 1 (the default is 100). By setting the variance
# to 1 instead of the default of 100, we are giving more weight to values of the abundance
# intercept close to the prior mean (in our case 0) as opposed to values farther away
# from the mean. This may potentially help in restricting the possibility of extremely
# large abundance values.
priors$beta.normal <- list(mean = 0, var = 1)
# Lastly, let's set the inverse-gamma prior for sigma.sq. In spatial N-mixture models,
# unreasonably large values of the spatial random effects can be generated if the
# prior on the spatial variance is too uninformative. Here we set the shape parameter
# to 2 and the scale parameter to 0.1, which implies a prior mean of sigma.sq = 0.1
# with an infinite variance. Thus, this is still a fairly weak prior, but it may help to
# "protect" us from extremely large values of the spatial variance unless there is a lot
# of support for it in the data.
priors$sigma.sq.ig <- c(2, 0.1)
out.sp.2 <- spNMix(abund.formula = ~ scale(elevation) + I(scale(elevation)^2) +
                                     scale(forest),
                   det.formula = ~ scale(date) + I(scale(date^2)) + scale(dur) +
                                   (1 | obs),
                   data = data.goldfinch, family = 'Poisson',
                   NNGP = TRUE, n.neighbors = 7, cov.model = 'exponential',
                   priors = priors, n.omp.threads = 1, n.report = 200,
                   n.batch = n.batch, batch.length = batch.length,
                   n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)

# Let's see if the new priors helped us achieve better convergence of the model.
summary(out.sp.2)
# Okay, looks slightly better, but still not fully converged. Let's look at some
# traceplots
# Abundance coefficients
plot(out.sp.2, 'beta', density = FALSE)
# Detection coefficients
plot(out.sp.2, 'alpha', density = FALSE)
# Are the detection intercept and abundance intercept still highly correlated?
cor(out.sp.2$beta.samples[, 1], out.sp.2$alpha.samples[, 1])
# Still correlated, but not nearly as much as the model using the default priors.

