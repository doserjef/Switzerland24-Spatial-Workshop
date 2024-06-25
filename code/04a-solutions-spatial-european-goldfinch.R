# 04a-solutions-spatial-european-goldfinch-ex.R: in this exercise, you will build on the script
#                                                1-swiss-european-goldfinch.R by extending
#                                                our occupancy model to a spatial
#                                                occupancy model that you will use to
#                                                predict the distribution of the European
#                                                goldfinch across Switzerland.
# Data source citations:
# Kéry, M. & Royle, J.A. (2016) _Applied Hierarchical Modeling in Ecology_ AHM1 - 11.3.
# Swiss Federal Statistical Office (http://www.bfs.admin.ch)
# Data were derived from objects included in the AHMBook and unmarked R packages.
rm(list = ls())
library(spOccupancy)
# For summarizing MCMC results
library(MCMCvis)
# For making species distribution maps
library(ggplot2)
library(stars)
library(pals)
library(patchwork)
# If not using the RStudio project, set working directory to the repository
# directory.
# setwd("../")
# Set seed for same results
set.seed(3802)

# In this example, your goal is to produce a species distribution map of the
# European goldfinch throughout Switzerland. Based on what we found in
# 1-swiss-european-goldfinch.R, the top performing non-spatial occupancy model
# included a linear effect of forest cover and a linear and quadratic effect of
# elevation. Here you will compare this model to a spatial model and interpret the
# results from the spatial model.

# NOTE: we label the tasks in exercises throughout the course that you need to
#       complete with the label TODO, followed by the associated number.

# Data prep ---------------------------------------------------------------
# Read in the data source (reads in an object called data.goldfinch)
load("data/europeanGoldfinchSwiss.rda")
str(data.goldfinch)
# Take a quick look at the spatial locations
plot(data.goldfinch$coords, pch = 19)

# TODO 4.1 ----------------------------------------------------------------
# Fit a non-spatial occupancy model that includes a linear effect
# of forest cover and linear and quadratic effects of elevation. # Recall that
# we fit this exact model in "1-swiss-european-goldfinch.R". Make sure that the
# model converged. How do you know? Save the resulting model results in an
# object called "out.non.sp". Use default initial values and prior distributions.
out.non.sp <- PGOcc(occ.formula = ~ scale(elevation) + I(scale(elevation)^2) + scale(forest),
                    det.formula = ~ scale(date) + I(scale(date)^2) + scale(dur),
                    data = data.goldfinch,
                    n.samples = 5000,
                    n.thin = 4,
                    n.burn = 3000,
                    n.chains = 3,
                    n.report = 500)
summary(out.non.sp)
plot(out.non.sp, 'beta', density = FALSE)

# TODO 4.2 ----------------------------------------------------------------
# Fit the same model as in TODO 4.1, but now include a spatial random effect
# (i.e., fit a spatial occupancy model with the same covariates). Use an
# NNGP with 10 neighbors and an exponential spatial covariance function.
# Set the initial value of the spatial variance
# (sigma.sq) to 3 and the initial value of the spatial decay parameter to
# 3 divided by the mean distance between sites in the data set (you can calculate
# the distance matrix of all sites using the dist() function). Run three chains
# of the model for 5000 samples each , split into batches each with 25
# MCMC samples. Set the initial tuning variance for phi to 0.5.
# Determine if the model has converged. Throw away the first 3000 samples as burn-in,
# and use a thinning rate of 4 to only keep 500 samples per chain. Call the model
# output out.sp.1. Has the model converged with these MCMC settings?
dist.mat <- dist(data.goldfinch$coords)
inits <- list(sigma.sq = 3, phi = 3 / mean(dist.mat))
tuning = list(phi = 0.5)
out.sp.1 <- spPGOcc(occ.formula = ~ scale(elevation) + I(scale(elevation)^2) +
                                    scale(forest),
                    det.formula = ~ scale(date) + I(scale(date)^2) + scale(dur),
                    data = data.goldfinch, inits = inits, tuning = tuning,
                    n.batch = 200, batch.length = 25,
                    NNGP = TRUE, n.neighbors = 10, n.burn = 3000, n.thin = 4,
                    n.chains = 3, n.report = 50)
summary(out.sp.1)
plot(out.sp.1, 'beta', density = FALSE)
plot(out.sp.1, 'theta', density = FALSE)
# The model has not converged.

# TODO 4.3 ----------------------------------------------------------------
# Refit the same model as before, but now manually specify the prior on phi to
# have an effective spatial range that falls between 20,000 and the maximum
# distance between two sites in the data set. Recall that the effective spatial
# range when using the exponential covariance function is 3 / phi. Save the
# model as out.sp.2.
priors <- list(phi.unif = c(3 / max(dist.mat), 3 / 20000))
out.sp.2 <- spPGOcc(occ.formula = ~ scale(elevation) + I(scale(elevation)^2) +
                                    scale(forest),
                    det.formula = ~ scale(date) + I(scale(date)^2) + scale(dur),
                    data = data.goldfinch, inits = inits, priors = priors,
                    tuning = tuning, n.batch = 200, batch.length = 25,
                    NNGP = TRUE, n.neighbors = 10, n.burn = 3000, n.thin = 4,
                    n.chains = 3, n.report = 50)
summary(out.sp.2)
plot(out.sp.2, 'beta', density = FALSE)
plot(out.sp.2, 'theta', density = FALSE)

# TODO 4.3 ----------------------------------------------------------------
# Fit the model a third time using the priors specified in
# TODO 4.3 where now you adjust the MCMC settings (e.g., increase
# the number of samples, potentially change the burn-in and thinning rates) to
# ensure that the model converges. Call the model "out.sp.3". Show evidence
# that the model has converged. Which parameters appear to be the parameters
# that are slowest to converge? Note that you will have to run the model for a
# couple of minutes in order to get it to converge.
out.sp.3 <- spPGOcc(occ.formula = ~ scale(elevation) + I(scale(elevation)^2) +
                                    scale(forest),
                    det.formula = ~ scale(date) + I(scale(date)^2) + scale(dur),
                    data = data.goldfinch, inits = inits, priors = priors,
                    tuning = tuning, n.batch = 1200, batch.length = 25,
                    NNGP = TRUE, n.neighbors = 10, n.burn = 10000, n.thin = 40,
                    n.chains = 3, n.report = 100)
summary(out.sp.3)
plot(out.sp.3, 'beta', density = FALSE)
plot(out.sp.3, 'theta', density = FALSE)

# TODO 4.4 ----------------------------------------------------------------
# Use the ppcOcc() function to perform a posterior predictive check that
# assesses the model fit for the spatial model. Use the freeman-tukey test statistic, and perform
# two posterior checks: one that groups the data by site prior to calculating
# the test statistic, and one that groups the data by visit prior to calculating
# the test statistic. See ?ppcOcc for details. Use the summary() function to
# calculate the Bayesian p-value and adapt code shown in the introductory
# package vignette (https://www.jeffdoser.com/files/spoccupancy-web/articles/modelfitting#posterior-predictive-checks)
# to visualize the Bayesian p-value.
# Grouping by sites -------------------
out.ppc.1 <- ppcOcc(out.sp.3, fit.stat = 'freeman-tukey', group = 1)
summary(out.ppc.1)
ppc.df <- data.frame(fit = out.ppc.1$fit.y,
                     fit.rep = out.ppc.1$fit.y.rep,
                     color = 'lightskyblue1')
ppc.df$color[ppc.df$fit.rep > ppc.df$fit] <- 'lightsalmon'
plot(ppc.df$fit, ppc.df$fit.rep, bg = ppc.df$color, pch = 21,
     ylab = 'Fit', xlab = 'True')
lines(ppc.df$fit, ppc.df$fit, col = 'black')
# Grouping by visits ------------------
out.ppc.2 <- ppcOcc(out.sp.3, fit.stat = 'freeman-tukey', group = 2)
ppc.df <- data.frame(fit = out.ppc.2$fit.y,
                     fit.rep = out.ppc.2$fit.y.rep,
                     color = 'lightskyblue1')
ppc.df$color[ppc.df$fit.rep > ppc.df$fit] <- 'lightsalmon'
plot(ppc.df$fit, ppc.df$fit.rep, bg = ppc.df$color, pch = 21,
     ylab = 'Fit', xlab = 'True')
lines(ppc.df$fit, ppc.df$fit, col = 'black')
summary(out.ppc.2)

# TODO 4.5 ----------------------------------------------------------------
# Use waicOcc() to compare the non-spatial model to the spatial model. Which
# model performs better according to the WAIC?
waicOcc(out.non.sp)
waicOcc(out.sp.3)
# The spatial model performs better

# TODO 4.6 ----------------------------------------------------------------
# Use four-fold cross validation to compare the non-spatial model to the spatial
# model. Run the same code you used to fit the two models previously, but now
# specify the arguments needed for k-fold cross-validation. Note that you can set
# the k.fold.only argument to solely perform cross-validation and not re-run the
# entire full model. Which model performs the best?
k.fold.non.sp <- PGOcc(occ.formula = ~ scale(elevation) + I(scale(elevation)^2) + scale(forest),
                       det.formula = ~ scale(date) + I(scale(date)^2) + scale(dur),
                       data = data.goldfinch,
                       n.samples = 5000,
                       n.thin = 4,
                       n.burn = 3000,
                       n.chains = 3,
                       n.report = 500,
                       k.fold = 4,
                       k.fold.threads = 4,
                       k.fold.only = TRUE)
k.fold.sp <- spPGOcc(occ.formula = ~ scale(elevation) + I(scale(elevation)^2) +
                                     scale(forest),
                    det.formula = ~ scale(date) + I(scale(date)^2) + scale(dur),
                    data = data.goldfinch, inits = inits, priors = priors,
                    tuning = tuning, n.batch = 1200, batch.length = 25,
                    NNGP = TRUE, n.neighbors = 10, n.burn = 10000, n.thin = 10,
                    n.chains = 3, n.report = 100, k.fold = 4, k.fold.threads = 4,
                    k.fold.only = TRUE)
k.fold.non.sp$k.fold.deviance
k.fold.sp$k.fold.deviance

# TODO 4.7 ----------------------------------------------------------------
# Using code from "1-swiss-european-goldfinch.R" predict occupancy of the
# European goldfinch across Switzerland. Call your resulting prediction
# object out.pred. Make sure to load in the necessary
# covariate data, which are available in "data/switzerlandPredData.rda". Make
# a species distribution map of the mean occupancy probability across the country.
# Also plot the standard deviation as a map of the uncertainty with our estimates.
# Predict occupancy probability across Switzerland
# Load prediction objects (loads objects pred.swiss and coords.0)
load("data/switzerlandPredData.rda")
str(pred.swiss)
# Standardize elevation and forest prediction values by values used to fit model
elevation.0 <- (pred.swiss[, 'elevation'] - mean(data.goldfinch$occ.covs$elevation)) /
  sd(data.goldfinch$occ.covs$elevation)
forest.0 <- (pred.swiss[, 'forest'] - mean(data.goldfinch$occ.covs$forest)) /
  sd(data.goldfinch$occ.covs$forest)
# Create prediction design matrix
X.0 <- cbind(1, elevation.0, elevation.0^2, forest.0)
# Predict at new locations
# NOTE: this will take a few minutes.
out.pred <- predict(out.sp.3, X.0, coords.0, verbose = TRUE, n.report = 100)
# Occupancy probability means
psi.0.mean <- apply(out.pred$psi.0.samples, 2, mean)
# Occupancy probability standard deviations
psi.0.sd <- apply(out.pred$psi.0.samples, 2, sd)
plot.df <- data.frame(psi.mean = psi.0.mean,
                      psi.sd = psi.0.sd,
                      x = coords.0[, 1],
                      y = coords.0[, 2])
pred.stars <- st_as_stars(plot.df, dims = c('x', 'y'))
psi.mean.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = psi.mean),interpolate = TRUE) +
  scale_fill_gradientn("", colors = ocean.tempo(1000), limits = c(0, 1),
                       na.value = NA) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "Easting", y = "Northing", title = 'Occupancy Mean')
psi.sd.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = psi.sd),interpolate = TRUE) +
  scale_fill_gradientn("", colors = ocean.tempo(1000), limits = c(0, 1),
                       na.value = NA) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "Easting", y = "Northing", title = 'Occupancy SD')
psi.mean.plot + psi.sd.plot

# TODO 4.8 ----------------------------------------------------------------
# In addition to making a map of the predicted occupancy probability, it can
# be very enlightening to plot the resulting predictions of the spatial random
# effects. The posterior predictive samples of the predicted spatial random effects
# are available in "out.pred$w.0.samples". Adapt your code from TODO 4.7 to now
# generate a plot of the spatial random effects across the country.
# Spatial process mean and sd
w.0.mean <- apply(out.pred$w.0.samples, 2, mean)
w.0.sd <- apply(out.pred$w.0.samples, 2, sd)
plot.df <- data.frame(w.mean = w.0.mean,
                      w.sd = w.0.sd,
                      x = coords.0[, 1],
                      y = coords.0[, 2])
pred.stars <- st_as_stars(plot.df, dims = c('x', 'y'))
w.mean.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = w.mean),interpolate = TRUE) +
  scale_fill_gradientn("", colors = ocean.tempo(1000),
                       na.value = NA) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "Easting", y = "Northing", title = 'Spatial Effect Mean')
w.sd.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = w.sd),interpolate = TRUE) +
  scale_fill_gradientn("", colors = ocean.tempo(1000),
                       na.value = NA) +
  theme_bw(base_size = 18) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "Easting", y = "Northing", title = 'Spatial Effect SD')
w.mean.plot + w.sd.plot
