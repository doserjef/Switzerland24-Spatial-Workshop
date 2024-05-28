# 1-swiss-european-goldfinch.R: this script fits a single-species occupancy model
#                               using data on the European Goldfinch
#                               from the Switzerland Breeding Bird Survey
#                               (Swiss MHB) in 2014.
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
set.seed(250)

# In this example, our goal is to produce a species distribution map of the
# European Goldfinch throughout Switzerland. We will compare two models:
# one that uses elevation to predict occupancy probability, and another that
# uses both elevation and forest cover.

# 1. Data prep ------------------------------------------------------------
# Read in the data source (reads in an object called data.goldfinch)
load("data/europeanGoldfinchSwiss.rda")
str(data.goldfinch)
# Take a quick look at the spatial locations
plot(data.goldfinch$coords, pch = 19)

# 2. Model fitting --------------------------------------------------------
# Fit a non-spatial, single-species occupancy model
# Only include linear and quadratic effects of elevation
out.elev <- PGOcc(occ.formula = ~ scale(elevation) + I(scale(elevation)^2),
                  det.formula = ~ scale(date) + I(scale(date^2)) + scale(dur),
                  data = data.goldfinch,
                  n.samples = 5000,
                  n.thin = 4,
                  n.burn = 3000,
                  n.chains = 3,
                  n.report = 500)
summary(out.elev)

# Notice the default priors and initial values specified in the
# "Preparing to run the model" section. Below we refit the same model, but
# now explicitly specify the priors and initial values
prior.list <- list(beta.normal = list(mean = 0, var = 2.72),
                   alpha.normal = list(mean = 0, var = 2.72))
inits.list <- list(beta = 0, alpha = 0,
                   z = apply(data.goldfinch$y, 1, max, na.rm = TRUE))
out.elev.2 <- PGOcc(occ.formula = ~ scale(elevation) + I(scale(elevation)^2),
                    det.formula = ~ scale(date) + I(scale(date^2)) + scale(dur),
                    data = data.goldfinch,
                    n.samples = 5000,
                    inits = inits.list,
                    priors = prior.list,
                    n.thin = 4,
                    n.burn = 3000,
                    n.chains = 3,
                    n.report = 500)
summary(out.elev.2)

# Now fit a model that also includes forest cover
out.full <- PGOcc(occ.formula = ~ scale(elevation) + I(scale(elevation)^2) + scale(forest),
                  det.formula = ~ scale(date) + I(scale(date^2)) + scale(dur),
                  data = data.goldfinch,
                  n.samples = 5000,
                  inits = inits.list,
                  priors = prior.list,
                  n.thin = 4,
                  n.burn = 3000,
                  n.chains = 3,
                  n.report = 500)

# 3. Model validation -----------------------------------------------------
# Perform a posterior predictive check to assess model fit.
ppc.out.elev <- ppcOcc(out.elev, fit.stat = 'freeman-tukey', group = 1)
ppc.out.full <- ppcOcc(out.full, fit.stat = 'freeman-tukey', group = 1)
# Calculate a Bayesian p-value as a simple measure of Goodness of Fit.
# Bayesian p-values between 0.1 and 0.9 indicate adequate model fit.
summary(ppc.out.elev)
summary(ppc.out.full)

# Simple plot to visualize the Bayesian p-value
ppc.df <- data.frame(fit = ppc.out.full$fit.y,
                     fit.rep = ppc.out.full$fit.y.rep,
                     color = 'lightskyblue1')
ppc.df$color[ppc.df$fit.rep > ppc.df$fit] <- 'lightsalmon'
plot(ppc.df$fit, ppc.df$fit.rep, bg = ppc.df$color, pch = 21,
     ylab = 'Fit', xlab = 'True')
abline(0, 1)

# Plot difference in discrepancy measure between the replicate and actual
# data across each of the sites
diff.fit <- ppc.out.full$fit.y.rep.group.quants[3, ] -
            ppc.out.full$fit.y.group.quants[3, ]
plot(diff.fit, pch = 19, xlab = 'Site ID', ylab = 'Replicate - True Discrepancy')
# No massive outliers. We could also make a plot across space to highlight variation
# in the effects.

# 4. Model comparison -----------------------------------------------------
# Compute Widely Applicable Information Criterion (WAIC)
# Lower values indicate better model fit.
# Model with elevation only
waicOcc(out.elev)
# Model with elevation and forest cover
waicOcc(out.full)
# Incorporating forest cover results in a substantial improvement in WAIC.

# 5. Posterior summaries --------------------------------------------------
# Concise summary of main parameter estimates
summary(out.full)
# Take a look at objects in resulting object
names(out.full)
# MCMC samples for the different parameters are stored in *.samples objects.
str(out.full$beta.samples)
# Create simple plot summaries using MCMCvis package.
# Occupancy covariate effects ---------
MCMCplot(out.full$beta.samples, ref_ovl = TRUE, ci = c(50, 95))
# Detection covariate effects ---------
MCMCplot(out.full$alpha.samples, ref_ovl = TRUE, ci = c(50, 95))
# Make a conditional effects plot showing relationship between occupancy
# probability and elevation -----------
# Create a set of values across the range of observed elevation values
elev.pred.vals <- seq(min(data.goldfinch$occ.covs$elev),
                      max(data.goldfinch$occ.covs$elev),
                      length.out = 100)

# Scale predicted values by mean and standard deviation used to fit the model
elev.pred.vals.scale <- (elev.pred.vals - mean(data.goldfinch$occ.covs$elev)) /
                        sd(data.goldfinch$occ.covs$elev)
# Predict occupancy across elev values at mean values of all other variables
pred.df <- as.matrix(data.frame(intercept = 1, elev = elev.pred.vals.scale,
                                elev.2 = elev.pred.vals.scale^2, forest = 0))
out.pred <- predict(out.full, pred.df)
str(out.pred)
psi.0.quants <- apply(out.pred$psi.0.samples, 2, quantile,
                      prob = c(0.025, 0.5, 0.975))
psi.plot.dat <- data.frame(psi.med = psi.0.quants[2, ],
                           psi.low = psi.0.quants[1, ],
                           psi.high = psi.0.quants[3, ],
                           elev = elev.pred.vals)
ggplot(psi.plot.dat, aes(x = elev, y = psi.med)) +
  geom_line(aes(x = elev, y = psi.low), linewidth = 0.9, lty = 2, col = 'gray') +
  geom_line(aes(x = elev, y = psi.high), linewidth = 0.9, lty = 2, col = 'gray') +
  # Or could use a ribbon if you like that better
  # geom_ribbon(aes(ymin = psi.low, ymax = psi.high), fill = 'grey70') +
  geom_line() +
  theme_bw() +
  scale_y_continuous(limits = c(0, 1)) +
  labs(x = 'Elevation (meters)', y = 'Occupancy Probability')

# 6. Prediction -----------------------------------------------------------
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
out.pred <- predict(out.full, X.0)
# Occupancy probability means
psi.0.mean <- apply(out.pred$psi.0.samples, 2, mean)
# Occupancy probability standard deviations
psi.0.sd <- apply(out.pred$psi.0.samples, 2, sd)

# Create a species distribution map with uncertainty ----------------------
plot.df <- data.frame(psi.mean = psi.0.mean,
                      psi.sd = psi.0.sd,
                      x = coords.0[, 1],
                      y = coords.0[, 2])
pred.stars <- st_as_stars(plot.df, dims = c('x', 'y'))
psi.mean.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = psi.mean),interpolate = TRUE) +
  scale_fill_gradientn("", colors = ocean.tempo(1000), limits = c(0, 1),
                       na.value = NA) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "Easting", y = "Northing", title = 'Occupancy Mean')
psi.sd.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = psi.sd),interpolate = TRUE) +
  scale_fill_gradientn("", colors = ocean.tempo(1000),
                       na.value = NA) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "Easting", y = "Northing", title = 'Occupancy SD')
psi.mean.plot + psi.sd.plot
