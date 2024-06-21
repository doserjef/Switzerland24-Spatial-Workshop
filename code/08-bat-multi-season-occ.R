# 08-bat-multi-season-occ.R: this script fits multiple multi-season occupancy models
#                            to estimate the distribution of the spotted bat across
#                            Oregon and Washington, USA.
# Data source citation:
#   Wright, W. J., Irvine, K. M., Rodhouse, T. J., & Litt, A. R. (2021).
#   Spatial Gaussian processes improve multi‐species occupancy models when range
#   boundaries are uncertain and nonoverlapping. Ecology and Evolution, 11(13), 8516-8527.
#
#   The data from the above citation come from the North American Bat Monitoring
#   Program.

rm(list = ls())
library(spOccupancy)
# For making plots
library(ggplot2)
library(dplyr)
library(sf)
library(patchwork)
set.seed(89421)

# Our goal here focuses primarily on prediction in which we seek to predict the
# distribution of the spotted bat (Euderma maculatum) across the states of Oregon and
# Washington, USA using data collected over a four year period from 2016-2019. We
# will compare a variety of multi-season occupancy models in accomplishing this task.
# Note that the models we fit here differ from those used in the actual publication,
# which fits a slightly different form of multi-season occupancy model (and do so
# in a multi-species framework). We will fit and compare the four models that we discussed in lecture:
# (1) Unstructured site random effect and unstructured year random effect
# (2) Unstructured site random effect and AR(1) year random effect
# (3) Spatial random effect and unstructured year random effect
# (4) Spatial random effect and AR(1) year random effect

# All models will include the same covariates on occupancy and detection, which are
# as follows
# Occupancy covariates (all static over time):
#    1. Log forest cover (standardized across study region)
#    2. 30-year average precipitation (standardized across study region)
#    3. Cliff and canyon cover (standardized across study region)
# Detection covariates:
#    1. Vegetation clutter
#    2. Binary variable indicating whether the sampling location was near water or not.
#    3. Average minimum temperature

# Load the data set -------------------------------------------------------
# Loads an spOccupancy formatted data list
load('data/wright2021BatData.rda')
str(data.list)
# Note that the detection-nondetection data are currently formatted for a
# multi-species, multi-season occupancy model. For those models, the detection
# nondetection data y are stored in a four-dimensional array with dimensions
# corresponding to species, site, primary time period, and secondary time period.
# Below, we select a single species (the spotted bat, EUMA) that we will work with here.
sp.indx <- which(dimnames(data.list$y)[[1]] == 'euma')
data.spotted.bat <- data.list
data.spotted.bat$y <- data.list$y[sp.indx, , , ]
str(data.spotted.bat)
# For single-species multi-season occupancy models the detection-nondetection data
# are formatted in a three-dimensional array with dimensions corresponding to
# site, primary time period (in this case year), and secondary time period.

# Exploratory data analysis -----------------------------------------------
# The detection-nondetection data come from 10 x 10km grid cells across the
# two states. Below we load the prediction data, which contains covariate information
# for grid cells across the entire region as well as an sf object that we will
# use to make a quick plot showing the grid cells that we have sampled across our
# region of interest.

# Loads three objects: (1) X.pred (design matrix for prediction); (2) coords.pred
# (grid cell coordinates for prediction); and (3) nw.grid (an sf object for
# generating maps of the region)
load('data/wright2021BatPredictionData.rda')
# The nw.grid$samp_all column indicates the grid cells that were sampled at least
# once over the four year period.
ggplot() +
  geom_sf(data = nw.grid, aes(fill = factor(samp_all))) +
  theme_bw() +
  scale_fill_viridis_d() +
  labs(fill = 'Sampled')
# Sample locations were determined based on a spatially balanced random sample of
# the plotted grid-based sampling scheme.

# Fit models --------------------------------------------------------------
# Next we fit the four candidate models. The models with an unstructured site
# random effect (non-spatial models) will be fit with tPGOcc(). The spatial models
# are fit with stPGOcc(). Note that ALL multi-season models use "n.batch" and
# "batch.length" arguments to specify the MCMC samples. This is because the
# AR1 temporal autocorrelation parameter is also a "difficult to sample" parameter,
# and so we need to use the Adaptive MCMC sampler we have discussed previously.

# MCMC criteria we will use for all four models
n.batch <- 1200
batch.length <- 25
n.burn <- 15000
n.thin <- 15
n.chains <- 3

# Note that when running all four candidate models below, we have included the
# k.fold and k.fold.threads arguments to perform four-fold cross-validation using
# four threads.

# Unstructured site and year ----------
# Specify slightly informative priors on the site and year random effects to restrict
# extremely large values. These priors give a prior mean of 0.5 to the random effect
# variances with an infinite variance, and so they are still fairly weakly informative,
# but perhaps a bit less so than the default values of 0.1 and 0.1 for the inverse-Gamma
# shape and scale parameters.
priors <- list(sigma.sq.psi.ig = list(a = 2, b = 0.5))
# Approx run time: 2 min
out.1 <- tPGOcc(occ.formula = ~ log.forest.s + ppt.s + cliff.s + (1 | site.id) +
                                (1 | year),
                det.formula = ~ scale(clutter) + scale(water) + scale(tmin),
                data = data.spotted.bat, priors = priors,
                n.batch = n.batch, batch.length = batch.length, ar1 = FALSE,
                n.burn = n.burn, n.thin = n.thin,
                n.chains = n.chains, n.report = 100, k.fold = 4,
                k.fold.threads = 4)
summary(out.1)
# Unstructured site and AR1 -----------
# Approx run time: 2 min
out.2 <- tPGOcc(occ.formula = ~ log.forest.s + ppt.s + cliff.s + (1 | site.id),
                det.formula = ~ scale(clutter) + scale(water) + scale(tmin),
                data = data.spotted.bat, priors = priors,
                n.batch = n.batch, batch.length = batch.length, ar1 = TRUE,
                n.burn = n.burn, n.thin = n.thin,
                n.chains = n.chains, n.report = 100,
                k.fold = 4, k.fold.threads = 4)
summary(out.2)

# Spatial and unstructured year -------
# Approx run time: 4.5 min
out.3 <- stPGOcc(occ.formula = ~ log.forest.s + ppt.s + cliff.s + (1 | year),
                 det.formula = ~ scale(clutter) + scale(water) + scale(tmin),
                 data = data.spotted.bat, priors = priors,
                 n.batch = n.batch, batch.length = batch.length,
                 cov.model = 'exponential', NNGP = TRUE, ar1 = FALSE,
                 n.neighbors = 8, n.burn = n.burn, n.thin = n.thin,
                 n.chains = n.chains, n.report = 100, k.fold = 4, k.fold.threads = 4)
summary(out.3)
# Spatial and AR1 ---------------------
# Approx run time: 5.5 min
out.4 <- stPGOcc(occ.formula = ~ log.forest.s + ppt.s + cliff.s,
                 det.formula = ~ scale(clutter) + scale(water) + scale(tmin),
                 data = data.spotted.bat, priors = priors,
                 n.batch = n.batch, batch.length = batch.length,
                 cov.model = 'exponential', NNGP = TRUE, ar1 = TRUE,
                 n.neighbors = 8, n.burn = n.burn, n.thin = n.thin,
                 n.chains = n.chains, n.report = 100, k.fold = 4,
                 k.fold.threads = 4)
summary(out.4)

# Model comparison -------------------
# Four-fold cross-validation
out.1$k.fold.deviance
out.2$k.fold.deviance
out.3$k.fold.deviance
out.4$k.fold.deviance
# WAIC
waicOcc(out.1)
waicOcc(out.2)
waicOcc(out.3)
waicOcc(out.4)

# What models are the top performing models? Why do you think the differences
# are more drastic for the cross-validation as opposed to WAIC?

# Prediction --------------------------------------------------------------
# Let's generate predictions for two models: Model 2 (unstructured site and AR1) and
# Model 4 (spatial and AR1).
# The objects we need for prediction are stored in X.pred (design matrix of covariates
# across all prediction cells) and coords.pred (the prediction coordinates)
str(X.pred)
str(coords.pred)
# Notice that X.pred is a three-dimensional array, with dimensions corresponding
# to the prediction sites, the prediction seasons, and the covariates

# Model 2 (unstructured site, AR1) ----
# First, we need to subset the prediction array to only include the covariates
# in the model
dimnames(X.pred)
# Remove year from X.pred
X.0.2 <- X.pred[, , c('(Intercept)', 'log.forest.s', 'ppt.s', 'cliff.s', 'site.id')]
# Specify t.cols, which indicates which of the primary time periods you want to
# predict for. In our case, we want to predict across all four years.
t.cols <- 1:dim(data.spotted.bat$y)[2]
out.pred <- predict(out.2, X.0 = X.0.2, t.cols = t.cols)
str(out.pred)
# Get the 2.5%, 50% (median), and 97.5% quantiles of the occupancy probability
# predictions in each grid cell
psi.0.quants.non.sp <- apply(out.pred$psi.0.samples, c(2, 3), quantile, c(0.025, 0.5, 0.975))

# Model 4 (spatial, AR1) --------------
# Subset design matrix to get covariates of interest
X.0.4 <- X.pred[, , c('(Intercept)', 'log.forest.s', 'ppt.s', 'cliff.s')]
# Specify t.cols, which indicates which of the primary time periods you want to
# predict for. In our case, we want to predict across all four years.
t.cols <- 1:dim(data.spotted.bat$y)[2]
# Note that we set n.omp.threads = 4 to do the prediction in parallel. Prediction
# in particular can be sped up fairly substantially by using this argument.
out.pred <- predict(out.4, X.0 = X.0.4, t.cols = t.cols, coords = coords.pred,
                    n.report = 100, verbose = TRUE,
                    n.omp.threads = 4)
str(out.pred)
# Get the 2.5%, 50% (median), and 97.5% quantiles of the occupancy probability
# predictions in each grid cell for each year
psi.0.quants.sp <- apply(out.pred$psi.0.samples, c(2, 3), quantile, c(0.025, 0.5, 0.975))

# Generate some species distribution maps ---------------------------------
# Number of years
n.years <- length(t.cols)
# Number of prediction sites
J <- nrow(X.0.4)
# Put things in a data frame for plotting with sf and ggplot.
plot.df <- data.frame(geometry = rep(nw.grid$geometry, times = n.years),
                      median.sp = c(psi.0.quants.sp[2, , ]),
                      ci.width.sp = c(psi.0.quants.sp[3, , ] - psi.0.quants.sp[1, , ]),
                      median.non.sp = c(psi.0.quants.non.sp[2, , ]),
                      ci.width.non.sp = c(psi.0.quants.non.sp[3, , ] - psi.0.quants.non.sp[1, , ]),
                      year = rep(2016:2019, each = J))
# Make it an sf object
plot.sf <- st_as_sf(plot.df)

# Plot median occupancy probability in each year
non.sp.plots <- ggplot() +
  geom_sf(data = plot.sf, aes(fill = median.non.sp)) +
  scale_fill_viridis_c("Occupancy\nprobability", limits = c(0, 1)) +
  theme_bw() +
  facet_wrap(vars(year))
non.sp.plots
sp.plots <- ggplot() +
  geom_sf(data = plot.sf, aes(fill = median.sp)) +
  scale_fill_viridis_c("Occupancy\nprobability", limits = c(0, 1)) +
  theme_bw() +
  facet_wrap(vars(year))
sp.plots

# Plot median and 95% CI width of occupancy probability in 2019 for both models
non.sp.med.plot <- plot.sf %>%
  dplyr::filter(year == 2019) %>%
  ggplot() +
    geom_sf(aes(fill = median.non.sp)) +
    scale_fill_viridis_c("Occupancy\nprobability", limits = c(0, 1)) +
    theme_bw() +
    labs(title = '(a) Non-spatial occupancy median')
non.sp.ci.plot <- plot.sf %>%
  dplyr::filter(year == 2019) %>%
  ggplot() +
    geom_sf(aes(fill = ci.width.non.sp)) +
    scale_fill_viridis_c("95% CI\nWidth", limits = c(0, 1)) +
    theme_bw() +
    labs(title = '(b) Non-spatial occupancy 95% CI width')
sp.med.plot <- plot.sf %>%
  dplyr::filter(year == 2019) %>%
  ggplot() +
    geom_sf(aes(fill = median.sp)) +
    scale_fill_viridis_c("Occupancy\nprobability", limits = c(0, 1)) +
    theme_bw() +
    labs(title = '(c) Spatial occupancy median')
sp.ci.plot <- plot.sf %>%
  dplyr::filter(year == 2019) %>%
  ggplot() +
    geom_sf(aes(fill = ci.width.sp)) +
    scale_fill_viridis_c("95% CI\nWidth", limits = c(0, 1)) +
    theme_bw() +
    labs(title = '(d) Spatial occupancy 95% CI width')

# Put them all together
non.sp.med.plot + non.sp.ci.plot + sp.med.plot + sp.ci.plot

# What differences are there when generating predictions with the spatial and
# non-spatial models?
