# 8-wood-thrush-spatial-trend-occ.R: this script fits a multi-season
#                                    spatially-varying occupancy model to estimate
#                                    a spatially varying trend in occupancy
#                                    probability of the wood thrush across the
#                                    eastern US from 2000-2009. The data used
#                                    here are a subset of data from the USGS
#                                    breeding bird survey.
# Data source citation:
#    Ziolkowski Jr., D.J., Lutmerding, M., English, W.B., Aponte, V.I.,
#    and Hudson, M-A.R., 2023, North American Breeding Bird Survey
#    Dataset 1966 - 2022: U.S. Geological Survey data release, https://doi.org/10.5066/P9GS9K64.

rm(list = ls())
library(spOccupancy)
library(ggplot2)
library(sf)
library(stars)
set.seed(3420)

# Our objectives in this exercise are to generate a spatially-explicit trend map
# of occupancy probability from 2000-2009 for the wood thrush across the
# eastern US. Note we will keep this exercise fairly short and will only fit a
# single model, but it is often useful to compare SVC models to simpler models
# to gauge how much support there is for including the spatially-varying
# coefficient (which in this case is a trend).

# Load the data set -------------------------------------------------------
# Loads an object called data.WOTH. Note that this is in the exact format
# required for multi-season spOccupancy models
load('data/wood-thrush-occ-bbs-data.rda')
str(data.WOTH)

# Exploratory data analysis -----------------------------------------------
# Plot of raw change in naive occurrence probability over the 10 year period
# Mean naive occurrence probability each year
y.naive.means <- apply(data.WOTH$y, 2, mean, na.rm = TRUE)
# Years in the data set
years <- 2000:2009
plot(years, y.naive.means, xlab = "Year",
     ylab = "Naive Occurrence Probability", pch = 19)
# Occupancy probability appears to be declining. If you're familiar with the wood thrush,
# this is not surprising, as it's populations are known to be substantially declining
# across most of the eastern US.

# Plot the number of times WOTH was detected at each site over the 10 years
# Proj4string for coordinate system
my.proj <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=km +no_defs"
# Create coordinates as sf object for plotting
coords.sf <- st_as_sf(as.data.frame(data.WOTH$coords),
                      coords = c('X', 'Y'),
                      crs = my.proj)
# Total number of years each site was sampled
coords.sf$y.sums <- apply(apply(data.WOTH$y, c(1, 2), max), 1, sum)
# Get map of the eastern US
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
usa.bbox <- st_bbox(usa)
usa.bbox[1] <- -100
usa.bbox <- as.vector(usa.bbox)
sf_use_s2(FALSE)
east.us <- st_crop(st_make_valid(usa), xmin = usa.bbox[1], ymin = usa.bbox[2],
                   xmax = usa.bbox[3], ymax = usa.bbox[4])
east.us <- east.us %>%
  st_transform(st_crs(coords.sf))
# Make a simple plot
ggplot() +
  geom_sf(data = east.us, alpha = 0, col = 'black') +
  geom_sf(data = coords.sf, aes(col = y.sums)) +
  scale_color_viridis_c() +
  theme_bw() +
  labs(col = '# of years\ndetected')

# Specify prior distributions ---------------------------------------------
# Here we will use a somewhat informative prior for the spatially-varying trend.
# We will restrict the lower bound of the effective spatial range of the
# spatially-varying trend to be 3 / q(20), where q(20) represents the
# distance of the 20% quantile across all the sites in the data set. In this case
# that corresponds to about 470km. Such a prior is often useful in SVC models
# when interest focuses on estimating the spatially-varying trend (or SVC for
# some other type of covariate), but you want to avoid extremely fine-scale spatial
# variation in the effect.

# Calculate distance matrix for all sites
dist.bbs <- dist(data.WOTH$coords)
# Get the 20 percent quantile of the distances to restrict the effective spatial
# range for the spatially-varying trend.
low.dist <- quantile(dist.bbs, 0.20)
low.dist
# The minimum distance, which we will use for the prior on the spatially-varying
# intercept.
min.dist <- min(dist.bbs)
min.dist
# Maximum distance
max.dist <- max(dist.bbs)
max.dist
# Specify the prior distributions
priors <- list(sigma.sq.ig = list(2, 1), # IG prior for spatial variances
               # Weakly informative prior for the SVI, more informative prior
               # for the SVC trend decay parameters.
               phi.unif = list(c(3 / max.dist, 3 / max.dist),
                               c(3 / min.dist, 3 / low.dist)),
               # Normal prior for occurrence coefficients
               beta.normal = list(mean = 0, var = 2.72),
               # Normal prior for observational (detection) coefficients
               alpha.normal = list(mean = 0, var = 2.72))

# Specify initial values and tuning variances -----------------------------
# SVC models can be sensitive to initial values, so it is often good to specify
# initial values instead of rely on the defaults.
# Here setting initial effective spatial ranges to be 500km for the intercept
# and 1000km for the trend.
# Initial values for the spatial variances uses a substantially higher value
# for the intercept than the trend.
inits <- list(phi = c(3 / 500, 3 / 1000),
              sigma.sq = c(5, 0.5),
              beta = c(1.5, 0), alpha = 0)
# We can also specify the tuning values for the spatial decay parameters. The defaults
# for this initial tuning variance is 1. Usually the more sites you're working with,
# the smaller the ideal tuning variance will be (but this is not always the case).
tuning.list <- list(phi = c(0.6, 0.5))

# Fit the model -----------------------------------------------------------
# MCMC criteria
n.batch <- 1200
batch.length <- 25
n.burn <- 20000
n.thin <- 5
# Only fitting 1 chain to speed things up a bit. Of course in a complete
# analysis we would want to fit the model for three chains.
n.chains <- 1

# Note that the only new argument in svcTPGOcc compared to stPGOcc is the
# svc.cols argument in which you specify the columns for which you want to include
# a spatially varying effect. This can be done by specifying the name of the
# variable themselves directly as how it is included in the model formula, or
# by specifying the corresponding number of the coefficient in order in the model
# formula (the intercept is always 1 unless explicitly fitting a model without
# an intercept).
# Approx run time for 1 chain: 3 min
out <- svcTPGOcc(occ.formula = ~ scale(years),
                 det.formula = ~ scale(day) + I(scale(day)^2) +
                                 scale(year.det) + I(scale(year.det)) +
                                 scale(rep.val) + I(scale(rep.val)^2),
                 # svc.cols = c(1, 2),
                 # Alternatively, an equivalent specification of svc.cols is
                 svc.cols = c('(Intercept)', 'scale(years)'),
                 data = data.WOTH,
                 inits = inits,
                 priors = priors,
                 n.batch = n.batch,
                 batch.length = batch.length,
                 accept.rate = 0.43,
                 cov.model = "exponential",
                 tuning = tuning.list,
                 n.omp.threads = 1,
                 verbose = TRUE,
                 NNGP = TRUE,
                 n.neighbors = 5,
                 n.report = 100,
                 n.burn = n.burn,
                 n.thin = n.thin,
                 n.chains = n.chains)

# Look at model results.
summary(out)
# NOTE: this model has not converged and we would want to run it longer. For context,
#       in our experience with SVC multi-season occupancy models, we will often need
#       to run models for at least 100,000 MCMC iterations, and often will manually
#       specify initial values as well.

# Predict trend across entire eastern US ----------------------------------
load('data/bbs-occ-pred-coordinates.rda')
# Design "matrix" should be a 3-D array with dimensions corresponding to site, year, and parameter
X.0 <- array(NA, dim = c(nrow(coords.0), ncol(out$y), dim(out$X)[3]))
# Intercept
X.0[, , 1] <- 1
# Year (make sure to standardize it accordingly)
unique.years <- unique(c(data.WOTH$occ.covs$years))
for (t in 1:ncol(out$y)) {
  X.0[, t, 2] <- unique.years[t]
}
X.0[, , 2] <- (X.0[, , 2] - mean(c(data.WOTH$occ.covs$years))) / sd(c(data.WOTH$occ.covs$years))
# The t.cols argument specifies that we are predicting for all ten years. Note our use
# of the n.omp.threads argument to speed up prediction a little bit.
# Note that this takes a substantial amount of memory, so there is a chance
# R aborts if you don't have enough RAM. If this happens, increase the thinning rate
# when fitting the model to save a smaller number of MCMC samples (e.g.,
# if you set n.thin = 10 instead of n.thin = 5).
out.pred <- predict(out, X.0, coords.0, t.cols = 1:10, n.report = 100,
                    verbose = TRUE, n.omp.threads = 4)

# Get the full spatially-varying coefficients
# The getSVCSamples() function can be used to extract the full spatially-varying
# coefficient at either the model fit locations, or the prediction locations. If you
# want to get the values at the prediction locations, you should provide the model
# fit object and the model prediction object. Note that this returns MCMC samples
# for the sum of the non-spatial (beta.samples) and spatial (w.samples) components
# for each SVC (and SVI).
svc.pred.samples <- getSVCSamples(out, out.pred)
str(svc.pred.samples)

plot.df <- data.frame(trend = apply(svc.pred.samples[[2]], 2, mean),
                      trend.prob.pos = apply(svc.pred.samples[[2]], 2, function(a) mean(a > 0)),
                      x = coords.0[, 1],
                      y = coords.0[, 2])
pred.stars <- st_as_stars(plot.df, dims = c('x', 'y'))
coords.sf <- st_as_sf(as.data.frame(data.WOTH$coords),
                      coords = c("X", "Y"),
                      crs = my.proj)
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
usa.no.states <- st_as_sf(maps::map("usa", fill = TRUE, plot = FALSE))
# Restrict to east of the 100th meridian
usa.bbox <- st_bbox(usa)
usa.bbox[1] <- -100
usa.bbox <- as.vector(usa.bbox)
sf_use_s2(FALSE)
# Full data
east.us <- st_crop(st_make_valid(usa), xmin = usa.bbox[1], ymin = usa.bbox[2],
                   xmax = usa.bbox[3], ymax = usa.bbox[4])
east.us <- east.us %>%
  st_transform(st_crs(coords.sf))

# Map of the mean effect
ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = trend), interpolate = FALSE) +
  geom_sf(data = east.us, alpha = 0, col = 'grey') +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC',
                       na.value = NA) +
  theme_bw(base_size = 14) +
  labs(x = "Longitude", y = "Latitude", fill = "", title = 'Wood Thrush Trend (2000-2009)') +
  theme(legend.position = c(0.87, 0.25),
        legend.background = element_rect(fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.title = element_text(size = 14),
        plot.title = element_text(size = 18),
        legend.text = element_text(size = 12))
# We see predominately negative trends across the region, but see a few pockets of
# locations where occupancy of the wood thrush may be increasing.
