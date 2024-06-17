# 09a-solutions-spatial-glmm-northern-cardinal.R: script to fit a variety of spatial and nonspatial GLMMs
#                                                 to estimate relative abundance of the Northern Cardinal in
#                                                 North Carolina, USA. These data come from the USGS
#                                                 North American Breeding Bird Survey. This script
#                                                 contains the solutions.
# Data source citation:
#    Ziolkowski Jr., D.J., Lutmerding, M., English, W.B., Aponte, V.I.,
#    and Hudson, M-A.R., 2023, North American Breeding Bird Survey
#    Dataset 1966 - 2022: U.S. Geological Survey data release, https://doi.org/10.5066/P9GS9K64.

rm(list = ls())
library(spAbundance)
library(ggplot2)
library(sf)
library(stars)
library(dplyr)
set.seed(293082)

# In this exercise, our objective is to get an understanding of how relative
# abundance of the Northern Cardinal varies across the state of North Carolina, USA.

# Load the data set -------------------------------------------------------
# Loads an object called data.list that is formatted for use in multi-species
# GLMMs in spAbundance.
load('data/nc-bbs-data.rda')
str(data.list)
# Notice the similarity in format compared to spOccupancy.
#    - y is a two-dimensional matrix with rows corresponding to species and
#      columns corresponding to sites.
#    - covs is a data frame or matrix of covariates for inclusion in the model.
#      Note that with all GLMM functions in spAbundance, there is no separation
#      of detection covariates from abundance covariates.
#    - coords is a matrix or data frame of the spatial coordinates for each site.

# Exploratory data analysis -----------------------------------------------
# Let's make a quick map of the data points across the state of North Carolina.

# This is the coordinate system used for the coordinates
my.crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"
# Convert coordinates to an sf object for plotting
coords.sf <- st_as_sf(as.data.frame(data.list$coords),
                      coords = c('X', 'Y'),
                      crs = my.crs)
# Get a map of NC used for the plots
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
nc <- usa %>%
  dplyr::filter(ID == 'north carolina') %>%
  st_transform(crs = my.crs)
# Make the plot.
ggplot() +
  geom_sf(data = nc, alpha = 0, col = 'black') +
  geom_sf(data = coords.sf, col = 'black') +
  theme_bw()
# Seems like we have a fairly good distribution of points across the state.

# Northern Cardinal (NOCA) Analysis ---------------------------------------
# First we will focus on estimating relative abundance of the Northern Cardinal
# (NOCA). We will compare four models that vary in whether they
# incorporate a spatial random effect, and whether they use a Poisson distribution
# or a negative binomial distribution. The following few lines of code
# subsets the multi-species data.list into a new object (data.NOCA) for use
# with modelling relative abundance of NOCA.
sp.names <- rownames(data.list$y)
# Determine the row in data.list$y that corresponds to Northern Cardinal
sp.indx <- which(sp.names == 'NOCA')
# Create a new data list that only contains the data for NOCA
data.NOCA <- data.list
data.NOCA$y <- data.list$y[sp.indx, ]
str(data.NOCA)

# TODO 9.1: Plot a histogram of the relative abundance count values for NOCA.
#           What do you notice about the distribution of counts? Based on this
#           distribution, do you think a Poisson or Negative Binomial model
#           will be preferred?

hist(data.NOCA$y)
# First, note the counts of NOCA are all quite high, indicating it is a pretty
# common species. Second, the distribution appears fairly normally distributed
# without much overdispersion, which could indicate that a Poisson distribution
# might be adequate for these data.

# TODO 9.2: Fit a non-spatial GLMM using the abund() function and a Poisson
#           distribution to the Northern Cardinal data set. Include the
#           following covariates in the model: tree canopy cover (tcc, linear),
#           elevation (elev, both linear and quadratic), precipitation (ppt),
#           day of year (day, both linear and quadratic), and a random effect
#           of observer (1 | obs). Note that the day of year and the random observer
#           effect are used to account for variation in detection probability.
#           Standardize all covariates to have a mean of 0 and sd of 1.
#           Ensure the model converges. You can use the default initial values,
#           priors, and tuning variances.
# MCMC criteria
n.batch <- 2000
batch.length <- 25
n.burn <- 30000
n.thin <- 20
n.chains <- 3

# Model 1: Non-spatial Poisson --------
out.non.sp.p <- abund(formula = ~ scale(tcc) + scale(elev) + I(scale(elev)^2) +
                                  scale(ppt) + scale(day) + I(scale(day)^2) +
                                  (1 | obs),
                      data = data.NOCA, family = 'Poisson', n.batch = n.batch,
                      batch.length = batch.length, n.burn = n.burn,
                      n.thin = n.thin, n.chains = n.chains, n.report = 400)
summary(out.non.sp.p)

# TODO 9.3: Fit the same non-spatial GLMM as in TODO 9.2 (same covariates),
#           but now use a Negative Binomial distribution. Ensure the model converges.
#           Default initial values, priors, and tuning variances should be adequate.
#           Look at the estimated value for the negative binomial overdispersion
#           parameter. Given the estimated value, do you think the Poisson or NB
#           will be favored by model selection criteria?

# Model 2: Non-spatial NB -------------
out.non.sp.nb <- abund(formula = ~ scale(tcc) + scale(elev) + I(scale(elev)^2) +
                                   scale(ppt) + scale(day) + I(scale(day)^2) +
                                   (1 | obs),
                       data = data.NOCA, family = 'NB', n.batch = n.batch,
                       batch.length = batch.length, n.burn = n.burn,
                       n.thin = n.thin, n.chains = n.chains, n.report = 400)
summary(out.non.sp.nb)

# TODO 9.4: Fit a spatial GLMM using the spAbund() function and a Poisson
#           distribution. Use the same covariates as in the previous two models.
#           Fit the model with 15 nearest neighbors in the NNGP approximation,
#           and use an exponential covariance function. Default inits, priors,
#           and tuning variances should again be adequate.

# Model 3: Spatial Poisson ------------
out.sp.p <- spAbund(formula = ~ scale(tcc) + scale(elev) + I(scale(elev)^2) +
                                scale(ppt) + scale(day) + I(scale(day)^2) +
                                (1 | obs),
                    data = data.NOCA, family = 'Poisson', n.neighbors = 15,
                    cov.model = 'exponential', n.batch = n.batch,
                    batch.length = batch.length, n.burn = n.burn,
                    n.thin = n.thin, n.chains = n.chains, n.report = 400)
summary(out.sp.p)

# TODO 9.5: What is the median estimate for the effective spatial range?
#           What is the 95% credible interval for the effective spatial range?
#           Recall the default prior for the spatial decay parameter
#           is Unif(3 / max.dist, 3 / min.dist), where max.dist and min.dist
#           are the max and min intersite distances. Calculate the distance
#           matrix and determine these maximum and minimum values. Note
#           that the coordinates in data.NOCA are specified in km.
#           What do you notice about the estimates of the effective spatial range
#           relative to the potential values specified by the prior distribution?
#           Why might this be the case?
esr.samples <- 3 / out.sp.p$theta.samples[, 2]
quantile(esr.samples, c(0.025, 0.5, 0.975))
dist.mat <- dist(data.NOCA$coords)
range(dist.mat)

# TODO: 9.6: As our last model for the NOCA example, fit a spatial negative
#            binomial GLMM with the spAbund function. Again use 15 nearest
#            neighbors and the same covariates as you have in the previous
#            three models.
# Model 4: Spatial NB -----------------
out.sp.nb <- spAbund(formula = ~ scale(tcc) + scale(elev) + I(scale(elev)^2) +
                       scale(ppt) + scale(day) + I(scale(day)^2) +
                       (1 | obs),
                     data = data.NOCA, family = 'NB', n.neighbors = 15,
                     cov.model = 'exponential', n.batch = n.batch,
                     batch.length = batch.length, n.burn = n.burn,
                     n.thin = n.thin, n.chains = n.chains, n.report = 400)
summary(out.sp.nb)

# TODO 9.7: Use the waicAbund() function to calculate the WAIC for the four
#           models we fit previously. Which model is the top performing model?
waicAbund(out.non.sp.p)
waicAbund(out.non.sp.nb)
waicAbund(out.sp.p)
waicAbund(out.sp.nb)
# Spatial Poisson model is the best.

# TODO 9.8: Each of the model output lists contains the "y.rep.samples" object,
#           which is an array of MCMC samples for the "fitted" values in our
#           model. These can be easily extracted from your saved model
#           object with the fitted() function. Using the spatial Poisson model,
#           extract the posterior MCMC samples with fitted(), calculate the mean
#           estimate for the fitted value at each of the 77 sites using the apply()
#           function, and then plot the true data points stored in data.NOCA$y
#           vs the estimated fitted values as a simple visual assessment
#           of the model fit. What do you notice?

fitted.out.sp.p <- fitted(out.sp.p)
str(fitted.out.sp.p)
y.rep.means <- apply(fitted.out.sp.p, 2, mean)
plot(data.NOCA$y, y.rep.means, pch = 19, xlab = 'True', ylab = 'Fitted')
abline(0, 1)
# Close correspondence between the observed and fitted value.

# TODO 9.9: Use the ppcAbund function with the top performing model to perform a
#           posterior predictive check. When working with count data, we are not
#           required to group the data points prior to performing the posterior
#           predictive check like we are when working with binary data. Perform a posterior
#           predictive check for the top performing model with a chi-squared fit
#           calculated based on the raw counts (i.e., group = 0). Use the
#           summary function to calculate a Bayesian p-value.

ppc.out.sp.p <- ppcAbund(out.sp.p, fit.stat = 'chi-squared', group = 0)
summary(ppc.out.sp.p)

# TODO 9.10: Use the prediction data loaded below to predict Northern Cardinal
#            abundance across the state using the top performing model. The prediction
#            data set consists of a grid of locations across North Carolina. Note that
#            when forming the design matrix X.0 for prediction, the first column
#            must be the intercept and all covariates should first be scaled
#            by the means and sds of the data used to fit the model. Code from previous exercises
#            should be helpful here. Because the variables day and obs are related
#            to the ability to detect a bird, set day = 0 in the design matrix
#            (and thus also set day^2 = 0). Set the ignore.RE argument to TRUE
#            when doing the prediction to set the random observer effect to 0.
#            Save the prediction object as out.pred.
load('data/bbs-glmm-pred-data.rda')
str(pred.df)

coords.0 <- as.matrix(pred.df[, c('x', 'y')])
# Standardize the covariates by their appropriate means and sds
tcc.pred <- (pred.df$tcc - mean(data.NOCA$covs$tcc)) / sd(data.NOCA$covs$tcc)
elev.pred <- (pred.df$tcc - mean(data.NOCA$covs$tcc)) / sd(data.NOCA$covs$tcc)
elev.2.pred <- elev.pred^2
ppt.pred <- (pred.df$ppt - mean(data.NOCA$covs$ppt)) / sd(data.NOCA$covs$ppt)
# Set day to 0
day.pred <- 0
day.2.pred <- 0
# Put together in a matrix
X.0 <- cbind(1, tcc.pred, elev.pred, elev.2.pred, ppt.pred, day.pred,
             day.2.pred)
# Note the names should match those used in the model formula
colnames(X.0) <- c('(Intercept)', 'scale(tcc)', 'scale(elev)', 'I(scale(elev)^2)',
                   'scale(ppt)', 'scale(day)', 'I(scale(day)^2)')
# Predict
out.pred <- predict(out.sp.p, X.0 = X.0, coords.0 = coords.0, ignore.RE = TRUE,
                    n.omp.threads = 4)
str(out.pred)

# Once you have completed TODO 9.10, you can run the below code to generate a
# map of Northern Cardinal Relative abundance in 2019 across North Carolina.

# Calculate the median expected relative abundance at each grid cell across the
# state.
mu.0.medians <- apply(out.pred$mu.0.samples, 2, median)

# Package the medians and coordinates together into a data frame
plot.df <- data.frame(abundance = mu.0.medians,
                      x = coords.0[, 1],
                      y = coords.0[, 2])
# Convert the data frame to a stars object for easy plotting.
pred.stars <- st_as_stars(plot.df, dims = c('x', 'y'))

# Make the map with ggplot (note we are using the map of North Carolina we used
# earlier in the exploratory data analysis).
ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = abundance), interpolate = FALSE) +
  geom_sf(data = nc, alpha = 0, col = 'black') +
  scale_fill_viridis_c(na.value = NA, option = 'plasma') +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude", fill = "", title = 'Northern Cardinal Relative Abundance 2019')
