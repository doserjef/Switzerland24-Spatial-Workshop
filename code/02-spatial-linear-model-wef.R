# 2-spatial-linear-model-wef.R: this script introduces us to fitting spatial
#                               models and associated geostatistical
#                               techniques. We use an example data set consisting
#                               of a survey of trees in a 10 ha forest stand
#                               in Oregon, USA, with the goal of predicting
#                               diameter at breast height (DBH) across the stand.
#                               This script will also serve as our introduction
#                               to the spAbundance R package.

rm(list = ls())
library(spAbundance)
# For plots
library(ggplot2)
library(patchwork)
# For working with and plotting spatial data
library(sf)
library(stars)
# For generating interpolation plots of spatial data and residuals
library(MBA)
# For calculating semivariograms
library(geoR)

# In this example, our goal is to produce a map of DBH across a 10ha forest
# stand (that is part of the Western Experimental Forest) in Oregon, USA.
# Set seed so we all get same exact results.
set.seed(100)

# Data prep and exploratory data analysis ---------------------------------
wef.dat <- read.csv('data/WEF.csv')
str(wef.dat)

# Generate a plot showing the species IDs across the study area
# First need to make an sf object.
wef.sf <- st_as_sf(wef.dat, coords = c('East_m', 'North_m'))
ggplot() +
  geom_sf(data = wef.sf, aes(col = Species)) +
  theme_bw(base_size = 14) +
  scale_color_viridis_d() +
  labs(x = 'Easting (m)', y = 'Northing (m)')

# Prep for fitting model in spAbundance -----------------------------------
# Fitting linear models in spAbundance requires a list with two elements:
# (1) y (the response variable); and (2) covs (a data frame or matrix
# with the covariates for inclusion in the model). For spatial models,
# we also need to supply the coordinates, so we will go ahead and do that
# here for later on.
wef.data.list <- list(y = wef.dat$DBH_cm,
                      covs = data.frame(Species = wef.dat$Species),
                      coords = wef.dat[, c('East_m', 'North_m')])
# Specify priors manually (these are also the defaults). Note we need a
# prior for the regression coefficients and the residual variance for a
# simple linear model
prior.list <- list(beta.normal = list(mean = 0, var = 100),
                   tau.sq.ig = c(0.01, 0.01))
# Some basic initial values
inits.list <- list(beta = 0, tau.sq = 1)

# Fit a simple linear model -----------------------------------------------
# Note that all model-fitting functions in spAbundance (and all spatial
# model-fitting functions in spOccupancy) split up the MCMC samples into a
# set of batches (n.batch), each with a set number of
# samples in it (batch.length). The total number of MCMC samples the model
# is run for is n.batch * batch.length
out.non.sp <- abund(formula = ~ factor(Species),
                    data = wef.data.list,
                    priors = prior.list,
                    inits = inits.list,
                    n.batch = 400,
                    family = 'Gaussian',
                    batch.length = 25,
                    n.burn = 5000,
                    n.thin = 5,
                    n.chains = 3)
# Quick summary of model results
summary(out.non.sp)

# Generate surface plot of DBH and residuals ------------------------------
# Get posterior distribution of fitted values from the model
y.rep.samples <- fitted(out.non.sp)
str(y.rep.samples)
# Fitted value medians
y.hat <- apply(y.rep.samples, 2, median)
# Residuals
wef.resids <- wef.data.list$y - y.hat
# Surface plot for DBH ----------------
n.x <- 100
n.y <- 100
dbh.surface <- mba.surf(cbind(wef.data.list$coords, wef.data.list$y),
                        no.X=n.x, no.Y=n.y, h=5, m=2, extend=FALSE)$xyz.est
plot.df <- data.frame(x = rep(dbh.surface$x, times = n.y),
                      y = rep(dbh.surface$y, each = n.x),
                      dbh = c(dbh.surface$z))
pred.stars <- st_as_stars(plot.df, dims = c('x', 'y'))
plot.dbh <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = dbh),interpolate = FALSE) +
  scale_fill_gradientn(colors = rainbow(10), na.value = NA) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "Easting", y = "Northing", title = 'Interpolated DBH', fill = 'DBH')
# Surface plot for residuals ----------
dbh.surface <- mba.surf(cbind(wef.data.list$coords, wef.resids),
                        no.X=n.x, no.Y=n.y, h=5, m=2, extend=FALSE)$xyz.est
plot.df <- data.frame(x = rep(dbh.surface$x, times = n.y),
                      y = rep(dbh.surface$y, each = n.x),
                      resids = c(dbh.surface$z))
pred.stars <- st_as_stars(plot.df, dims = c('x', 'y'))
plot.resid <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = resids),interpolate = FALSE) +
  scale_fill_gradientn(colors = rainbow(10), na.value = NA) +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank()) +
  labs(x = "Easting", y = "Northing", title = 'Interpolated Residuals', fill = 'Residual')

plot.dbh + plot.resid
# Do you still see residual spatial patterns in the data?

# Visualize empirical semivariogram ---------------------------------------
# Define the distances used to calculate the empirical variogram
max.dist <- quantile(dist(wef.data.list$coords), 0.75)
bins <- 20
dbh.variog <- variog(coords = wef.data.list$coords, data = wef.data.list$y,
                     uvec = seq(0, max.dist, length = bins))
resids.variog <- variog(coords = wef.data.list$coords, data = wef.resids,
                     uvec = seq(0, max.dist, length = bins))
plot.df <- data.frame(distance = dbh.variog$u,
                      dbh = dbh.variog$v,
                      resids = resids.variog$v)
dbh.vario.plot <- ggplot(data = plot.df, aes(x = distance, y = dbh)) +
  geom_point() +
  scale_y_continuous(limits = c(0, max(plot.df$dbh))) +
  theme_bw() +
  labs(x = "Distance", y = 'Semivariance for DBH')
resids.vario.plot <- ggplot(data = plot.df, aes(x = distance, y = resids)) +
  geom_point() +
  scale_y_continuous(limits = c(0, max(plot.df$resids))) +
  theme_bw() +
  labs(x = "Distance", y = 'Semivariance for Residuals')
dbh.vario.plot + resids.vario.plot


dist.mat <- dist(wef.data.list$coords)
# Fit a spatial linear model ----------------------------------------------
# Set priors and initial values for the two spatial parameters.
# NOTE: here we change the default priors for the two spatial parameters.
#       For phi, we manually set the bounds of the effective spatial range
#       to fall within 10 and 100 m. We also run the model for a single
#       chain just to speed things up a bit. Also note this model uses an
#       approximation to a full GP. We will talk about that soon...
prior.list$phi.unif <- c(3 / 100, 3 / 10)
prior.list$sigma.sq.ig <- c(2, 50)
inits.list$sigma.sq <- 50
inits.list$phi <- 3 / 50
# Approx run time: 1.5 min
out.sp <- spAbund(formula = ~ factor(Species),
                  data = wef.data.list,
                  priors = prior.list,
                  inits = inits.list,
                  n.batch = 2000,
                  NNGP = TRUE,
                  n.neighbors = 8,
                  cov.model = 'exponential',
                  family = 'Gaussian',
                  batch.length = 25,
                  n.burn = 10000,
                  n.thin = 20,
                  n.chains = 1)
# Quick summary of model results.
summary(out.sp)

# Visualize empirical semivariogram ---------------------------------------
# Get residuals from spatial model
y.rep.samples.sp <- fitted(out.sp)
str(y.rep.samples.sp)
# Fitted value medians
y.hat.sp <- apply(y.rep.samples.sp, 2, median)
# Residuals
wef.resids.sp <- wef.data.list$y - y.hat.sp
# Define the distances used to calculate the empirical variogram
resids.variog.sp <- variog(coords = wef.data.list$coords, data = wef.resids.sp,
                           uvec = seq(0, max.dist, length = bins))
plot.df <- data.frame(distance = dbh.variog$u,
                      dbh = dbh.variog$v,
                      resids = resids.variog.sp$v)
dbh.vario.plot <- ggplot(data = plot.df, aes(x = distance, y = dbh)) +
  geom_point() +
  scale_y_continuous(limits = c(0, max(plot.df$dbh))) +
  theme_bw() +
  labs(x = "Distance", y = 'Semivariance for DBH')
resids.vario.plot <- ggplot(data = plot.df, aes(x = distance, y = resids)) +
  geom_point() +
  scale_y_continuous(limits = c(0, max(plot.df$resids))) +
  theme_bw() +
  labs(x = "Distance", y = 'Semivariance for Spatial Model Residuals')
dbh.vario.plot + resids.vario.plot

# Predict DBH for 1353 new trees ------------------------------------------
# Read in prediction data
wef.pred.dat <- read.csv("data/WEFpred.csv")
str(wef.pred.dat)
# Format covariates for prediction
# When fitting a model in spOccupancy/spAbundance with a categorical variable
# and then doing prediction, we have to supply the "design matrix" to the
# predict() function instead of just the categorical variable itself. We can
# see what this should look like by looking at the "X" component of the model
# output list. This gives the design matrix for the model fit.
str(out.sp$X)
head(out.sp$X)
# We could form this manually by hand, or we can use the model.matrix() function
# to help us out.
# By hand
X.0 <- as.matrix(data.frame(`(Intercept)` = 1,
                            `factor(Species)GF` = ifelse(wef.pred.dat$Species == 'GF',
                                                         1, 0),
                            `factor(Species)SF` = ifelse(wef.pred.dat$Species == 'SF',
                                                         1, 0),
                            `factor(Species)WH` = ifelse(wef.pred.dat$Species == 'WH',
                                                         1, 0)))
# Using model.matrix
X.0.2 <- model.matrix(~ factor(Species), wef.pred.dat)
head(X.0)
head(X.0.2)
# Form covariates
coords.0 <- as.matrix(wef.pred.dat[, c('East_m', 'North_m')])
# Predict DBH for the new trees. Can look at the manual page by typing ?predict.spAbund
pred.out <- predict(out.sp, X.0, coords.0, verbose = TRUE)
str(pred.out)

# Generate map of predicted DBH for new trees
y.0.hat <- apply(pred.out$y.0.samples, 2, mean)
pred.df <- data.frame(dbh = y.0.hat,
                      East_m = coords.0[, 1],
                      North_m = coords.0[, 2])
pred.sf <- st_as_sf(pred.df, coords = c('East_m', 'North_m'))
ggplot() +
  geom_sf(data = pred.sf, aes(col = dbh)) +
  theme_bw(base_size = 14) +
  scale_color_viridis_c() +
  labs(x = 'Easting (m)', y = 'Northing (m)', col = 'DBH')

