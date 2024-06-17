# 09b-spatial-multi-species-glmm.R: in this script we will fit a multi-species
#                                   spatial GLMM to estimate relative abundance of 88 bird
#                                   species across the state of North Carolina, USA in 2019.
#                                   Using predictions of abundance across the state
#                                   for each species, we will then generate a map
#                                   of Shannon's diversity across the state as a
#                                   derived quantity from the multi-species Bayesian
#                                   model. These data come from the North American
#                                   Breeding Bird Survey.

# Data source citation:
#    Ziolkowski Jr., D.J., Lutmerding, M., English, W.B., Aponte, V.I.,
#    and Hudson, M-A.R., 2023, North American Breeding Bird Survey
#    Dataset 1966 - 2022: U.S. Geological Survey data release, https://doi.org/10.5066/P9GS9K64.
rm(list = ls())
library(spAbundance)
library(ggplot2)
library(sf)
library(stars)
library(corrplot)
library(patchwork)

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
# Let's plot observed species richness across the data points.

# This is the coordinate system used for the coordinates
my.crs <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"
# Convert coordinates to an sf object for plotting
coords.sf <- st_as_sf(as.data.frame(data.list$coords),
                      coords = c('X', 'Y'),
                      crs = my.crs)
# Get a map of NC used for the
usa <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE))
nc <- usa %>%
  dplyr::filter(ID == 'north carolina') %>%
  st_transform(crs = my.crs)
# Add observed species richness to coordinates sf object
coords.sf$richness <- apply(ifelse(data.list$y > 0, 1, 0), 2, sum)
# Make the plot.
ggplot() +
  geom_sf(data = nc, alpha = 0, col = 'black') +
  geom_sf(data = coords.sf, aes(col = richness)) +
  theme_bw() +
  scale_color_viridis_c() +
  labs(col = 'Richness')
# It looks like species richness is highest in the northeastern portion of the state.

# Fit a spatial multi-species GLMM ----------------------------------------
# Fit a multi-species spatial GLMM using the sfMsAbund function. We will include
# the following covariates in the model: tree canopy cover (tcc, linear),
# elevation (elev, both linear and quadratic), precipitation (ppt),
# day of year (day, both linear and quadratic), and a random effect
# of observer (1 | obs). Note that the day of year and the random observer
# effect are used to account for variation in detection probability. We will use a
# Poisson distribution and use 5 neighbors in the NNGP approximation. Let's use
# 5 spatial factors to keep things (relatively) fast, although with such a large
# community of species we may want to explore using a larger number of factors.

# MCMC criteria
n.batch <- 2000
batch.length <- 25
n.burn <- 30000
n.thin <- 20
# NOTE: only fitting the model for one chain in this example, but for a real
#       example we would use multiple chains when ensuring convergence.
n.chains <- 1
# Approx run time: 20 min (you can cut back the MCMC criteria above if you want to
#                          run for a shorter period of time).
out <- sfMsAbund(formula = ~ scale(tcc) + scale(elev) + I(scale(elev)^2) +
                             scale(ppt) + scale(day) + I(scale(day)^2) +
                             (1 | obs),
                 data = data.list, family = 'Poisson', n.neighbors = 5,
                 n.factors = 5, cov.model = 'exponential', n.batch = n.batch,
                 batch.length = batch.length, n.burn = n.burn,
                 n.thin = n.thin, n.chains = n.chains, n.report = 100)
# NOTE: this model has not converged and we would likely need to run it for
#       substantially longer.
summary(out)

# Generate a residual interspecies correlation matrix ---------------------
# Number of posterior samples
n.post.samples <- (n.batch * batch.length - n.burn) / n.thin * n.chains
# Number of species
N <- nrow(data.list$y)
# Number of factors
q <- out$q
# Array to hold posterior samples of the interspecies correlation matrices
cor.mat.samples <- array(NA, dim = c(n.post.samples, N, N))
# Form the matrix
for (l in 1:n.post.samples) {
  lambda.mat <- matrix(out$lambda.samples[l, ], nrow = N, ncol = q)
  cor.mat.samples[l, , ] <- cov2cor(lambda.mat %*% t(lambda.mat))
}
# Get median and 95% Bayesian credible intervals for the correlation matrix
cor.mat.quants <- apply(cor.mat.samples, c(2, 3), quantile, c(0.025, 0.5, 0.975))
sp.names <- dimnames(data.list$y)[[1]]
dimnames(cor.mat.quants)[[2]] <- sp.names
dimnames(cor.mat.quants)[[3]] <- sp.names
# Median residual correlation matrix
cor.mat.quants[2, , ]

corrplot(cor.mat.quants[2, , ], method = 'square', type = 'lower',
         order = 'hclust')

# Predict relative abundance across the state -----------------------------
# Use the prediction data loaded below to predict abundance across the state. The prediction
# data set consists of a grid of locations across North Carolina. Note that
# when forming the design matrix X.0 for prediction, the first column
# must be the intercept and all covariates should first be scaled
# by the means and sds of the data used to fit the model. Because the
# variables day and obs are related to the ability to detect a bird,
# set day = 0 in the design matrix (and thus also set day^2 = 0).
# Set the ignore.RE argument to TRUE when doing the prediction to set the
# random observer effect to 0. After we predict abundance, we will derive two
# biodiversity metrics: species richness and Shannon's diversity.
# Load the data
load('data/bbs-glmm-pred-data.rda')
str(pred.df)

# Extract coordinates for prediction
coords.0 <- as.matrix(pred.df[, c('x', 'y')])
# Standardize the covariates by their appropriate means and sds
tcc.pred <- (pred.df$tcc - mean(data.list$covs$tcc)) / sd(data.list$covs$tcc)
elev.pred <- (pred.df$tcc - mean(data.list$covs$tcc)) / sd(data.list$covs$tcc)
elev.2.pred <- elev.pred^2
ppt.pred <- (pred.df$ppt - mean(data.list$covs$ppt)) / sd(data.list$covs$ppt)
# Set day to 0
day.pred <- 0
day.2.pred <- 0
# Put together in a matrix
X.0 <- cbind(1, tcc.pred, elev.pred, elev.2.pred, ppt.pred, day.pred,
             day.2.pred)
# Note the names should match those used in the model formula
colnames(X.0) <- c('(Intercept)', 'scale(tcc)', 'scale(elev)', 'I(scale(elev)^2)',
                   'scale(ppt)', 'scale(day)', 'I(scale(day)^2)')
# Predict. Note the use of n.omp.threads here can really speed things up.
out.pred <- predict(out, X.0 = X.0, coords.0 = coords.0, ignore.RE = TRUE,
                    n.omp.threads = 4)
str(out.pred)

# Derived biodiversity metrics --------------------------------------------
# Species richness --------------------
# Posterior distribution of richness. Species richness is calculated using the
# predicted relative abundance values (y.0.samples). We count each species in the
# richness calculation that has a relative abundance > 0.
rich.samples <- apply(ifelse(out.pred$y.0.samples > 0, 1, 0), c(1, 3), sum)
str(rich.samples)
hist(rich.samples)
# Shannon's diversity -----------------
# Recall the formula for Shannon's diversity is -sum(p_i * log(p_i)), where p_i
# is the proportion of each species (in terms of abundance) across the entire community.

# First calculate the total abundance across all species at each site for each
# MCMC iteration
total.abund.samples <- apply(out.pred$y.0.samples, c(1, 3), sum)
# Now calculate the proportion of the total abundance samples for each species
prop.abund.samples <- array(NA, dim = c(nrow(out.pred$y.0.samples),
                                        ncol(out.pred$y.0.samples),
                                        dim(out.pred$y.0.samples)[3]))
for (l in 1:nrow(out.pred$y.0.samples)) {
  for (i in 1:N) {
    prop.abund.samples[l, i, ] <- out.pred$y.0.samples[l, i, , 1] / total.abund.samples[l, ]
  }
}
# Calculate Shannon's diversity index for each MCMC sample
shannon.div.samples <- matrix(NA, nrow(out.pred$y.0.samples), nrow(pred.df))
for (l in 1:nrow(out.pred$y.0.samples)) {
  for (j in 1:nrow(pred.df)) {
  curr.indx <- which(prop.abund.samples[l, , j] > 0)
  shannon.div.samples[l, j] <- -1 * sum(prop.abund.samples[l, curr.indx, j] *
                                        log(prop.abund.samples[l, curr.indx, j]))
  }
}
# Make plots of means and 95% Credible Intervals
richness.meds <- apply(rich.samples, 2, median)
richness.ci.width <- apply(rich.samples, 2, quantile, 0.975) -
                     apply(rich.samples, 2, quantile, 0.025)
shannon.meds <- apply(shannon.div.samples, 2, median)
shannon.ci.width <- apply(shannon.div.samples, 2, quantile, 0.975) -
                    apply(shannon.div.samples, 2, quantile, 0.025)
plot.df <- data.frame(rich.med = richness.meds,
                      rich.ci = richness.ci.width,
                      shannon.med = shannon.meds,
                      shannon.ci = shannon.ci.width,
                      x = coords.0[, 1],
                      y = coords.0[, 2])
# Convert the data frame to a stars object for easy plotting.
pred.stars <- st_as_stars(plot.df, dims = c('x', 'y'))

# Make the maps with ggplot (note we are using the map of North Carolina we used
# earlier in the exploratory data analysis).
rich.med.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = rich.med), interpolate = FALSE) +
  geom_sf(data = nc, alpha = 0, col = 'black') +
  scale_fill_viridis_c(na.value = NA, option = 'plasma') +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude", fill = "", title = '(a) Richness medians')
rich.ci.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = rich.ci), interpolate = FALSE) +
  geom_sf(data = nc, alpha = 0, col = 'black') +
  scale_fill_viridis_c(na.value = NA, option = 'plasma') +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude", fill = "", title = '(b) Richness 95% CI width')
shannon.med.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = shannon.med), interpolate = FALSE) +
  geom_sf(data = nc, alpha = 0, col = 'black') +
  scale_fill_viridis_c(na.value = NA, option = 'viridis') +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude", fill = "", title = '(c) Shannon Diversity medians')
shannon.ci.plot <- ggplot() +
  geom_stars(data = pred.stars, aes(x = x, y = y, fill = shannon.ci), interpolate = FALSE) +
  geom_sf(data = nc, alpha = 0, col = 'black') +
  scale_fill_viridis_c(na.value = NA, option = 'viridis') +
  theme_bw() +
  labs(x = "Longitude", y = "Latitude", fill = "", title = '(d) Shannon diversity 95% CI width')

rich.med.plot + rich.ci.plot + shannon.med.plot + shannon.ci.plot
