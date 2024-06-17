# 06-plant-spatial-multi-species-occ.R: in this script, we fit a multi-species
#                                       spatial factor occupancy models to a simulated
#                                       data set of 10 plant species. Our goal is to
#                                       explore the contribution of elevation as a
#                                       niche partitioning mechanism.
rm(list = ls())
library(spOccupancy)
library(MCMCvis)
library(ggplot2)
library(dplyr)
library(tidyr)
library(corrplot)
# For plotting coordinates
library(sf)
set.seed(8736)

# Our goal in this script is to understand how elevation and other factors
# contribute to niche partitioning in a simulated community of 10 plants. We will
# do this by fitting two spatially-explicit multi-species occupancy models in which
# we simultaneously account for imperfect detection, species correlations, and
# spatial autocorrelation. The first model will not include any covariates on
# occupancy probability and will solely estimate spatial variation in occupancy
# probability based on the estimated spatial factors. The second model will be
# identical to the first model, but we will also include elevation as a predictor
# variable. We can get a sense of the importance that elevation has in niche
# partitioning across the community by looking at how the estimated residual
# correlations between species changes in the model with and without elevation.

# Data prep ---------------------------------------------------------------
# Load the data source in an object called data.list
load('data/sim-plant-elevation-data.rda')
str(data.list)

# Take a look at how many detections there are of each species
apply(data.list$y, 1, sum, na.rm = TRUE)

# Plot the coordinates in space
plot(data.list$coords, pch = 19)

# Model fitting -----------------------------------------------------------
# We first fit the model without elevation. We will use four spatial factors.
# Recall we are using a Nearest Neighbor Gaussian Process for
# each of the spatial factors, and so we also need to specify the
# number of neighbors. Here we will set n.neighbors = 5 to speed things up.
# We could use WAIC or k-fold cross-validation to compare this choice to using 15 neighbors.

# Manually specify a few priors -------
# Will use default priors for all community level parameters.
# For phi (spatial decay), we will set a prior that
# restricts the effective spatial range to fall between 0.1 and 1. This is a
# reasonably vague prior as the coordinates all fall within a 1x1km square
priors <- list(phi.unif = list(a = 3 / 1,  b = 3 / .1))

# Set MCMC criteria -------------------
n.batch <- 400
batch.length <- 25
n.burn <- 5000
n.thin <- 5
n.chains <- 3
# Run the models ----------------------
# Approx run time: <2 min
out.no.elev <- sfMsPGOcc(occ.formula = ~ 1,
                         det.formula = ~ effort.s,
                         data = data.list,
                         cov.model = 'exponential', NNGP = TRUE,
                         n.neighbors = 5, n.factors = 4, priors = priors,
                         verbose = TRUE, n.omp.threads = 1, n.report = 100,
                         n.batch = n.batch, batch.length = batch.length,
                         n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)

# Now let's fit our second model, where we include a linear and quadratic
# effect of elevation
# Approx run time: <2 min
out.elev <- sfMsPGOcc(occ.formula = ~ elev.s + elev.s.sq,
                         det.formula = ~ effort.s,
                         data = data.list,
                         cov.model = 'exponential', NNGP = TRUE,
                         n.neighbors = 5, n.factors = 4, priors = priors,
                         verbose = TRUE, n.omp.threads = 1, n.report = 100,
                         n.batch = n.batch, batch.length = batch.length,
                         n.burn = n.burn, n.thin = n.thin, n.chains = n.chains)

# Assess convergence of both models and get a quick summary. Have things converged?
# How do the ESS values look?
summary(out.no.elev)
summary(out.elev)

# Model comparison --------------------------------------------------------
# Model without elevation and 4 factors
waicOcc(out.no.elev)
# Model with elevation and 4 factors
waicOcc(out.elev)
# Including elevation results in a substantial decrease in WAIC.

# Generate estimates of residual species correlations for the two models --
# Number of posterior samples
n.post.samples <- (n.batch * batch.length - n.burn) / n.thin * n.chains
# Number of species
N <- nrow(data.list$y)
# Number of factors
q <- out.elev$q
# Array to hold posterior samples of the interspecies correlation matrices
cor.mat.samples.elev <- array(NA, dim = c(n.post.samples, N, N))
cor.mat.samples.no.elev <- array(NA, dim = c(n.post.samples, N, N))
# Form the matrix
for (l in 1:n.post.samples) {
  lambda.mat.elev <- matrix(out.elev$lambda.samples[l, ], nrow = N, ncol = q)
  cor.mat.samples.elev[l, , ] <- cov2cor(lambda.mat.elev %*% t(lambda.mat.elev))
  lambda.mat.no.elev <- matrix(out.no.elev$lambda.samples[l, ], nrow = N, ncol = q)
  cor.mat.samples.no.elev[l, , ] <- cov2cor(lambda.mat.no.elev %*% t(lambda.mat.no.elev))
}
# Get median and 95% Bayesian credible intervals for the correlation matrix
cor.mat.quants.no.elev <- apply(cor.mat.samples.no.elev, c(2, 3), quantile,
                                c(0.025, 0.5, 0.975))
cor.mat.quants.elev <- apply(cor.mat.samples.elev, c(2, 3), quantile,
                             c(0.025, 0.5, 0.975))
sp.names <- dimnames(data.list$y)[[1]]
dimnames(cor.mat.quants.no.elev)[[2]] <- sp.names
dimnames(cor.mat.quants.no.elev)[[3]] <- sp.names
dimnames(cor.mat.quants.elev)[[2]] <- sp.names
dimnames(cor.mat.quants.elev)[[3]] <- sp.names
# Median residual correlation matrix
cor.mat.quants.no.elev[2, , ]
cor.mat.quants.elev[2, , ]

par(mfrow = c(1, 2))
corrplot(cor.mat.quants.no.elev[2, , ], method = 'square', type = 'lower',
         main = 'Without Elevation')
corrplot(cor.mat.quants.elev[2, , ], method = 'square', type = 'lower',
         main = 'With Elevation')
par(mfrow = c(1, 1))
# By comparing how the values in these two maps change when including elevation
# versus not including elevation, we can get a sense of how important elevation is
# in determining the locations across the landscape each species uses. If
# all residual correlations are much closer to 0 when including elevation as opposed
# to not including elevation, that suggests elevation is an important factor in
# niche partitioning across the spatial domain. This sort of idea is discussed
# more in depth in Chapter 8 of Applied Hierarchical Modeling in Ecology version 2.

# Visualizing the spatial factors -----------------------------------------
# The spatial factors are stored in "w.samples" components of the model
# objects.
str(out.elev$w.samples)
# Extract medians
w.medians <- apply(out.elev$w.samples, c(2, 3), median)
# Convert coordinates to an sf object for easy plotting
coords.sf <- st_as_sf(data.frame(x = data.list$coords[, 1],
                                 y = data.list$coords[, 2]),
                      coords = c('x', 'y'))
# Add spatial factor medians to coords.sf object
coords.sf$factor.1 <- w.medians[1, ]
coords.sf$factor.2 <- w.medians[2, ]
coords.sf$factor.3 <- w.medians[3, ]
coords.sf$factor.4 <- w.medians[4, ]

# Make a long table for plotting in ggplot
plot.df.long <- coords.sf %>%
  pivot_longer(cols = factor.1:factor.4, names_to = 'parameter',
               values_to = 'estimate')
ggplot() +
  geom_sf(data = plot.df.long, aes(col = estimate)) +
  theme_bw() +
  scale_color_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC',
                        na.value = NA) +
  facet_wrap(vars(parameter), nrow = 1)

# Plot of species richness ------------------------------------------------
# The latent occupancy posterior samples for each species are stored in the
# z.samples component of the data list.
str(out.elev$z.samples)
# Can derive a full posterior distribution of species richness by summing the
# z.samples for each species at each site and posterior sample.
rich.samples <- apply(out.elev$z.samples, c(1, 3), sum)
# Add mean richness to the coordinates sf object for plotting
coords.sf$richness <- apply(rich.samples, 2, mean)
# Make the plot
ggplot() +
  geom_sf(data = coords.sf, aes(col = richness), size = 2) +
  theme_bw() +
  scale_color_viridis_c() +
  labs(col = 'Species\nRichness')
# Add can easily get a 95% credible interval for the derived richness estimates
# as well. This is the beauty of Bayesian analysis.
coords.sf$richness.ci <- apply(rich.samples, 2, quantile, 0.975) -
                         apply(rich.samples, 2, quantile, 0.025)
# Make the plot
ggplot() +
  geom_sf(data = coords.sf, aes(col = richness.ci), size = 2) +
  theme_bw() +
  scale_color_viridis_c() +
  labs(col = 'Species\nRichness\n95%CI')
