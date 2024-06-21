# 14-multi-species-hds-harv-forest.R: script to compare the three types of
#                                     multi-species HDS models in spAbundance
#                                     to estimate density of 20 bird species
#                                     across 90 point count locations in Harvard
#                                     Forest in 2022. Harvard Forest is a long
#                                     term ecological research center in the
#                                     northeastern USA.
rm(list = ls())
library(spAbundance)
# Load unmarked for the gxhn function.
library(unmarked)
# The usual packages for making plots
library(ggplot2)
library(dplyr)
library(tidyr)
library(sf)
set.seed(33278)

# Our objective here is to estimate abundance/density of a community of 10 bird
# species across 90 survey locations in a northeastern USA forest. We will
# explicitly compare the four types of multi-species HDS models in spAbundance:
# (1) a multi-species HDS model; (2) a latent factor multi-species HDS model; and
# (3) a spatial factor multi-species HDS model.

# Load the multi-species HDS data -----------------------------------------
# Loads a spAbundance list object called data.harv.
load('data/harvard-forest-HDS-data.rda')
# We see there are 20 species in the current data set, but to speed things up
# let's just work with the 10 most common species
# How many observations of each species are there?
apply(data.harv$y, 1, sum)
sp.indx <- order(apply(data.harv$y, 1, sum))[-c(1:10)]
data.list <- data.harv
data.list$y <- data.harv$y[sp.indx, , ]
# Plot the coordinates to see arrangement of the point count locations
plot(data.list$coords, pch = 19)
# Notice the plots are arranged in a set of 10 clusters, each with 9 points. Also
# notice that the coordinates are in km.

# Fit the models ----------------------------------------------------------
# For all models we will include a linear and quadratic effect of day on
# detection, as well as a linear effect of wind on detection. For abundance/density,
# we will include a categorical variable of land cover (from the National
# Land Cover Database (nlcd) in the US) as well as a continuous variable
# (percent tree canopy cover). Use a Poisson distribution for all models.

# Fit a spatial factor HDS model ------------------------------------------
# We will fit the models in reverse order of their complexity (spatial factor,
# latent factor, then regular multi-species).
# Prior distributions -----------------
# Here, as we often do, we will use an informative prior on the spatial decay
# parameter. This time, we will base the prior on the cluster sampling design.
# The 9 point counts within a given cluster are distributed 250m (.25km) apart,
# such that the maximum distance across any two points within a cluster is about .706km
# (you can calculate that distance yourself using the Pythagorean Theorem). To account
# for the sampling design, we will use a prior that forces the effective spatial
# range for each spatial factor to be at minimum .706km. In other words, our
# prior will ensure there is at least a small amount of spatial correlation
# between points within a cluster.
dists.harv <- dist(data.list$coords)
prior.list <- list(beta.comm.normal = list(mean = 0, var = 10),
                   tau.sq.beta.ig = list(a = 0.1, b = 0.1),
                   # Remember the community detection parameters in HDS are on the log
                   # scale, which is different from both N-mixture and occupancy models.
                   alpha.comm.normal = list(mean = 0, var = 10),
                   tau.sq.alpha.ig = list(a = 0.1, b = 0.1),
                   phi.unif = list(a = 3 / max(dists.harv), b = 3 / 0.706))
# Initial values ----------------------
inits.list <- list(beta.comm = 0, alpha.comm = 0, beta = 0, alpha = 0,
                   tau.sq.beta = 1, tau.sq.alpha = 1, phi = 3 / 3)
# Fit the model -----------------------
# Let's fit the model with 2 spatial factors
n.factors <- 2
# MCMC criteria
n.batch <- 800
batch.length <- 25
n.burn <- 10000
n.thin <- 10
n.chains <- 3
# Approx run time: 3.5 minutes
out.sfMsDS <- sfMsDS(abund.formula = ~ scale(tcc) + nlcd,
                     det.formula = ~ scale(day) + I(scale(day)^2) + scale(wind),
                     data = data.list, cov.model = 'exponential', NNGP = TRUE,
                     n.neighbors = 10, n.factors = n.factors, n.batch = n.batch,
                     batch.length = batch.length, family = 'Poisson',
                     priors = prior.list, inits = inits.list,
                     transect = 'point', det.func = 'halfnormal', n.thin = n.thin,
                     n.burn = n.burn, n.chains = n.chains)

# Has the model converged?
summary(out.sfMsDS)

# Fit a latent factor HDS model -------------------------------------------
# We will use the same prior distributions and initial values as before.
out.lfMsDS <- lfMsDS(abund.formula = ~ scale(tcc) + nlcd,
                     det.formula = ~ scale(day) + I(scale(day)^2) + scale(wind),
                     data = data.list, n.factors = n.factors, n.batch = n.batch,
                     batch.length = batch.length, family = 'Poisson',
                     priors = prior.list, inits = inits.list,
                     transect = 'point', det.func = 'halfnormal', n.thin = n.thin,
                     n.burn = n.burn, n.chains = n.chains)

# Has the model converged?
summary(out.lfMsDS)

# Fit a regular multi-species HDS model -----------------------------------
# We will use the same prior distributions and initial values as before.
out.msDS <- msDS(abund.formula = ~ scale(tcc) + nlcd,
                 det.formula = ~ scale(day) + I(scale(day)^2) + scale(wind),
                 data = data.list, n.batch = n.batch,
                 batch.length = batch.length, family = 'Poisson',
                 priors = prior.list, inits = inits.list,
                 transect = 'point', det.func = 'halfnormal', n.thin = n.thin,
                 n.burn = n.burn, n.chains = n.chains)

# Has the model converged
summary(out.msDS)

# Compare models with WAIC ------------------------------------------------
# Note it takes a minute or two to calculate WAIC for each of the models
waic.sfMsDS <- waicAbund(out.sfMsDS)
waic.lfMsDS <- waicAbund(out.lfMsDS)
waic.msDS <- waicAbund(out.msDS)
# Overall WAIC across all species is the sum of the individual species WAIC values.
sum(waic.sfMsDS$WAIC)
sum(waic.lfMsDS$WAIC)
sum(waic.msDS$WAIC)

# Interestingly, it seems like the spatial model performs the best, yet the
# basic multi-species model outperforms the latent factor model.

# Plot detection probability as a function of distance for all 10 species and
# the overall community ---------------------------------------------------
# Number of species we modeled.
n.sp <- nrow(data.list$y)
# Species names
sp.names <- dimnames(data.list$y)[[1]]
# Species-specific detection intercepts.
det.int.samples <- out.sfMsDS$alpha.samples[, 1:n.sp]
# Mean values of sigma at average covariate values for each species.
sigma.means <- apply(exp(det.int.samples), 2, mean)
# Mean value of sigma for the overall community
sigma.comm.means <- mean(exp(out.sfMsDS$alpha.comm.samples[, 1]))
sigma.comm.quants <- quantile(exp(out.sfMsDS$alpha.comm.samples[, 1]), c(0.025, 0.975))
# Specify a vector of distances across which you want to make the plot (this is the
# range of the x-axis of the plot). Remember that we specified the units in km and that
# each point count survey has a radius of 250km.
x.vals <- seq(0, .250, length.out = 200)
n.vals <- length(x.vals)
# Create a data frame to hold detection probability values.
det.plot.df <- data.frame(val = NA,
                          x.val = rep(x.vals, n.sp),
                          sp = rep(sp.names, each = n.vals))
# Loop through the species (probably a faster way to do this, but I'll be lazy
# here and use a for loop)
for (i in 1:n.sp) {
  indx <- ((i - 1) * n.vals + 1):(i * n.vals)
  # gxhn(x, sigma) is the half-normal detection function value at
  # distance x and scale parameter sigma.
  det.plot.df$val[indx] <- gxhn(x.vals, det.means[i])
}

comm.plot.df <- data.frame(mean = gxhn(x.vals, sigma.comm.means),
                           x.val = x.vals,
                           low = gxhn(x.vals, sigma.comm.quants[1]),
                           high = gxhn(x.vals, sigma.comm.quants[2]))
# Black is the community-level mean, gray region is the 95% credible interval
# for the community-level effect. Colored lines represent the means for each species.
ggplot(data = comm.plot.df) +
  geom_ribbon(aes(x = x.val, ymin = low, ymax = high), fill = 'grey',
              alpha = 0.5) +
  geom_line(data = det.plot.df, aes(x = x.val, y = val, col = sp), lwd = 1, lty = 1) +
  theme_bw(base_size = 14) +
  geom_line(aes(x = x.val, y = mean), col = 'black', lwd = 1.3) +
  labs(x = 'Distance (m)', y = 'Detection Probability', col = 'Species')

# Generate maps of the two latent spatial factors across the 90 sites -----
# The spatial factors are stored in "w.samples" components of the model
# objects.
str(out.sfMsDS$w.samples)
# Extract medians
w.medians <- apply(out.sfMsDS$w.samples, c(2, 3), median)
# Convert coordinates to an sf object for easy plotting
coords.sf <- st_as_sf(data.frame(x = data.list$coords[, 1],
                                 y = data.list$coords[, 2]),
                      coords = c('x', 'y'))
# Add spatial factor medians to coords.sf object
coords.sf$factor.1 <- w.medians[1, ]
coords.sf$factor.2 <- w.medians[2, ]

# Make a long table for plotting in ggplot
plot.df.long <- coords.sf %>%
  pivot_longer(cols = factor.1:factor.2, names_to = 'parameter',
               values_to = 'estimate')
ggplot() +
  geom_sf(data = plot.df.long, aes(col = estimate)) +
  theme_bw() +
  scale_color_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC',
                        na.value = NA) +
  facet_wrap(vars(parameter), nrow = 1)
