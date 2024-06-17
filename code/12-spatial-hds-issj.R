# 12-spatial-hds-issj.R: script to fit single-species spatial and non-spatial
#                        hierarchical distance sampling models to estimate the
#                        density of the island scrub jay.
# Data source citation
#    Sillett, T. S., Chandler, R. B., Royle, J. A., Kéry, M., & Morrison, S. A.
#    (2012). Hierarchical distance‐sampling models to estimate population
#    size and habitat‐specific abundance of an island endemic.
#    Ecological applications, 22(7), 1997-2006.
rm(list = ls())
library(spAbundance)
# Load the unmarked package where we will get the data and use the gxhn function
library(unmarked)
# The usual packages for making plots
library(ggplot2)
library(sf)
library(stars)
library(pals)
set.seed(232)

# Our objective here is to predict density of the island scrub jay, an
# endemic bird to the island of Santa Cruz off the western coast of the
# continental US. A distance sampling protocol was used at 307 point count surveys
# across Santa Cruz Island during the fall of 2018. The distance data were binned
# into 3 distance intervals (0-100m, 100-200m, and 200-300m). This is the full
# range of the species, as it solely occurs on this small island (250km2). Given
# the restricted range of the species and the large historical
# impact of agriculture on the island to native habitat, there is widespread
# interest in estimating it's population size.

# Load and format data for spAbundance ------------------------------------
# The "issj" data frame in unmarked contains the issj data. We will need
# to reformat this data frame into the spAbundance format.
data(issj)
str(issj)

# There are five pieces of data we will need to fit HDS models in spAbundance:
#    (1) y: the distance sampling count data, which should be formatted in a
#           two-dimensional matrix or data frame where the rows correspond to sites and
#           the columns correspond to the distance bins.
#    (2) covs: a data frame or matrix of site-level covariates for inclusion in the
#              model for abundance/density and/or detection.
#    (3) dist.breaks: a vector of distances that denotes the breakpoints of the
#                     distance bands.
#    (4) offset: an area offset that can be used to scale estimates from abundance
#                per transect/point count survey to density per unit area.
#    (5) coords: a matrix of spatial coordinates that is required for spatially-explicit
#                models.

# Distance sampling data
y <- as.matrix(issj[, 1:3])
str(y)
table(y)
# As with most count data, we see a lot of zeros, and a few larger counts.
# Covariates (elevation, % forest cover within the 300m radius circle,
#             % chaparral woodland cover within the 300m radius circle)
# More on chaparral plant communities: https://en.wikipedia.org/wiki/Chaparral
covs <- issj[, c('elevation', 'forest', 'chaparral')]
str(covs)
# Distance break points. Note these are in km.
dist.breaks <- c(0, 0.1, 0.2, 0.3)
# Offset to convert the point-level (.3km radius) abundance to abundance per hectare
# (1 ha = .01km2)
area.km2 <- pi * .3^2
area.ha.offset <- area.km2 / .01
area.ha.offset
# Thus, when we convert our estimates to density (i.e., individuals per ha), each
# 300m radius point count survey is representative of 28.27ha.
# Coordinates (divide by 1000 to convert to km)
coords <- as.matrix(issj[, c('x', 'y')]) / 1000
# Proj4string: "+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"
plot(coords, pch = 19)
# Put everything together into a list for spAbundance
data.issj <- list(y = y, covs = covs, offset = area.ha.offset,
                  dist.breaks = dist.breaks, coords = coords)
str(data.issj)
# This list is now in the required format for fitting HDS models in spAbundance.

# Model prep --------------------------------------------------------------
# Priors
dist.issj <- dist(data.issj$coords)
summary(dist.issj)
priors.list <- list(beta.normal = list(mean = 0, var = 10),
                    alpha.normal = list(mean = 0, var = 10),
                    phi.unif = c(3 / 35, 3 / 5),
                    sigma.sq.ig = c(2, 1),
                    kappa.unif = c(0, 100))
# What is the maximum and minimum effective spatial range allowed by the
# prior on phi.unif?
# We'll just use the default initial values
# Tuning variances
tuning.list <- list(beta = 0.1, alpha = 0.1, kappa = 0.2, phi = 0.5, w = 0.5)

# Run the models ----------------------------------------------------------
# Let's run two models that both account for overdispersion in the count data,
# but in different ways.
# Let's run two models: (1) Non-spatial NB; (2) Spatial Poisson.
# For both models include the following covariates:
#   Abundance/density: elevation (linear), chaparral (linear and quadratic)
#   Detection: chaparral.
# These were the covariates that Sillett et al. found to yield the best model for
# this specific portion of the data that were analyzed.
# We will use a half-normal detection function for all models.
n.batch <- 1800
batch.length <- 25
n.burn <- 20000
n.thin <- 5
n.chains <- 3
# NB non-spatial ----------------------
# Approx run time: < 2 min
out.nb <- DS(abund.formula = ~ scale(elevation) + scale(chaparral) + I(scale(chaparral)^2),
             det.formula = ~ scale(chaparral),
             data = data.issj, n.batch = n.batch, batch.length = batch.length,
             family = 'NB', det.func = 'halfnormal',
             transect = 'point', tuning = tuning.list, priors = priors.list,
             accept.rate = 0.43, n.report = 200, n.burn = n.burn, n.thin = n.thin,
             n.chains = n.chains)
summary(out.nb)
plot(out.nb, 'beta', density = FALSE)
# Spatial Poisson ---------------------
# Approx run time: 3 min
out.sp.p <- spDS(abund.formula = ~ scale(elevation) + scale(chaparral) + I(scale(chaparral)^2),
                 det.formula = ~ scale(chaparral),
                 data = data.issj, n.batch = n.batch, batch.length = batch.length,
                 family = 'Poisson', det.func = 'halfnormal',
                 transect = 'point', tuning = tuning.list, priors = priors.list,
                 accept.rate = 0.43, n.report = 200, n.burn = n.burn, n.thin = n.thin,
                 n.chains = n.chains, NNGP = TRUE, n.neighbors = 5, cov.model = 'exponential')
summary(out.sp.p)
plot(out.sp.p, 'beta', density = FALSE)

# Compare models with WAIC ------------------------------------------------
# Each of these calls takes a few seconds. Speeding this up is something I'm
# hoping to implement in the future.
waic.nb <- waicAbund(out.nb)
waic.sp.p <- waicAbund(out.sp.p)
waic.nb
waic.sp.p
# Which is the best performing model?

# Assess model goodness of fit --------------------------------------------
# Informal plot of the fitted values vs. the true values. These fitted values indicate
# the number of individuals our model predicts to observe within each distance band
# at each site.
out.sp.p.fitted <- fitted(out.sp.p)
y.rep.means <- apply(out.sp.p.fitted$y.rep.samples, c(2, 3), mean)
plot(data.issj$y, y.rep.means, pch = 19, xlab = 'True', ylab = 'Fitted')
abline(0, 1)
# Model fit doesn't look great, in particular on the left side of the plot.
# Perhaps a zero-inflated model could help (but that is not yet implemented
# in spAbundance...)
# Posterior predictive check
out.sp.p.ppc <- ppcAbund(out.sp.p, fit.stat = 'freeman-tukey', group = 1)
summary(out.sp.p.ppc)
# Note that the Bayesian p-value does not suggest any lack of fit, while our
# visual assessment of the fitted values vs. the true values suggest a limited
# ability to predict very high values.

# Interpret the best performing model -------------------------------------
summary(out.sp.p)
# What is the relationship between elevation and density?
# Remember that our coefficients are on the log scale
beta.means <- apply(out.sp.p$beta.samples, 2, mean)
beta.means
# To interpret the effect of elevation, we can calculate (exp(beta) - 1) * 100.
# This will give us the percent change in ISSJ density for a one standard
# deviation incresae in elevation.
(exp(beta.means[2]) - 1) * 100
# Thus, for a 1 standard deviation increase in elevation, we would expect ISSJ
# density/ha to decrease by 33.8 percent. What is 1 standard deviation in elevation?
sd(data.issj$covs$elevation) # About 125 m.

# Generate plot of detection probability vs. distance ------------------------
# Direct interpretation of the detection parameters can be a bit difficult given the
# way that the covariates enter into the HDS model. We see there is a negative
# effect of percent chaparral on detection, but what exactly that means is not
# super straightforward. The intercept is even more challenging to understand,
# making it difficult to assess how detection probability varies with distance.

# Let's first visualize the relationship between detection probability and distance.

# To do this, we need to generate estimates of "sigma", the scale parameter of the
# half-normal detection function. Here we will make a plot of detection vs. distance
# at mean values of chapparal (0 since we standardized it when fitting the model)
# Remember that log(sigma) = alpha_0 + alpha_1 * scale(chapparal)
sigma.samples <- exp(out.sp.p$alpha.samples[, 1] + out.sp.p$alpha.samples[, 2] * 0)
# Get median and 95% credible interval
sigma.quants <- quantile(sigma.samples, c(0.025, 0.5, 0.975))
# Specify a vector of distances across which you want to make the plot (this is the
# range of the x-axis of the plot). Remember that we specified the units in km.
x.vals <- seq(0, 0.3, length.out = 200)
# gxhn(x, sigma) is the half-normal detection function value at
# distance x and scale parameter sigma.
det.plot.df <- data.frame(med = gxhn(x.vals, sigma.quants[2]),
                          low = gxhn(x.vals, sigma.quants[1]),
                          high = gxhn(x.vals, sigma.quants[3]),
                          x.val.m = x.vals * 1000) # Notice the conversion to meters
# We see detection probability is quite high across most of the 300m radius
# point count circle at mean levels of chapparal (~27% cover).
ggplot(data = det.plot.df) +
  geom_ribbon(aes(x = x.val.m, ymin = low, ymax = high), fill = 'grey',
              alpha = 0.5) +
  theme_bw(base_size = 14) +
  geom_line(aes(x = x.val.m, y = med), col = 'black', linewidth = 1.3) +
  labs(x = 'Distance (m)', y = 'Detection Probability')

# Now let's make the same plot, but now plot the curve at three different levels
# of chapparal cover: 6%, 50%, and 85%
chapparal.vals <- c(.06, .5, .85)
# Standardize by mean and sd used to fit the model
chapparal.vals.s <- (chapparal.vals - mean(data.issj$covs$chaparral)) / sd(data.issj$covs$chaparral)
# Calculate the median sigma value for each of the three chaparral values
sigma.med.1 <- median(exp(out.sp.p$alpha.samples[, 1] +
                          out.sp.p$alpha.samples[, 2] * chapparal.vals.s[1]))
sigma.med.2 <- median(exp(out.sp.p$alpha.samples[, 1] +
                          out.sp.p$alpha.samples[, 2] * chapparal.vals.s[2]))
sigma.med.3 <- median(exp(out.sp.p$alpha.samples[, 1] +
                          out.sp.p$alpha.samples[, 2] * chapparal.vals.s[3]))
det.plot.df.2 <- data.frame(value = c(gxhn(x.vals, sigma.med.1),
                                      gxhn(x.vals, sigma.med.2),
                                      gxhn(x.vals, sigma.med.3)),
                            x.val.m = rep(x.vals * 1000, times = 3),
                            chaparral = rep(factor(c('6%', '50%', '85%'),
                                                   levels = c('6%', '50%', '85%')),
                                            each = length(x.vals)))
ggplot(data = det.plot.df.2) +
  geom_line(aes(x = x.val.m, y = value, lty = chaparral), col = 'black', linewidth = 1.1) +
  theme_bw(base_size = 14) +
  labs(x = 'Distance (m)', y = 'Detection Probability', lty = 'Percent\nChaparral')

# Prediction --------------------------------------------------------------
# Finally, let's predict density of the island scrub jay across the island.
# The prediction data are available in the "cruz" data frame, which also comes
# from the unmarked package. We will predict solely with the top model.
str(cruz)
# Note that this prediction data set represents the 2787 9ha cells that make up
# the entirety of Santa Cruz Island. Also note that the elevation in the
# prediction data set is provided in feet, while we we used values in meters
# when fitting the model, so we need to first multiply the elevation values
# by the appropriate conversion factor (0.3048).
# Format design matrix for predicting density with the standardized covariates
elev.pred <- (cruz$elevation * .3048 - mean(data.issj$covs$elevation)) / sd(data.issj$covs$elevation)
chaparral.pred <- (cruz$chaparral - mean(data.issj$covs$chaparral)) / sd(data.issj$covs$chaparral)
X.0 <- cbind(1, elev.pred, chaparral.pred, chaparral.pred^2)
colnames(X.0) <- c("(Intercept)", "scale(elevation)", "scale(chaparral)",
                   "I(scale(chaparral)^2)")
# Divide prediction coordinates in m by 1000 to get coordinates in km like we used
# to fit the model.
coords.0 <- cruz[, c('x', 'y')] / 1000
# Predict
# Approx run time: <3 min
out.pred <- predict(out.sp.p, X.0, coords.0, n.omp.threads = 4)
# Note that for distance sampling, the predicted mu.0.samples and N.0.samples
# values will be scaled according to any offset that was sent into the model
# when fitting it with spDS or DS. So, in our case, our estimates are density values
# (i.e., individuals / ha)
str(out.pred)

# Generate an abundance-based species distribution map across the species' range.
# Note because we predicted in 9ha cells and our predictions are per ha, we will
# multiply our estimates by 9.
cell.size <- 9
plot.dat <- data.frame(x = cruz$x, y = cruz$y,
                       mean.mu = apply(out.pred$mu.0.samples * cell.size, 2, mean),
                       sd.mu = apply(out.pred$mu.0.samples * cell.size, 2, mean))
plot.sf <- st_as_sf(plot.dat,
                    coords = c('x', 'y'),
                    crs = 26911)
dat.stars <- st_as_stars(plot.dat, dims = c('x', 'y'))
ggplot() +
  geom_stars(data = dat.stars, aes(x = x, y = y, fill = mean.mu)) +
  scale_fill_viridis_c("Individuals\nper 9 ha", na.value = NA) +
  theme_bw() +
  labs(x = 'Easting', y = 'Northing', title = 'Mean expected abundance') +
  coord_equal()
# There are a few cells with very high abundance that is making the plot a bit
# boring to look at. Let's discretize the plot to make it a little prettier.
dat.stars.cut <- cut(dat.stars, c(0, 0.2, 0.5, 1.0, 2.0, 5.0, ceiling(max(plot.dat$mean.mu))))
ggplot() +
  geom_stars(data = dat.stars.cut, aes(x = x, y = y, fill = mean.mu)) +
  scale_fill_viridis_d("Individuals\nper 9 ha", na.value = NA,
                       na.translate = FALSE) +
  theme_bw() +
  labs(x = 'Easting', y = 'Northing', title = 'Mean expected abundance') +
  coord_equal()
# Such discretization can often help in visualizing spatial variation in maps, but you
# should also be careful in how the discretization is done, as it has the potential to
# lead to misinformed interpretations of the map...

# Get an estimate of total population size --------------------------------
# We can get an estimate of total population size as a derived quantity from our
# predictions across the entire Santa Cruz Island. Remember that we predicted
# density on a per ha basis, but each prediction location corresponds to a 9ha
# grid cell, so we should multiple the per ha values by 9, then sum them all
# together (for each MCMC iteration) to get total abundance.
# across the island.
N.total.samples <- apply(out.pred$N.0.samples * cell.size, 1, sum)
# Median, lower, and upper bound of total population estimate.
quantile(N.total.samples, c(0.025, 0.5, 0.975))
# Note that this corresponds pretty closely with the estimate Sillet et al. (2012) reported
# using frequentist approaches (2267, 95% bootstrap CI 1613-3007).
