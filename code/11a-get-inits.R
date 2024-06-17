# 11a-get-inits.R: script to fit a multi-species spatial N-mixture model
#                  with data from the Hubbard Brook Experimental Forest.
#                  Here our goal is to predict abundance of four species across the forest.
#                  In this script, we will fit an initial multi-species spatial
#                  N-mixture model for 1 chain and a fairly small number of
#                  iterations to extract starting values for use in a larger
#                  model run.
# Clear out the workspace as usual
rm(list = ls())
library(spAbundance)
library(ggplot2)
library(sf)
library(stars)
set.seed(74823)

# Load the data set -------------------------------------------------------
# The hbefCount2015 data object in the spAbundance package contains repeated
# count data for 12 foliage-gleaning bird species across Hubbard Brook Experimental
# Forest in New Hampshire, USA.
data(hbefCount2015)
# Total individuals observed across all locations and replicates for each species
apply(hbefCount2015$y, 1, sum, na.rm = TRUE)
# We see a clear mix of rare and abundant species.

# Exploratory data analysis -----------------------------------------------
coords.sf <- st_as_sf(x = as.data.frame(hbefCount2015$coords),
                      coords = c('X', 'Y'),
                      crs = "+proj=utm +zone=19 +units=m +datum=NAD83")
# Add the total individuals observed at each site to the sf object
coords.sf$total.abund <- apply(hbefCount2015$y, 2, sum, na.rm = TRUE)
# Read in a shapefile of the HBEF watersheds to get the outline of the forest.
hbef <- st_read('data/hbef-spatial-layer/', 'hbef_wsheds')
# Join all watersheds together
hbef <- st_union(hbef)
ggplot() +
  geom_sf(data = hbef, alpha  = 0) +
  geom_sf(data = coords.sf, aes(col = total.abund)) +
  scale_color_viridis_c() +
  theme_bw() +
  labs(col = 'Total Raw\nAbundance', x = 'Longitude', y = 'Latitude') +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Subset the data ---------------------------------------------------------
# Let's work with a subset of 4 species to speed things up a bit for this
# example.
sp.names <- dimnames(hbefCount2015$y)[[1]]
my.sp <- c('OVEN', 'REVI', 'BTBW', 'BTNW')
sp.indx <- which(sp.names %in% my.sp)
data.small.hbef <- hbefCount2015
data.small.hbef$y <- hbefCount2015$y[sp.indx, , ]
str(data.small.hbef)

# Model preparation -------------------------------------------------------
# Priors ------------------------------
dist.hbef <- dist(data.small.hbef$coords)
# Note the coordinates are in meters
summary(dist.hbef)
# Note my choice of prior here. I am restricting the effective spatial
# range to occur between 500m and 7000m, which spans most of the range of the
# intersite distance values.
prior.list <- list(beta.normal = list(mean = 0, var = 10),
                   alpha.normal = list(mean = 0, var = 2.72),
                   phi.unif = list(a = 3 / 7000, b = 3 / 500))
# Initial values ----------------------
n.sp <- nrow(data.small.hbef$y)
inits.list <- list(beta = matrix(c(0, runif(2, -1, 1)), n.sp, 3,
                                 byrow = TRUE),
                   alpha = matrix(runif(4 * n.sp, -1, 1), n.sp),
                   phi = 3 / median(dist.hbef))

# Run the model -----------------------------------------------------------
# Here we will run the model using 1 latent spatial factor and will only run the model
# for a total of 10000 iterations. This model will likely take substantially
# longer to reach convergence. Here we are simply using this initial model
# run to extract initial values for all parameters, which we will then be used
# as initial values for a subsequent longer model fit.
n.factors <- 1
n.batch <- 400
batch.length <- 25
n.burn <- 5000
n.thin <- 5
n.chains <- 1
out <- sfMsNMix(abund.formula = ~ scale(elev) + I(scale(elev)^2),
                det.formula = ~ scale(tod) + scale(day) + I(scale(day)^2),
                data = data.small.hbef, inits = inits.list, priors = prior.list,
                cov.model = 'exponential', NNGP = TRUE, n.neighbors = 5,
                n.factors = n.factors, n.batch = n.batch,
                batch.length = batch.length, n.burn = n.burn,
                n.thin = n.thin, n.chains = n.chains, n.report = 50)

# Notice the very slow mixing of the abundance coefficients
plot(out, 'beta', density = FALSE)

# Extract posterior summaries for use as future initial values ------------
# Community-level abundance mean parameters
beta.comm.means <- apply(out$beta.comm.samples, 2, mean)
# Community-level abundance variance parameters
tau.sq.beta.means <- apply(out$tau.sq.beta.samples, 2, mean)
# Community-level detection mean parameters
alpha.comm.means <- apply(out$alpha.comm.samples, 2, mean)
# Community-level detection variance parameters
tau.sq.alpha.means <- apply(out$tau.sq.alpha.samples, 2, mean)
# Number of sites
J <- ncol(out$y)
# Number of abundance regression coefficients
p.occ <- ncol(out$X)
# Number of detection regression coefficients
p.det <- ncol(out$X.p)
# Species-level abundance parameters
beta.means <- matrix(apply(out$beta.samples, 2, mean), n.sp, p.occ)
# Species-level detection parameters
alpha.means <- matrix(apply(out$alpha.samples, 2, mean), n.sp, p.det)
# Spatial decay parameters
phi.means <- apply(out$theta.samples, 2, mean)
# Spatial factor loadings
lambda.means <- matrix(apply(out$lambda.samples, 2, mean), n.sp, n.factors)
# Spatial factors
w.means <- apply(out$w.samples, c(2, 3), mean)

# Package all the mean estimates into a list for the initial values
# for future model runs
inits.list <- list(beta.comm = beta.comm.means, tau.sq.beta = tau.sq.beta.means,
                   alpha.comm = alpha.comm.means, tau.sq.alpha = tau.sq.alpha.means,
                   beta = beta.means, alpha = alpha.means, phi = phi.means,
                   lambda = lambda.means, w = w.means)

# Save initial values to hard drive ---------------------------------------
save(inits.list, file = 'data/hbef-sfMsNMix-inits.rda')
