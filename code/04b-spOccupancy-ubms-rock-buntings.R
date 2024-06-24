

# Code for 'A tutorial on the fitting of spatial occupancy models
# to Swiss Bird Atlas data in SE Switzerland with spOccupancy and ubms'

# MK, February-June 2024


# 1 Introduction
# --------------
# No code


# 2 Preparation and summary of the Rock bunting and habitat data
# --------------------------------------------------------------

# Note that we start working with the R workspace named
# 'Data_RockBunting_SE_Switzerland.Rdata'


# Look at contents of the R workspace
ls()

str(AtlasData)

head(AtlasData$occ.covs, 3)

str(AtlasData$det.covs)

str(CovarData)
head(CovarData)    # not shown

# Compute extents of the prediction domain
diff(range(CovarData$x))
diff(range(CovarData$y))

# Load needed packages
require(raster)     # Later have to use something else, e.g., sf
mapPalette1 <- colorRampPalette(c("grey", "yellow", "orange", "red"))

# Buildings, elevation, farmland and forest
range(CovarData[, "buildings"])
range(CovarData[, "elevation"])
range(CovarData[, "farmland"])
range(CovarData[, "forest"])

par(mfrow = c(2, 2), mar = c(2,2,4,6), cex.main = 1.2)
r1 <- rasterFromXYZ(data.frame(x = CovarData[,1], y = CovarData[,2], 
  z = CovarData[, "buildings"]))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Buildings")
r2 <- rasterFromXYZ(data.frame(x = CovarData[,1], y = CovarData[,2], 
  z = CovarData[, "elevation"]))
plot(r2, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Elevation")
r3 <- rasterFromXYZ(data.frame(x = CovarData[,1], y = CovarData[,2], 
  z = CovarData[, "farmland"]))
plot(r3, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Farmland")
r4 <- rasterFromXYZ(data.frame(x = CovarData[,1], y = CovarData[,2], 
  z = CovarData[, "forest"]))
plot(r4, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Forest")

# Glacier, grassland, lakes and nitrogen
range(CovarData[, "glacier"])
range(CovarData[, "grassland"])
range(CovarData[, "lakes"])
range(CovarData[, "nitrogen"])

par(mfrow = c(2, 2), mar = c(2,2,4,6), cex.main = 1.2)
r1 <- rasterFromXYZ(data.frame(x = CovarData[,1], y = CovarData[,2], 
  z = CovarData[, "glacier"]))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Glacier")
r2 <- rasterFromXYZ(data.frame(x = CovarData[,1], y = CovarData[,2], 
  z = CovarData[, "grassland"]))
plot(r2, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Grassland")
r3 <- rasterFromXYZ(data.frame(x = CovarData[,1], y = CovarData[,2], 
  z = CovarData[, "lakes"]))
plot(r3, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Lakes")
r4 <- rasterFromXYZ(data.frame(x = CovarData[,1], y = CovarData[,2], 
  z = CovarData[, "nitrogen"]))
plot(r4, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Nitrogen")

# Northness, Rivers, Roads and Rocks
range(CovarData[, "northness"])
range(CovarData[, "rivers"])
range(CovarData[, "roads"])
range(CovarData[, "rocks"])

par(mfrow = c(2, 2), mar = c(2,2,4,6), cex.main = 1.2)
r1 <- rasterFromXYZ(data.frame(x = CovarData[,1], y = CovarData[,2], 
  z = CovarData[, "northness"]))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Northness")
r2 <- rasterFromXYZ(data.frame(x = CovarData[,1], y = CovarData[,2], 
  z = CovarData[, "rivers"]))
plot(r2, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Rivers")
r3 <- rasterFromXYZ(data.frame(x = CovarData[,1], y = CovarData[,2], 
  z = CovarData[, "roads"]))
plot(r3, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Roads")
r4 <- rasterFromXYZ(data.frame(x = CovarData[,1], y = CovarData[,2], 
  z = CovarData[, "rocks"]))
plot(r4, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Rocks")

# Shoreline, Slope, Structures and wetlands
range(CovarData[, "shoreline"])
range(CovarData[, "slope"])
range(CovarData[, "structures"])
range(CovarData[, "wetlands"])

r1 <- rasterFromXYZ(data.frame(x = CovarData[,1], y = CovarData[,2], 
  z = CovarData[, "shoreline"]))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Shoreline")
r2 <- rasterFromXYZ(data.frame(x = CovarData[,1], y = CovarData[,2], 
  z = CovarData[, "slope"]))
plot(r2, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Slope")
r3 <- rasterFromXYZ(data.frame(x = CovarData[,1], y = CovarData[,2], 
  z = CovarData[, "structures"]))
plot(r3, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Structures")
r4 <- rasterFromXYZ(data.frame(x = CovarData[,1], y = CovarData[,2], 
  z = CovarData[, "wetlands"]))
plot(r4, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Wetlands")

# kfrivers (no idea what that is...)
range(CovarData[, "kfrivers"])

rX <- rasterFromXYZ(data.frame(x = CovarData[,1], y = CovarData[,2], 
  z = CovarData[, "kfrivers"]))
plot(rX, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "KFRivers")

# Compute the observed presence/absence state (zobs)
zobs <- apply(AtlasData$y, 1, max, na.rm = TRUE)
table(zobs)

# Observed SDM in SE Switzerland
plot(AtlasData$coords[,1], AtlasData$coords[,2], pch = 16, asp = 1, main = paste('Locations of data and of detections (red) of \n Rock Buntings in SE Switzerland'), cex = 0.5, col = 'grey', xlab = 'x coordinate', ylab = 'y coordinate')
points(AtlasData$coords[,1][zobs == 1], AtlasData$coords[,2] [zobs == 1], pch = 16, cex=0.6, col = 'red')


# 3 Fitting nonspatial models in unmarked, ubms and spOccupancy
# -------------------------------------------------------------

library(unmarked)
library(ubms)
library(spOccupancy)

# Calculate proportion missing responses
mean(is.na(AtlasData$y))


# Start with unmarked
# Prepare environmental variates at site-level
covs <- as.data.frame(AtlasData$occ.covs)  # Grab all covariates

# Prepare detection variates
detcovs <- list(date1 = AtlasData$det.covs$date1, date2 = AtlasData$det.covs$date2)

# Package all the data in an unmarked data frame for occu()
umf <- unmarkedFrameOccu(
  y = AtlasData$y,    # Pres/Abs measurements
  siteCovs = covs,    # Environmental covariates at site-level
  obsCovs = detcovs)  # Observation-specific covariates
summary(umf)

# Fitting a large model with effects of 9 covariates (ART = 40 sec)
system.time(
  summary(
  fm1 <- occu(
   ~ date1 + date2
   ~ z.buildings + z.buildings2 + z.elevation + z.elevation2 + 
     z.northness + z.northness2 + z.rivers + z.rivers2 + 
     z.rocks + z.rocks2 + z.slope + z.slope2 + 
     z.structures + z.structures2 + z.wetlands + z.wetlands2 +
     z.kfrivers + z.kfrivers2,
     se = TRUE,   # Could set to FALSE for exploration
     data=umf, control=list(trace=T, REPORT=5, maxit = 500))) )

# We continue with ubms, which is a wrapper for Stan. 
# It takes a hellufalot of time to produce the results…

# Fitting the model in Stan using ubms (ART = 51 min)
system.time(
  fm2 <- stan_occu(
   ~ date1 + date2
   ~ z.buildings + z.buildings2 + z.elevation + z.elevation2 + 
     z.northness + z.northness2 + z.rivers + z.rivers2 + 
     z.rocks + z.rocks2 + z.slope + z.slope2 + 
     z.structures + z.structures2 + z.wetlands + z.wetlands2 +
     z.kfrivers + z.kfrivers2,
     data=umf) )
ubms::traceplot(fm2)       # Check convergence of chains
print(fm2)                 # Print posterior summaries

# Finally, we fit the same model with spOccupancy.
str(AtlasData)

# Specify model formulas
# Occupancy
occ.formula <- 
  ~ z.buildings + z.buildings2 + z.elevation + z.elevation2 + 
    z.northness + z.northness2 + z.rivers + z.rivers2 + 
    z.rocks + z.rocks2 + z.slope + z.slope2 + 
    z.structures + z.structures2 + z.wetlands + z.wetlands2 +
    z.kfrivers + z.kfrivers2

# Next would be for all 17 * 2 terms …
full.occ.formula <- ~ z.buildings + z.buildings2 + 
  z.elevation + z.elevation2 + z.farmland + z.farmland2 + 
  z.forest + z.forest2 + z.glacier + z.glacier2 + 
  z.grassland + z.grassland2 + z.lakes + z.lakes2 + 
  z.nitrogen + z.nitrogen2 + z.northness + z.northness2 + 
  z.rivers + z.rivers2 + z.roads + z.roads2 + z.rocks + z.rocks2 +
  z.shoreline + z.shoreline2 + z.slope + z.slope2 + 
  z.structures + z.structures2 + z.wetlands + z.wetlands2 + 
  z.kfrivers + z.kfrivers2

# Detection
det.formula <- ~ date1 + date2

# Fit the nonspatial occupancy model
fm3 <- PGOcc(occ.formula = occ.formula,
             det.formula = det.formula,
             data = AtlasData,
             n.samples = 11000,
             n.omp.threads = 6,
             n.burn = 1000,
             n.thin = 20,
             n.chains = 4,
             n.report = 1000, verbose = TRUE)
summary(fm3)

Almost everything has converged fine according to Rhat. We check out the traceplots as well.

# Traceplots
plot(fm3$beta.samples)       # Occupancy params (not shown)
plot(fm3$alpha.samples)      # Detection params

# We compare the estimates from the three packages. 
# For this, we produce the posterior means as our Bayesian point estimators 
# for the solutions from spOccupancy (i.e., for fm3).

# Get posterior means from spOccupancy output
coef_fm3 <- c(apply(fm3$beta.samples, 2, mean), apply(fm3$alpha.samples, 2, mean))

# We compare them in a plot.

off <- 0.2
par(mar = c(6,6,6,3), cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5)
plot((1:22)-off, coef(fm1), pch = 16, cex = 1.5, col = 'blue', xlab = 'Parameter Number', ylab = 'Estimate', frame = FALSE, xlim = c(0, 25))
points((1:22), coef(fm2), pch = 16, cex = 1.5, col = 'brown')
points((1:22)+off, coef_fm3, pch = 16, cex = 1.5, col = 'green')
abline(h = 0, lwd = 1, lty = 3, col = 'grey')
legend('bottomright', pch = 16, cex = 1.2, col = c('blue', 'brown', 'green'), legend = c("unmarked (MLE)", "ubms/Stan (Bayes)", "spOccupancy (Bayes)"))  # bty = 'n'



# 3.3 Spatial predictions from the nonspatial model: first species distribution maps
# ----------------------------------------------------------------------------------

# Make a data frame containing the covariates for the prediction data set
cbind(names(CovarData))     # not shown
newdata <- data.frame(CovarData[20:53])

pred.um <- unmarked::predict(fm1, newdata = newdata, type = 'state')
options(scipen = 10)
head(pred.um)

occ.formula <- 
  ~ z.buildings + z.buildings2 + z.elevation + z.elevation2 + 
    z.northness + z.northness2 + z.rivers + z.rivers2 + 
    z.rocks + z.rocks2 + z.slope + z.slope2 + 
    z.structures + z.structures2 + z.wetlands + z.wetlands2 +
    z.kfrivers + z.kfrivers2

# Same by hand for the point prediction
selected.covs <- c(1:4, 17:20, 23:24, 27:34)
link.scale.pred <- as.matrix(cbind(1, newdata[,selected.covs])) %*% unmarked::coef(fm1)[1:19]
pred.by.hand <- plogis(link.scale.pred)    # Apply inverse link

# Test that we get the same as does unmarked with predict()
plot(pred.um[,1], pred.by.hand, xlab = 'by unmarked::predict()', ylab = 'by hand', pch = 16, cex = 0.6, col = rgb(0,0,0,0.2), frame = FALSE)


# Predictions with ubms
pred.ubms <- ubms::predict(fm2, newdata = newdata, submodel = 'state')
str(pred.ubms)
head(pred.ubms)


# With spOccupancy
# Form predictions from the non-spatial model (takes 2 sec)
# (using the spOccupancy fitted model object)
spo.cov <- as.matrix(cbind(1, newdata[,selected.covs]))
system.time(
  pred.spo <- predict(fm3, spo.cov)
)
str(pred.spo)

# Get posterior means and SDs for the model fit from spOccupancy
pm.psi <- apply(pred.spo$psi.0.samples, 2, mean)
psd.psi <- apply(pred.spo$psi.0.samples, 2, sd)

summary(pm.psi)
summary(psd.psi)

# We can now plot these predictions in geographic space to obtain our species distribution maps.

# Load needed packages
require(raster)
mapPalette1 <- colorRampPalette(c("grey", "yellow", "orange", "red"))

# Bunting species distribution map in SE Switzerland from unmarked:

par(mfrow = c(1, 2), mar = c(1,3,4,6), cex.main = 1.5)
# Point predictions (based on MLEs of the parameters)
r1 <- rasterFromXYZ(data.frame(x = CovarData$x, y = CovarData$y, 
  z = pred.um[,1]))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Rock bunting in SE Switzerland:\nOccupancy probability (unmarked)", zlim = c(0, 1))

# Prediction standard errors
r1 <- rasterFromXYZ(data.frame(x = CovarData$x, y = CovarData$y, 
  z = pred.um[,2]))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Occupancy uncertainty\n(prediction SE from unmarked)", zlim = c(0, 0.3))

# And the Rock Bunting species distribution map from ubms:

par(mfrow = c(1, 2), mar = c(1,3,4,6), cex.main = 1.5)
# Point predictions (posterior means)
r1 <- rasterFromXYZ(data.frame(x = CovarData$x, y = CovarData$y, 
  z = pred.ubms[,1]))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Rock bunting in SE Switzerland:\nOccupancy probability (ubms)", zlim = c(0, 1))

# Prediction uncertainty (posterior standard deviations)
r1 <- rasterFromXYZ(data.frame(x = CovarData$x, y = CovarData$y, 
  z = pred.ubms[,2]))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Occupancy uncertainty\n(posterior SDs from ubms)", zlim = c(0, 0.3))

# And, finally, the Rock Bunting species distribution map from spOccupancy:

par(mfrow = c(1, 2), mar = c(1,3,4,6), cex.main = 1.5)
# Point predictions (posterior means)
r1 <- rasterFromXYZ(data.frame(x = CovarData$x, y = CovarData$y, 
  z = pm.psi))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Rock bunting in SE Switzerland:\nOccupancy probability (spOccupancy)", zlim = c(0, 1))

# Prediction uncertainty (posterior standard deviations)
r1 <- rasterFromXYZ(data.frame(x = CovarData$x, y = CovarData$y, 
  z = psd.psi))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Occupancy uncertainty\n(posterior SDs from spOccupancy)", zlim = c(0, 0.3))

# They should all look almost the same.
#  To test this, we subtract the prediction from one of the other and plot these differences quickly.
# Compute differences and plot them

diff_um_ubms <- pred.um[,1] - pred.ubms[,1]
diff_ubms_spo <- pred.ubms[,1] - pm.psi
range(diff_um_ubms)
range(diff_ubms_spo)

par(mfrow = c(2, 2), mar = c(1,3,4,8), cex.main = 1.5)
# Difference in point predictions between unmarked and ubms
r1 <- rasterFromXYZ(data.frame(x = CovarData$x, y = CovarData$y, 
  z = diff_um_ubms))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Difference: unmarked minus ubms", zlim = c(-0.1, 0.05))

# Difference in point predictions between ubms and spOccupancy
r1 <- rasterFromXYZ(data.frame(x = CovarData$x, y = CovarData$y, 
  z = diff_ubms_spo))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Difference: ubms minus spOccupancy", zlim = c(-0.005, 0.02))

par(mar = c(6,6,4,4))

hist(diff_um_ubms, xlab = "unmarked minus ubms", col = 'grey', breaks = 500, xlim = c(-0.01, 0.01))
abline(v = 0, col = 'red', lwd = 3, lty = 3)

hist(diff_ubms_spo, xlab = "ubms minus spOccupancy", col = 'grey', breaks = 500, xlim = c(-0.005, 0.005))
abline(v = 0, col = 'red', lwd = 3, lty = 3)

# Compute range size

# Compute range size from unmarked (would have to bootstrap for SE)
rangeSize.um <- sum(pred.um[,1])

# Compute range size from ubms FOR NOW
rangeSize.ubms <- sum(pred.ubms[,1])

# Compute posterior distribution of range size from spOccupancy
rangeSize.spo <- apply(pred.spo[[1]], 1, sum)
hist(rangeSize.spo, breaks = 30, col = 'grey', main = 'Estimated range size of the Rock Bunting in SE Switzerland\n(Posterior distribution from non-spatial model in spOccupancy)\n (point estimates from unmarked (red) and ubms (blue) also shown', xlab = 'Number of sq.kms', freq = FALSE)
abline(v = rangeSize.um, lwd = 5, lty = 1, col = rgb(1,0,0,0.4))
abline(v = rangeSize.ubms, lwd = 5, lty = 1, col = rgb(0,0,1,0.4))




# 4 Fitting spatial models to the SE Swiss Rock Bunting data with spOccupancy
# ---------------------------------------------------------------------------


# 4.1 Spatial exponential model with 5-Nearest-neighbour Gaussian Process (NNGP 5)
# -------------------------------------------------------------------------------

# Compute distances between sites
distMat <- dist(AtlasData$coords)

# Select exponential covariance model
cov.model <- "exponential"

# Choose inits
# as in the model-fitting package vignette
spo.inits <- list(alpha = 0, 
                 beta = 0, 
                 z = apply(AtlasData$y, 1, max, na.rm = TRUE), 
                 sigma.sq = 2, 
                 phi = 3 / mean(distMat), 
                 w = rep(0, nrow(AtlasData$y)))

# Set some tuning param
spo.tuning <- list(phi = 1)

# Define some priors (as in vignette)
min.dist <- min(distMat)
max.dist <- max(distMat)
spo.priors <- list(beta.normal = list(mean = 0, var = 2.72), 
                   alpha.normal = list(mean = 0, var = 2.72),
                   sigma.sq.ig = c(2, 1), 
                   phi.unif = c(3/max.dist, 3/min.dist))

# We launch the model a first time, with just short chains.

# Resultant number of draws for each chain
4 * ((50 * 100) - 2000) / 12

# Run a quickie to get ballpark estimates (ART 3.4 mins)
fm4a <- 
    spPGOcc(occ.formula = occ.formula, 
       det.formula = det.formula, 
       data = AtlasData, inits = spo.inits, 
       n.batch = 100, batch.length = 50, 
       priors = spo.priors, cov.model = cov.model, 
       NNGP = TRUE, n.neighbors = 5, tuning = spo.tuning, 
       n.omp.threads = 6,
       n.burn = 2000, n.thin = 12, n.chains = 4, 
       n.report = 100, verbose = TRUE)
summary(fm4a)

plot(fm4a$alpha.samples)   # Traceplots of detection params, not shown
plot(fm4a$beta.samples)    # Traceplots of occupancy params, not shown
plot(fm4a$theta.samples)   # Traceplots of covariance params

# Longer chains

# Resultant number of draws for each chain
4 * ((500 * 200) - 50000) / 200

# Choose inits, parly on solution from short run
# as in the model-fitting package vignette
spo.inits <- list(alpha = 0, 
                 beta = 0, 
                 z = apply(AtlasData$y, 1, max, na.rm = TRUE), 
                 sigma.sq = 7.0934, 
                 phi = 0.0452, 
                 w = rep(0, nrow(AtlasData$y)))

# Launch model in spOccupancy
fm4b <- 
    spPGOcc(occ.formula = occ.formula, 
       det.formula = det.formula, 
       data = AtlasData, inits = spo.inits, 
       n.batch = 500, batch.length = 200, 
       priors = spo.priors, cov.model = cov.model, 
       NNGP = TRUE, n.neighbors = 5, tuning = spo.tuning, 
       n.omp.threads = 6,
       n.burn = 50000, n.thin = 200, n.chains = 4, 
       n.report = 100, verbose = TRUE)
summary(fm4b)

plot(fm4b$alpha.samples)   # Traceplots of detection params, not shown
plot(fm4b$beta.samples)    # Traceplots of occupancy params, not shown
plot(fm4b$theta.samples)   # Traceplots of covariance params

# This looks much better now.

# GoF

system.time(              # 56 sec
  ppc.out <- ppcOcc(fm4b, fit.stat = 'freeman-tukey', group = 1)
  )
summary(ppc.out)

ppc.df <- data.frame(fit = ppc.out$fit.y, 
                     fit.rep = ppc.out$fit.y.rep, 
                     color = 'lightskyblue1')
ppc.df$color[ppc.df$fit.rep > ppc.df$fit] <- 'lightsalmon'
xylim <- range(c(ppc.df$fit, ppc.df$fit.rep))
plot(ppc.df$fit, ppc.df$fit.rep, bg = ppc.df$color, pch = 21, 
     ylab = 'Fit', xlab = 'True', xlim = xylim, ylim = xylim)
abline(0, 1)
#lines(ppc.df$fit, ppc.df$fit, col = 'black')

# Model does not fit. Quantify magnitude of lack of fit

# Compute informal lack of fit ratio
hist(ppc.df$fit/ppc.df$fit.rep, xlim = c(0.9, 1.7), main = '"Lack of fit" ratio')
abline(v = 1, col = 'red', lwd = 3)
abline(v = mean(ppc.df$fit/ppc.df$fit.rep), col = 'blue', lwd = 3)


# Try to see where the model breaks

# Compute a residual-like quantity from the PPC results
resi <- ppc.out$fit.y.group.quants[3, ] - ppc.out$fit.y.rep.group.quants[3, ]
plot(resi, pch = 19, cex = 1, col = rgb(0,0,0,0.3), 
xlab = 'Site ID', ylab = 'Observed - Replicate Discrepancy', main = 'Serial plot of PPC residuals (Observed - Expected point-wise discrepancy)')

# Make a map of the posterior means of predicted occupancy 
# with 'residuals' plotted over it
par(mfrow = c(1, 1), mar = c(4,4,6,8), cex.main = 1.5)
r1 <- rasterFromXYZ(data.frame(x = CovarData[1], y = CovarData[2], 
  z = pm.psi4b))

plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Rock bunting distribution and 'PPC residuals'\n(green: OK residuals, black: high residuals)", zlim = c(0, 1))

points(x = AtlasData$coords[,1][abs(resi) < 0.5],
   y = AtlasData$coords[,2][abs(resi) < 0.5],
   cex = 0.5, pch = 16, col = 'green')
points(x = AtlasData$coords[,1][resi > 0.5],
   y = AtlasData$coords[,2][resi > 0.5],
   cex = 0.5, pch = 16, col = 'black')
points(x = AtlasData$coords[,1][resi < -0.5],
   y = AtlasData$coords[,2][ resi < -0.5],
   cex = 0.5, pch = 16, col = 'purple')

# Plot 'PPC residuals' against the predictions directly
quad1 <- paste(AtlasData$coords[,1], AtlasData$coords[,2], sep = '.')
quad2 <- paste(CovarData$x, CovarData$y, sep = '.')
idx <- pmatch(quad1, quad2)
plot(pm.psi4b[idx], resi, xlab = 'Prediction of occupancy prob.', ylab = 'PPC residual', frame = FALSE)

# Next, we plot the PPC residuals against all covariates.

# Plot 'PPC residuals' against the covariates
par(mfrow = c(2, 2), mar = c(4,4,6,2), cex.main = 1.5)
plot(CovarData$buildings[idx], resi, xlab = '', ylab = 'PPC residual', main = 'Buildings', frame = FALSE)
plot(CovarData$elevation[idx], resi, xlab = '', ylab = 'PPC residual', main = 'Elevation', frame = FALSE)
plot(CovarData$farmland[idx], resi, xlab = '', ylab = 'PPC residual', main = 'Farmland', frame = FALSE)
plot(CovarData$forest[idx], resi, xlab = '', ylab = 'PPC residual', main = 'Forest', frame = FALSE)

plot(CovarData$glacier[idx], resi, xlab = '', ylab = 'PPC residual', main = 'Glacier', frame = FALSE)
plot(CovarData$grassland[idx], resi, xlab = '', ylab = 'PPC residual', main = 'Grassland', frame = FALSE)
plot(CovarData$lakes[idx], resi, xlab = '', ylab = 'PPC residual', main = 'Lakes', frame = FALSE)
plot(CovarData$slope[idx], resi, xlab = '', ylab = 'PPC residual', main = 'Slope', frame = FALSE)

plot(CovarData$northness[idx], resi, xlab = '', ylab = 'PPC residual', main = 'Northness', frame = FALSE)
plot(CovarData$rivers[idx], resi, xlab = '', ylab = 'PPC residual', main = 'Rivers', frame = FALSE)
plot(CovarData$roads[idx], resi, xlab = '', ylab = 'PPC residual', main = 'Roads', frame = FALSE)
plot(CovarData$rocks[idx], resi, xlab = '', ylab = 'PPC residual', main = 'Rocks', frame = FALSE)

plot(CovarData$shoreline[idx], resi, xlab = '', ylab = 'PPC residual', main = 'Shoreline', frame = FALSE)
plot(CovarData$slope[idx], resi, xlab = '', ylab = 'PPC residual', main = 'Slope', frame = FALSE)
plot(CovarData$structures[idx], resi, xlab = '', ylab = 'PPC residual', main = 'Structures', frame = FALSE)
plot(CovarData$wetlands[idx], resi, xlab = '', ylab = 'PPC residual', main = 'Wetlands', frame = FALSE)


# Model selection by WAIC

(waic1 <- waicOcc(fm3))           # nonspatial model
(waic2b <- waicOcc(fm4b))          # spatial model

# So, the spatial model is preferred by a large margin over the corresponding spatial model.
# We go on forming and the plotting predictions now.
#  This proceeds very similarly as what we did above for the non-spatial model, 
# but now we must also supply the coordinates for the prediction function.

spo.cov <- as.matrix(cbind(1, newdata[,selected.covs]))
str(spo.cov)

spo.coords <- CovarData[1:2]
head(spo.coords)

# Form predictions from the spatial model (NNGP with 5 neighbours)
system.time(
  pred.fm4b <- predict(fm4b, spo.cov, spo.coords, verbose = TRUE)
)

# Make species distribution map

# Load needed packages
require(raster)

# Map posterior mean and posterior SD of occupancy probability
# We also compute posterior means of w
pm.psi4b <- apply(pred.fm4b$psi.0.samples, 2, mean)  # Post mean
psd.psi4b <- apply(pred.fm4b$psi.0.samples, 2, sd)   # Post sd
pm.w4b <- apply(pred.fm4b$w.0.samples, 2, mean)      # Post mean process

summary(pm.psi4b)
summary(psd.psi4b)
summary(pm.w4b)


mapPalette1 <- colorRampPalette(c("grey", "yellow", "orange", "red"))

par(mfrow = c(1, 3), mar = c(4,4,6,8), cex.main = 2)
# Posterior means
r1 <- rasterFromXYZ(data.frame(x = CovarData[1], y = CovarData[2], 
  z = pm.psi4b))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Rock bunting in SE Switzerland:\nOcc. prob. (post. means)", zlim = c(0, 1))

# Posterior sds
r1 <- rasterFromXYZ(data.frame(x = CovarData[1], y = CovarData[2], 
  z = psd.psi4b))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Occupancy uncertainty\n(post. SDs)", zlim = c(0, 0.36))

# Posterior means of spatial process w
r1 <- rasterFromXYZ(data.frame(x = CovarData[1], y = CovarData[2], 
  z = pm.w4b))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Spatial process in\n occupancy (post. means)", zlim = c(-5, 3.5))


# We want to compare the two maps from the spatial and nonspatial models 

# Compute difference map: non-spatial psi minus spatial psi
diff.psi <- pm.psi - pm.psi4b

par(mfrow = c(1, 3), mar = c(2, 4, 6, 10), cex.main = 2)
# Posterior means nonspatial model
r1 <- rasterFromXYZ(data.frame(x = CovarData[1], y = CovarData[2], 
  z = pm.psi))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Non-spatial Rock bunting SDM:\nOcc. prob. (post. means)", zlim = c(0, 1))

# Posterior means spatial model
r1 <- rasterFromXYZ(data.frame(x = CovarData[1], y = CovarData[2], 
  z = pm.psi4b))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Spatial Rock bunting SDM:\nOcc. prob. (post. means)", zlim = c(0, 1))

# Occupancy difference map
range(diff.psi)
r1 <- rasterFromXYZ(data.frame(x = CovarData[1], y = CovarData[2], 
  z = diff.psi))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Non-spatial psi - spatial psi", zlim = c(-0.7, 0.7))


# Let's compute the expected range size of the rock bunting in SE Switzerland. 
# For this we can simply compute the sum(psi) as a derived quantity and 
# then summarize its posterior distribution.

# Compute range size and plot (not shown)
rangeSize2 <- apply(pred.fm4b[[1]], 1, sum)
hist(rangeSize2, breaks = 40, col = 'grey', main = 'Estimated range size of the Rock Bunting in SE Switzerland\n(Posterior distribution from spatial model in spOccupancy)', xlab = 'Number of sq.kms', freq = FALSE)


# Plot posterior distributions from non-spatial and spatial next to each other
rangeSize1 <- rangeSize.spo
hist(rangeSize1, breaks = 60, col = rgb(1,0,0,0.3), main = 'Estimated range size of the Rock Bunting in SE Switzerland\n(Posterior distribution)', xlab = 'Number of sq.kms', freq = FALSE, xlim = c(1400, 2000), ylim = c(0, 0.009), cex.main = 1.6)
hist(rangeSize2, breaks = 40, col = rgb(0,0,1,0.3), freq = FALSE, add = TRUE)
legend('topright', pch = 15, col = c(rgb(0,0,1,0.3), rgb(1,0,0,0.3)), legend = c('Spatial model', 'Nonspatial model'), bty = 'n', cex = 1.5)



#  4.2 Spatial exponential model with 15-Nearest-neighbour Gaussian Process (NNGP 15)
# ----------------------------------------------------------------------------------

# We launch the model again, using much larger MCMC settings than before.

# Resultant number of draws for each chain
4 * ((500 * 400)-50000) / 600

# Launch model in spOccupancy (key change marked in grey)
fm4c <- 
    spPGOcc(occ.formula = occ.formula, 
       det.formula = det.formula, 
       data = AtlasData, inits = spo.inits, 
       n.batch = 500, batch.length = 400, 
       priors = spo.priors, cov.model = cov.model, 
       NNGP = TRUE, n.neighbors = 15, tuning = spo.tuning, 
       n.omp.threads = 6,
       n.burn = 50000, n.thin = 600, n.chains = 4, 
       n.report = 100, verbose = TRUE)

summary(fm4c)

plot(fm4c$beta.samples)    # For the occupancy params
plot(fm4c$alpha.samples)   # For the detection params
plot(fm4c$theta.samples)   # For the spatial params

# Compare the 5NNGP and the 15NNGP estimates

library(abind)

# Compute point estimates
m4bpm <- c(apply(fm4b$beta.samples, 2, mean), apply(fm4b$alpha.samples, 2, mean), apply(fm4b$theta.samples, 2, mean))
m4cpm <- c(apply(fm4c$beta.samples, 2, mean), apply(fm4c$alpha.samples, 2, mean), apply(fm4c$theta.samples, 2, mean))

# Compute 95% CRIs
qt <- function(x) quantile(x, c(0.025, 0.975))
m4bCRI <- abind(apply(fm4b$beta.samples, 2, qt), apply(fm4b$alpha.samples, 2, qt), apply(fm4b$theta.samples, 2, qt))
m4cCRI <- abind(apply(fm4c$beta.samples, 2, qt), apply(fm4c$alpha.samples, 2, qt), apply(fm4c$theta.samples, 2, qt))

# Make plot
xylim <- c(-5, 1.2)
par(mar = c(5,5,3,2), cex.lab = 1.5, cex.axis = 1.5)
plot(m4bpm, m4cpm, xlab = 'Estimates (NNGP with 5 NN)', ylab = 'Estimates (NNGP with 15 NN)', frame = FALSE, col = 'grey', pch = 1, cex = 1, xlim = xylim, ylim = xylim)
abline(0, 1, col = 'grey', lty = 2)
segments(m4bCRI[1,], m4cpm, m4bCRI[2,], m4cpm, lwd = 1, col = 'grey')
segments(m4bpm, m4cCRI[1,], m4bpm, m4cCRI[2,], lwd = 1, col = 'grey')
points(m4bpm[1:19], m4cpm[1:19], col = rgb(1,0,0, 0.4), pch = 16, cex = 1.2)
points(m4bpm[20:22], m4cpm[20:22], col = rgb(0,1,0, 0.4), pch = 16, cex = 1.2)
points(m4bpm[23:24], m4cpm[23:24], col = rgb(0,0,1, 0.4), pch = 16, cex = 1.2)

# In this figure, red is for the occupancy, green for detection and blue for spatial params.

# Form predictions from the spatial model (NNGP with 15 neighbours)
system.time(
  pred.fm4c <- predict(fm4c, spo.cov, spo.coords, verbose = TRUE)
)

# Load needed packages
require(raster)

# Map posterior mean and posterior SD of occupancy probability
# We also compute posterior means of w
pm.psi4c <- apply(pred.fm4c$psi.0.samples, 2, mean)  # Post mean
psd.psi4c <- apply(pred.fm4c$psi.0.samples, 2, sd)   # Post sd
pm.w4c <- apply(pred.fm4c$w.0.samples, 2, mean)      # Post mean process

summary(pm.psi4c)
summary(psd.psi4c)
summary(pm.w4c)

mapPalette1 <- colorRampPalette(c("grey", "yellow", "orange", "red"))

par(mfrow = c(1, 3), mar = c(4,4,6,8), cex.main = 2)
# Posterior means
r1 <- rasterFromXYZ(data.frame(x = CovarData[1], y = CovarData[2], 
  z = pm.psi4c))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Rock bunting in SE Switzerland:\nOcc. prob. (post. means, 15 NN)", zlim = c(0, 1))

# Posterior sds
r1 <- rasterFromXYZ(data.frame(x = CovarData[1], y = CovarData[2], 
  z = psd.psi4c))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Occupancy uncertainty\n(post. SDs, 15 NN)", zlim = c(0, 0.36))

# Posterior means of spatial process w
r1 <- rasterFromXYZ(data.frame(x = CovarData[1], y = CovarData[2], 
  z = pm.w4c))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Spatial process in\n occupancy (post. means, 15 NN)", zlim = c(-6, 3.5))

# We want to compare the two maps of the spatial process from the NNGP 5 and the NNGP 15 runs.

# Compute difference map: NNGP 15 minus NNGP 5
diff.w <- pm.w4c - pm.w4b
summary(diff.w)

par(mfrow = c(1, 3), mar = c(2, 4, 6, 10), cex.main = 2)
# Posterior means of spatial process w with 5 nearest neighbours
r1 <- rasterFromXYZ(data.frame(x = CovarData[1], y = CovarData[2], 
  z = pm.w4b))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Spatial process with NNGP 5", zlim = c(-6, 6))

# Posterior means of spatial process w with 15 nearest neighbours
r1 <- rasterFromXYZ(data.frame(x = CovarData[1], y = CovarData[2], 
  z = pm.w4c))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Spatial process with NNGP 15", zlim = c(-6, 6))

# Spatial process difference map
r1 <- rasterFromXYZ(data.frame(x = CovarData[1], y = CovarData[2], 
  z = diff.w))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Difference NNGP 5 minus NNGP 15", zlim = c(-2, 1))

# Not much difference








# 5 Fitting restricted spatial regression (RSR) models with ubms
# --------------------------------------------------------------

library(ubms)   # Also automatically loads unmarked

# Remind ourselves of format
str(AtlasData)
str(CovarData)

# Make unique quadrat identifier in both data sets
ID_y <- paste(AtlasData$coords[,1], AtlasData$coords[,2])
ID_cv <- paste(CovarData[,1], CovarData[,2])

# Create new data frames for reponse (detection/nondetection) and for detection covs
y_ubms <- matrix(NA, nrow = 12757, ncol = 21)
detcov_ubms <- matrix(NA, nrow = 12757, ncol = 21)
# Fill them
for(i in 1: length(ID_y)){
  sel.quad <- which(ID_cv == ID_y[i])
  y_ubms[sel.quad,] <- AtlasData$y[i,]
  detcov_ubms[sel.quad,] <- AtlasData$det.covs[[1]][i,]
}

# Sum tests ... look good (not shown)
sum(y_ubms, na.rm = TRUE) ; sum(AtlasData$y, na.rm = TRUE)
sum(detcov_ubms, na.rm = TRUE) ; sum(AtlasData$det.covs[[1]], na.rm = TRUE)

# Repackage all the data in a new unmarked data frame for occu()
spatial_umf <- unmarkedFrameOccu(
  y = y_ubms,         # Reformatted Presence/Absence measurements
  siteCovs = CovarData[-(3:19)],  # Environmental covariates at site-level
  obsCovs = list(date1 = detcov_ubms,
    date2 = detcov_ubms^2)) # Observation-specific covariates
summary(spatial_umf)

# A rook neighbourhood
site_cov <- CovarData[-(3:19)]
with(site_cov, RSR(x, y, threshold=1, plot_site = 3050))


# A queen's neighbourhood
with(site_cov, RSR(x, y, threshold=1.5, plot_site = 3050))


# A larger neighbourhood
with(site_cov, RSR(x, y, threshold=5, plot_site = 3050))


# A much larger neighbourhood
with(site_cov, RSR(x, y, threshold=10, plot_site = 3050))


# And EVEN larger neighbourhood
with(site_cov, RSR(x, y, threshold=20, plot_site = 3050))


# Double formula: first part is for detection, second for occupancy
form1 <- ~ date1 + date2 ~ z.buildings + z.buildings2 + 
  z.elevation + z.elevation2 + z.northness + z.northness2 + 
  z.rivers + z.rivers2 + z.rocks + z.rocks2 + z.slope + z.slope2 + 
  z.structures + z.structures2 + z.wetlands + z.wetlands2 +
  z.kfrivers + z.kfrivers2 + RSR(x, y, threshold = 1.5)

form2 <- ~ date1 + date2 ~ z.buildings + z.buildings2 + 
  z.elevation + z.elevation2 + z.northness + z.northness2 + 
  z.rivers + z.rivers2 + z.rocks + z.rocks2 + z.slope + z.slope2 + 
  z.structures + z.structures2 + z.wetlands + z.wetlands2 +
  z.kfrivers + z.kfrivers2 + RSR(x, y, threshold = 5)

form3 <- ~ date1 + date2 ~ z.buildings + z.buildings2 + 
  z.elevation + z.elevation2 + z.northness + z.northness2 + 
  z.rivers + z.rivers2 + z.rocks + z.rocks2 + z.slope + z.slope2 + 
  z.structures + z.structures2 + z.wetlands + z.wetlands2 +
  z.kfrivers + z.kfrivers2 + RSR(x, y, threshold = 10)

form4 <- ~ date1 + date2 ~ z.buildings + z.buildings2 + 
  z.elevation + z.elevation2 + z.northness + z.northness2 + 
  z.rivers + z.rivers2 + z.rocks + z.rocks2 + z.slope + z.slope2 + 
  z.structures + z.structures2 + z.wetlands + z.wetlands2 +
  z.kfrivers + z.kfrivers2 + RSR(x, y, threshold = 20)

# Unfortunately, on trying this, we find that ubms never gets started
# We try to reduce the computational problem by only modeling in ubms part of previous domain.

# Here's the map again of the spatial effect from spOccupancy.

# Posterior means of spatial process w, with possible geo-limits
library(raster)
r1 <- rasterFromXYZ(data.frame(x = CovarData[1], y = CovarData[2], 
  z = pm.w4c))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Spatial process in\n occupancy (post. means, 15 NN)", zlim = c(-6, 3.5))
abline(v = 750, col = 'black', lwd = 2) # Cut it here
abline(v = 790, col = 'black', lwd = 2) # Cut it here

# Current data for ubms
str(y_ubms)
str(CovarData)
str(detcov_ubms)

# Check out a new limit for defining smaller domain West of limit
plot(CovarData[,1:2], pch = '.')      # Current domain
abline(v = 750, col = 'red', lwd = 3) # Cut it here
abline(v = 790, col = 'red', lwd = 3) # Cut it here


# Restrict domain to Eastern half, with x coord >= 750
condition <- CovarData[, 1] > 749.5 & CovarData[, 1] < 790
sum(condition)               # Results in 3022 quadrats

# Make new unmarked/ubms data frame
spatial_umf2 <- unmarkedFrameOccu(
  y = y_ubms[condition,],   # Presence/Absence measurements
  siteCovs = CovarData[-(3:19)][condition,],# Site-level covariates
  obsCovs = list(date1 = detcov_ubms[condition,],
    date2 = detcov_ubms[condition,]^2)) # Obs.level covariates
summary(spatial_umf2)

# Try ubms again
system.time(
  fm51 <- stan_occu(form = form1, data = spatial_umf2, 
  chains = 3, cores = 4) )
traceplot(fm51)       # Check convergence of chains
print(fm51)           # Print posterior summaries


From here on: have to run it for much longer

system.time(
  fm52 <- stan_occu(form = form2, data = spatial_umf, 
  chains = 3, cores = 3) )
traceplot(fm52)       # Check convergence of chains
print(fm52)           # Print posterior summaries

system.time(
  fm53 <- stan_occu(form = form3, data = spatial_umf, 
  chains = 3, cores = 3) )
traceplot(fm53)       # Check convergence of chains
print(fm53)           # Print posterior summaries


system.time(
  fm54 <- stan_occu(form = form4, data = spatial_umf, 
  chains = 3, cores = 3) )
traceplot(fm54)       # Check convergence of chains
print(fm54)           # Print posterior summaries

fm2


# Repeat predictions from nonspatial model fm2X
cbind(names(CovarData))     # not shown
newdata <- data.frame(CovarData[20:53])
pred.ubms <- ubms::predict(fm2, newdata = newdata, submodel = 'state')
str(pred.ubms)
head(pred.ubms)


# Use plot_spatial()
# par(mfrow = c(2, 2))    # Does not work
plot_spatial(fm51)      # top left in 2x2 plot: 1.5
plot_spatial(fm52)      # top right: 5
plot_spatial(fm53)      # bottom left: 10
plot_spatial(fm54)      # bottom right: 20


# Use plot_spatial() for the spatial random field in the models
tmp1 <- plot_spatial(fm51, "eta") # top left in 2x2 plot: 1.5
tmp2 <- plot_spatial(fm52, "eta")      # top right: 5
tmp3 <- plot_spatial(fm53, "eta")      # bottom left: 10
tmp4 <- plot_spatial(fm54, "eta")      # bottom right: 20

# Are these any different at all ?

# Predictions from spatial models fm51 – fm54
pairs(cbind(tmp1$data$est, tmp2$data$est, tmp3$data$est, tmp4$data$est))
# OK, they are…


# Predictions from spatial models fm51 – fm54
pred.fm51 <- ubms::predict(fm51, submodel = 'state')
pred.fm52 <- ubms::predict(fm52, submodel = 'state')
pred.fm53 <- ubms::predict(fm53, submodel = 'state')
pred.fm54 <- ubms::predict(fm54, submodel = 'state')


# Rock Bunting species distribution map in SE Switzerland from ubms:

par(mfrow = c(3, 2), mar = c(2,1,4,5), cex.main = 1.3)
# Nonspatial model
r1 <- rasterFromXYZ(data.frame(x = CovarData$x, y = CovarData$y, 
  z = pred.ubms[,1]))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Nonspatial model", zlim = c(0, 1))

# Spatial model with RSR threshold 1.5
r1 <- rasterFromXYZ(data.frame(x = CovarData$x, y = CovarData$y, 
  z = pred.fm51[,1]))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Spatial model (RSR Thresh. 1.5)", zlim = c(0, 1))

# Spatial model with RSR threshold 5
r1 <- rasterFromXYZ(data.frame(x = CovarData$x, y = CovarData$y, 
  z = pred.fm52[,1]))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Spatial model (RSR Thresh. 5)", zlim = c(0, 1))

# Spatial model with RSR threshold 10
r1 <- rasterFromXYZ(data.frame(x = CovarData$x, y = CovarData$y, 
  z = pred.fm53[,1]))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Spatial model (RSR Thresh. 10)", zlim = c(0, 1))

# Spatial model with RSR threshold 20
r1 <- rasterFromXYZ(data.frame(x = CovarData$x, y = CovarData$y, 
  z = pred.fm54[,1]))
plot(r1, col = mapPalette1(100), axes = FALSE, box = FALSE, main = "Spatial model (RSR Thresh. 20)", zlim = c(0, 1))

# These look fairly similar.