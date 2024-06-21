# 11a-solutions-crimson-mantled-woodpecker-nmix.R: script to fit a single-species
#                                                  N-mixture model to estimate how
#                                                  abundance of the crimson-mantled woodpecker
#                                                  varies across an elevation gradient in
#                                                  the Andes mountain in western Bolivia.
#                                                  These data are part of a very nice study
#                                                  by Montaño‐Centellas et al. that looks at
#                                                  drivers of avian community assembly across
#                                                  a tropical elevation gradient.
# Data source citation:
#   Montaño‐Centellas, F. A., Loiselle, B. A., & Tingley, M. W. (2021).
#   Ecological drivers of avian community assembly along a tropical elevation
#   gradient. Ecography, 44(4), 574-588.

rm(list = ls())
library(spAbundance)
library(ggplot2)
set.seed(3829)

# Read in the data set ----------------------------------------------------
# Loads in an object called data.list
load('data/montano-centellas-bird-data.rda')
str(data.list)

# The data we are using here are a subset of data used in the above referenced
# manuscript. Here our goal is to fit and compare three models: (1) null model that
# assumes abundance is constant across the elevation gradient; (2) linear only
# relationship with elevation; and (3) a quadratic relationship with elevation.
# In all models, we will include linear and quadratic effects of julian day and
# survey start time on detection probability, as well as a linear effect of effort
# and the direction of the transect (up or down slope) on detection probability.

# TODO 10a.1: use the NMix() function to fit a single-species N-mixture model
#             assuming abundance is constant. Use a Poisson distribution. Manually set the
#             prior distribution for the abundance regression coefficients to
#             Normal(0, var = 10) and for the detection coefficients Normal(0, 2.72).
n.batch <- 4000
batch.length <- 25
n.burn <- 50000
n.thin <- 20
n.chains <- 3
priors <- list(beta.normal = list(mean = 0, var = 10),
               alpha.normal = list(mean = 0, var = 2.72))
# Approx run time: <1 minutes
out.1 <- NMix(abund.formula = ~ 1,
              det.formula = ~ jday.s + I(jday.s^2) + tod.s + I(tod.s^2) + effort.s + factor(direction),
              data = data.list, priors = priors, family = 'Poisson',
              n.batch = n.batch, batch.length = batch.length, n.burn = n.burn,
              n.thin = n.thin, n.chains = n.chains, n.report = 400)
# Check convergence
summary(out.1)
plot(out.1, 'beta')

# TODO 10a.2: fit the same model, but now include a linear effect of elevation.
#             Make sure to standardize elevation.
# Approx run time: <1 minute
out.2 <- NMix(abund.formula = ~ scale(elev),
              det.formula = ~ jday.s + I(jday.s^2) + tod.s + I(tod.s^2) + effort.s + factor(direction),
              data = data.list, priors = priors, family = 'Poisson',
              n.batch = n.batch, batch.length = batch.length, n.burn = n.burn,
              n.thin = n.thin, n.chains = n.chains, n.report = 1000)
# Check convergence
summary(out.2)
plot(out.2, 'beta')

# TODO 10a.3: fit the same model, but now include a linear and quadratic effect of elevation.
#             Make sure to standardize elevation.
# Approx run time: <1 minute
out.3 <- NMix(abund.formula = ~ scale(elev) + I(scale(elev)^2),
              det.formula = ~ jday.s + I(jday.s^2) + tod.s + I(tod.s^2) + effort.s + factor(direction),
              data = data.list, priors = priors, family = 'Poisson',
              n.batch = n.batch, batch.length = batch.length, n.burn = n.burn,
              n.thin = n.thin, n.chains = n.chains, n.report = 1000)
# Check convergence
summary(out.3)
plot(out.3, 'beta')

# TODO 10a.4: use the waicAbund() function to calculate the WAIC and compare the three models.
#             Which model performs the best?
waicAbund(out.1)
waicAbund(out.2)
waicAbund(out.3)

# TODO 10a.5: use the ppcAbund() function to perform a posterior predictive check
#             for the top-performing model. Do the models adequately fit the data?
# Don't group the data before calculation of the fit statistic
ppc.out.0 <- ppcAbund(out.3, fit.stat = 'chi-square', group = 0)
# Group the data by site before calculation of the fit statistic
ppc.out.1 <- ppcAbund(out.3, fit.stat = 'chi-square', group = 1)
# Group the data by visit before calculation of the fit statistic
ppc.out.2 <- ppcAbund(out.3, fit.stat = 'chi-square', group = 2)
summary(ppc.out.0)
summary(ppc.out.1)
summary(ppc.out.2)
par(mfrow = c(1, 1))
plot(ppc.out.1$fit.y, ppc.out.1$fit.y.rep, pch = 19)
abline(0, 1)

# TODO 10a.6: use the fitted() function to extract the fitted values for the
#             observed counts and generate a plot of the mean fitted values vs.
#             the observed count values. How does the plot look?
fitted.vals <- fitted(out.3)
y.rep.means <- apply(fitted.vals$y.rep.samples, c(2, 3), mean, na.rm = TRUE)
plot(c(data.list$y), c(y.rep.means), pch = 19)
abline(0, 1)

# TODO 10a.7: using the predict() function, generate a conditional effects plot
#             showing the relationship between elevation and abundance.
# Get prediction elevation values across the range of the observed data.
elev.pred.vals <- seq(min(data.list$abund.covs$elev),
                      max(data.list$abund.covs$elev),
                      length.out = 100)
# Scale predicted values by mean and standard deviation used to fit the model
elev.pred.vals.scale <- (elev.pred.vals - mean(data.list$abund.covs$elev)) /
                        sd(data.list$abund.covs$elev)
# Predict abundance across elevation gradient
pred.df <- as.matrix(data.frame(intercept = 1, elev = elev.pred.vals.scale,
                                elev.2 = elev.pred.vals.scale^2))
out.pred <- predict(out.3, pred.df)
str(out.pred)
# Quantiles and medians of predicted abundance across the elevation gradient.
mu.0.quants <- apply(out.pred$mu.0.samples, 2, quantile, prob = c(0.025, 0.5, 0.975))
mu.plot.dat <- data.frame(mu.med = mu.0.quants[2, ],
                          mu.low = mu.0.quants[1, ],
                          mu.high = mu.0.quants[3, ],
                          elev = elev.pred.vals)
# Make a plot
ggplot(mu.plot.dat, aes(x = elev, y = mu.med)) +
  geom_ribbon(aes(ymin = mu.low, ymax = mu.high), fill = 'grey70') +
  geom_line() +
  theme_bw() +
  labs(x = 'Elevation (m)', y = 'Expected abundance')
