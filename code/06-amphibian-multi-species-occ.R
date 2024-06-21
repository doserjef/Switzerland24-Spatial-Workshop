# 06-amphibian-multi-species-occ.R: this script showcases how to fit multi-species
#                                   occupancy models in spOccupancy. Here we focus
#                                   on a community of tropical amphibians and how
#                                   their occupancy relates to agriculture and
#                                   topography. These data come from Ribeiro Jr.
#                                   et al. (2018). We fit a standard multi-species
#                                   occupancy model as well as a multi-species
#                                   occupancy model that accounts for residual
#                                   correlations between species in a joint
#                                   species distribution modeling framework.
# Data source citation:
#   Ribeiro Jr, J. W., Siqueira, T., Brejão, G. L., & Zipkin, E. F. (2018).
#   Effects of agriculture and topography on tropical amphibian species
#   and communities. Ecological Applications, 28(6), 1554-1564.
rm(list = ls())
library(spOccupancy)
# For plotting and summarizing results
library(MCMCvis)
library(ggplot2)
library(dplyr)
library(corrplot)
# If not using the RStudio project, set working directory to the repository
# directory.
# setwd("../")
set.seed(9857)

# Our goal in this script is to understand the effects of five landscape variables
# on the occupancy of 36 amphibian species: (1) forest cover; (2) agricultural cover;
# (3) stream catchment area (proxy for many local habitat characteristics of streams);
# (4) stream density within a 200m radius buffer; and (5) slope. Note that here
# we are using data that solely come from passive acoustic monitoring, while in
# the full manuscript, they also leverage data from active transect surveys. Here
# we only work with the subset of species that were ever detected with an ARU.

# Data prep ---------------------------------------------------------------
# Read in the data source (reads in an object called data.list)
load("data/multiSpeciesRibeiroJr2018EcoApps.rda")
# Check out the data list structure. Detection-nondetection data y
# are currently formatted for multiple species models.
str(data.list)
# The structure of the detection-nondetection data y is a three-dimensional
# array, where the first dimension corresponds to the species, second dimension
# corresponds to site, and third dimension corresponds to repeat visits. Formatting
# the data into this array can be one of the trickier things when fitting multi-species
# models in spOccupancy. See the below vignette for an example of how to do that.
# https://www.jeffdoser.com/files/spoccupancy-web/articles/dataformatting.

# How many detections of each species?
apply(data.list$y, 1, sum, na.rm = TRUE)
# Lots of rare species in this data set!

# Model fitting -----------------------------------------------------------
# Fit a non-spatial, multi-species occupancy model without residual species
# correlations. Notice msPGOcc has exactly the same arguments as PGOcc. We will
# use all the default priors and initial values.
out <- msPGOcc(occ.formula = ~ scale(forest) + scale(agriculture) +
                               scale(catchment) + scale(density) +
                               scale(slope),
               det.formula = ~ scale(date) + I(scale(date)^2) + scale(rain),
               data = data.list,
               n.samples = 8000,
               n.thin = 4,
               n.burn = 4000,
               n.chains = 3,
               n.report = 500)
# Quick summary of community-level parameters
summary(out, level = 'community')
# Quick summary of species-level parameters
summary(out, level = 'species')
# Fit a multi-species occupancy model that accounts for residual species correlations
# using a latent factor modeling approach. All arguments are the same as with
# msPGOcc, just also need to specify the number of latent factors to use.
out.lf <- lfMsPGOcc(occ.formula = ~ scale(forest) + scale(agriculture) +
                                    scale(catchment) + scale(density) +
                                    scale(slope),
                    det.formula = ~ scale(date) + I(scale(date)^2) + scale(rain),
                    data = data.list,
                    n.factors = 4,
                    n.samples = 8000,
                    n.thin = 4,
                    n.burn = 4000,
                    n.chains = 3,
                    n.report = 500)
summary(out.lf)

# Model validation --------------------------------------------------------
# Perform a posterior predictive check to assess model fit.
ppc.out <- ppcOcc(out, fit.stat = 'freeman-tukey', group = 1)
ppc.out.lf <- ppcOcc(out.lf, fit.stat = 'freeman-tukey', group = 1)
# Calculate a Bayesian p-value as a simple measure of Goodness of Fit.
# Bayesian p-values between 0.1 and 0.9 indicate adequate model fit.
summary(ppc.out)
summary(ppc.out.lf)

# Model comparison --------------------------------------------------------
# Compute Widely Applicable Information Criterion (WAIC)
# Lower values indicate better model fit.
waicOcc(out)
waicOcc(out.lf)
# For multi-species models, we can also extract the WAIC individually for
# each species (the overall WAIC is just the sum of all species-specific WAIC values)
waicOcc(out, by.sp = TRUE)
waicOcc(out.lf, by.sp = TRUE)

# Posterior summaries -----------------------------------------------------
# Concise summary of main parameter estimates
summary(out.lf, level = 'community')
summary(out.lf, level = 'species')
summary(out.lf, level = 'both')
# Create simple plot summaries using MCMCvis package.
# Occupancy community-level effects
MCMCplot(out.lf$beta.comm.samples, ref_ovl = TRUE, ci = c(50, 95))
# Detection covariate effects ---------
MCMCplot(out.lf$alpha.comm.samples, ref_ovl = TRUE, ci = c(50, 95))
# Extract all occupancy parameters for one species of interest
dimnames(data.list$y)[[1]]
MCMCplot(out.lf$beta.samples, ref_ovl = TRUE, ci = c(50, 95),
         param = 'tADEMAR', exact = FALSE)

# Plot median effect of agriculture for all species and community ---------
# Species-level agriculture effect
beta.ag.samples <- MCMCchains(out.lf$beta.samples, param = 'agriculture', exact = FALSE)
# Community-level agriculture effect
beta.ag.comm.samples <- MCMCchains(out.lf$beta.comm.samples, param = 'agriculture',
                                   exact = FALSE)
# Species codes, plus the community code (COMM)
sp.names <- dimnames(data.list$y)[[1]]
sp.names.factor <- factor(sp.names, levels = sp.names)
sp.codes <- as.numeric(sp.names.factor)
N <- length(sp.names)
cov.plot.df <- data.frame(ag.mean = apply(beta.ag.samples, 2, mean),
                          ag.low = apply(beta.ag.samples, 2, quantile, 0.25),
                          ag.lowest = apply(beta.ag.samples, 2, quantile, 0.025),
                          ag.high = apply(beta.ag.samples, 2, quantile, 0.75),
                          ag.highest = apply(beta.ag.samples, 2, quantile, 0.975),
                          sp = sp.codes)
# Rearrange and add things to get the plot to display effects in increasing order.
cov.plot.df <- cov.plot.df %>%
  arrange(ag.mean)
cov.plot.df$sp.factor <- as.character(sp.names.factor[cov.plot.df$sp])
cov.plot.df$sort.sp <- 1:N

# Add in the community level covariate
comm.plot.df <- data.frame(ag.mean = mean(beta.ag.comm.samples),
                           ag.low = quantile(beta.ag.comm.samples, 0.25),
                           ag.lowest = quantile(beta.ag.comm.samples, 0.025),
                           ag.high = quantile(beta.ag.comm.samples, 0.75),
                           ag.highest = quantile(beta.ag.comm.samples, 0.975),
                           sp = N + 1,
                           sp.factor = 'COMM',
                           sort.sp = N + 1)
cov.plot.df <- rbind(cov.plot.df, comm.plot.df)
cov.plot.df$sp.factor <- factor(cov.plot.df$sp.factor, levels = c(unique(cov.plot.df$sp.factor)))

# Make the figure
ag.cov.plot <- ggplot(data = cov.plot.df, aes(x = sort.sp, fill = ag.mean, group = sp.factor)) +
  geom_hline(yintercept = 0, col = 'black', linewidth = 0.75, lty = 2) +
  geom_boxplot(aes(ymin = ag.lowest, lower = ag.low, middle = ag.mean,
                   upper = ag.high, ymax = ag.highest), stat = 'identity', col = 'black') +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC',
                       na.value = NA) +
  theme_bw(base_size = 14) +
  guides(fill = "none") +
  labs(x = "Species", y = "Effect of Agriculture") +
  scale_x_continuous(breaks = 1:(N+1), labels = cov.plot.df$sp.factor) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ag.cov.plot

# Generate conditional probability plots ----------------------------------
# Plot relationship with forest cover for each species
# Predict occupancy along a gradient of forest cover.
# Create a set of values across the range of observed forest values
forest.pred.vals <- seq(min(data.list$occ.covs$forest),
                        max(data.list$occ.covs$forest),
                        length.out = 100)
# Scale predicted values by mean and standard deviation used to fit the model
forest.pred.vals.scale <- (forest.pred.vals - mean(data.list$occ.covs$forest)) /
  sd(data.list$occ.covs$forest)
# Predict occupancy across forest values at mean values of all other variables
pred.df <- as.matrix(data.frame(intercept = 1, forest = forest.pred.vals.scale,
                                agriculture = 0, catchment = 0, density = 0,
                                slope = 0))
# NOTE: to generate marginal effects plots with a latent factor model without
#       needing to specify any coordinates, we need to do a bit of a workaround
#       and first convert the class of the model to "msPGOcc", then use predict().
#       I am in the process of adding in an argument to predict() to make this
#       more streamlined and less clunky.
# Convert class to msPGOcc
class(out.lf) <- 'msPGOcc'
out.pred <- predict(out.lf, pred.df)
# Convert class back to lfMsPGOcc to avoid any problems down the road.
class(out.lf) <- 'lfMsPGOcc'
str(out.pred)
psi.0.quants <- apply(out.pred$psi.0.samples, c(2, 3), quantile,
                      prob = c(0.025, 0.5, 0.975))
sp.codes <- attr(data.list$y, "dimnames")[[1]]
psi.plot.dat <- data.frame(psi.med = c(t(psi.0.quants[2, , ])),
                           psi.low = c(t(psi.0.quants[1, , ])),
                           psi.high = c(t(psi.0.quants[3, , ])),
                           forest = forest.pred.vals,
                           sp.codes = rep(sp.codes,
                                          each = length(forest.pred.vals)))
ggplot(psi.plot.dat, aes(x = forest, y = psi.med)) +
  geom_ribbon(aes(ymin = psi.low, ymax = psi.high), fill = 'grey70') +
  geom_line() +
  facet_wrap(vars(sp.codes)) +
  theme_bw() +
  labs(x = 'Forest (% cover)', y = 'Occupancy Probability')

# Generate residual correlation matrix ------------------------------------
# Species-species residual covariance matrix ---------
# Factor loadings
lambda.means <- matrix(apply(out.lf$lambda.samples, 2, mean), nrow(out.lf$y), out.lf$q)
lambda.means
# Calculate lambda %*% lambda as a singular species-species covariance matrix
Sigma.means <- lambda.means %*% t(lambda.means)
# Singular correlation matrix
species.cor.mat <- cov2cor(Sigma.means)
rownames(species.cor.mat) <- dimnames(data.list$y)[[1]]
colnames(species.cor.mat) <- dimnames(data.list$y)[[1]]
# Remember, this is not not indicative of true biotic interactions! It is a
# residual species correlation matrix. Species with a negative relationship tend
# to not occur together after accounting for the effects of the covariates, while
# species with a positive correlation tend to occurr together after accounting
# for the effects of the covariates.
corrplot(species.cor.mat, method = 'square', type = 'lower', order = 'hclust')
