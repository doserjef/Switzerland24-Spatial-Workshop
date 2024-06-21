# 12b-hbef-main-sfMsNMix.R: script to fit multi-species spatial N-mixture models
#                           with data from the Hubbard Brook Experimental Forest.
#                           This script is used to run models using initial
#                           values that were obtained from a shorter model fit
#                           in 12a-get-inits.R. This script also shows one way
#                           (which I use very often) to run multiple MCMC chains
#                           in parallel using spAbundance/spOccupancy. In particular,
#                           we will see how to run three separate instances of this
#                           R script simultaneously by running the scripts through the
#                           command line. We can access the command line via the
#                           "Terminal" window of RStudio.
# Clear out the workspace as usual
rm(list = ls())
library(spAbundance)

# Get chain number from command line run ----------------------------------
# The commandArgs() function allows us to take information provided to R when
# running the script through the command line and subsequently access it in R.
# Here we will send in the chain number through the command line to avoid having
# to run three separate scripts. This is also a good way to run R scripts on
# a server or High Performance Computer. If just running this through R Studio,
# the "chain" variable will be empty
chain <- as.numeric(commandArgs(trailingOnly = TRUE))
# If not running the script from the command line, set the chain number manually:
# NOTE: make sure this line is commented out if running from the command line.
# chain <- 1
# A useful check to make sure the chain is properly provided.
if (length(chain) == 0) base::stop('Need to tell spAbundance the chain number')

# Set seed ----------------------------------------------------------------
# Let's set a seed to make sure we get the same results. Note that we need to
# set a different seed for each chain, otherwise we will get exactly the same
# MCMC results.
if (chain == 1) {
  set.seed(88839)
} else if (chain == 2) {
  set.seed(3245)
} else if (chain == 3) {
  set.seed(11884)
}

# Data prep ---------------------------------------------------------------
# This is all the same as 11a-get-inits.R
data(hbefCount2015)

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
prior.list <- list(beta.comm.normal = list(mean = 0, var = 10),
                   alpha.comm.normal = list(mean = 0, var = 2.72),
                   phi.unif = list(a = 3 / 7000, b = 3 / 500))
# Initial values ----------------------
# Load the initial values from our initial model fit (loads inits.list).
load('../data/hbef-sfMsNMix-inits.rda')
# Add a small amount of noise to the occupancy and detection coefficients. This
# is a useful technique to use for models that are sensitive to initial values,
# but you still want to allow for some amount of randomness in the initial values
# that are used across different chains. Note that I use the same initial values for
# phi, lambda, and N (the latent abundance values) across all three chains.
n.sp <- nrow(data.small.hbef$y)
inits.list$beta.comm <- inits.list$beta.comm + runif(length(inits.list$beta.comm), -0.2, 0.2)
inits.list$tau.sq.beta <- inits.list$tau.sq.beta + runif(length(inits.list$tau.sq.beta), 0.0, 0.2)
inits.list$alpha.comm <- inits.list$alpha.comm + runif(length(inits.list$alpha.comm), -0.2, 0.2)
inits.list$tau.sq.alpha <- inits.list$tau.sq.alpha + runif(length(inits.list$tau.sq.alpha), 0, 0.2)
inits.list$beta <- inits.list$beta + runif(length(inits.list$beta), -0.2, 0.2)
inits.list$alpha <- inits.list$alpha + runif(length(inits.list$alpha), -0.2, 0.2)
inits.list$N <- NULL

# Specify tuning values
tuning.list <- list(beta = 0.1, alpha = 0.1, phi = 1, lambda = 0.5, w = 0.5)

# Run the model -----------------------------------------------------------
n.factors <- 1
n.batch <- 1600
batch.length <- 25
n.burn <- 20000
n.thin <- 10
# Note that we are setting n.chains = 1 here, and will run multiple chains
# sequentially by running this script from the command line.
n.chains <- 1
out <- sfMsNMix(abund.formula = ~ scale(elev) + I(scale(elev)^2),
                det.formula = ~ scale(tod) + scale(day) + I(scale(day)^2),
                data = data.small.hbef, inits = inits.list, priors = prior.list,
                cov.model = 'exponential', NNGP = TRUE, n.neighbors = 5,
                n.factors = n.factors, n.batch = n.batch, tuning = tuning.list,
                batch.length = batch.length, n.burn = n.burn,
                n.thin = n.thin, n.chains = n.chains, n.report = 100)

# Save the resulting spAbundance object to your hard drive ----------------
# Note the use of paste0 to save an object with the chain number in the file name.
save(out, file = paste0('../results/hbef-sfMsNMix-results-chain-', chain, '.rda'))
