
# Demo of simple N-mixture model of Royle, Biometrics, 2004
# ---------------------------------------------------------

# Spatial hierarchical modeling workshop, Sempach, 24-27 June 2024

# Assume we're studying asp vipers (a snake species) in a nature reserve 
# that measures 1 km2, and we partition it into 100 "sites" of 1 ha each.

# Then, we go out a couple of times and count every viper seen and from this want to
# estimate the average density (per 1 ha), the number in each site, and 
# the total population size.
# We also wonder about environmental drives of density variation in space.

# We simulate data under this model and this should clarify the definition of the model.
# To demonstrate the model, we will use maximum likelihood with package 'unmarked',
# but we also show how to fit the simplest Nmix model with package 'spAbundance'

library(unmarked)
library(spAbundance)

set.seed(1)

# Select constants of data simulation
nSites <- 100       # Number of sites
lambda <- 4         # Average number of asps per ha
p <- 0.4            # Detection probability for single snake during single visit

# Simulate true latent state
# Note sites are independent
# ?rpois          # Look up function
N <- rpois(n = nSites, lambda = lambda)         # Spatial variability in abundance

# Look at the truth
N
mean(N)  ;  var(N)  ;  cat('Disp.index' = var(N) / mean(N))
table(N)
plot(sort(N))
sum(N)              # Total abundance in NR

# Simulate observation process with only one of the two measurement errors: false negatives
# There must not be any false positives (i.e., double counts during the same visit)
# ?rbinom         # Look up function
C1 <- rbinom(n = nSites, size = N, prob = p)

# Look at simulated counts and compare with truth
C1
table(C1)
table(N)
mean(C1/N, na.rm = TRUE)
cbind('true' = N, 'observed' = C1)
sum(N)
sum(C1)

# Try to fit an Nmix model to single-visit counts
data.file <- unmarkedFramePCount(y = as.matrix(C1))
summary(data.file)
fm1 <- pcount(formula = ~ 1 ~ 1, data = data.file)
summary(fm1)
backTransform(fm1, "state")
backTransform(fm1, "det")
# We see these estimates are crap. Indeed, this model is not estimable for
# single-visit data. (Although, to really study identifiability, we would 
# have to repeat the data simulation/data analysis cycle many times, e.g., 1000 times,
# and then look at the distribution of these estimates.)


# Now simulate data from a second visit to each site
C2 <- rbinom(nSites, size = N, prob = p)

# Combine the two counts per site and compare stuff
C <- cbind(C1, C2)          # Counts in a nSite x nVisits matrix
C
cbind('true' = N, C)
sum(N)  ;  sum(C[,1])  ;  sum(C[,2])

# Look at maximum count per site
Cmax <- apply(C, 1, max)
sum(Cmax)
cor(C)


# Try again to fit the Nmix now (with unmarked)
summary(data.file <- unmarkedFramePCount(y = C))
summary(fm2 <- pcount(formula = ~ 1 ~ 1, data = data.file))
backTransform(fm2, "state")
backTransform(fm2, "det")
# Now this looks good. Indeed, with two visits at each site, there should not be
# any problems unless we try to fit a too complex model.



# .... and with spAbundance
# -------------------------

# Data set
data.list <- list(y = C)

# Priors
prior.list <- list(beta.normal = list(mean = rep(0, 1), var = rep(100, 1)),
                   alpha.normal = list(mean = rep(0, 1), var = rep(2.72, 1))) 
# Starting values
inits.list <- list(alpha = 0, beta = 0, N = apply(C, 1, max, na.rm = TRUE))

# MCMC settings
n.batch <- 1000
batch.length <- 25
n.burn <- 500
n.thin <- 5
n.chains <- 3

out2 <- NMix(abund.formula = ~ 1,
            det.formula = ~ 1, 
            data = data.list, 
            n.batch = n.batch, 
            batch.length = batch.length, 
            inits = inits.list, 
            priors = prior.list, 
            accept.rate = 0.43, 
            n.omp.threads = 1, 
            verbose = TRUE, 
            n.report = 50,
            n.burn = n.burn,
            n.thin = n.thin,
            n.chains = n.chains) 

# Check convergence visually
plot(out2$alpha.samples)
plot(out2$beta.samples)
plot(out2, param = 'alpha')
plot(out2, param = 'beta')

# Posterior summaries
summary(out2)

# Compare with MLEs
fm2

# There are slight differences in the estimates. We think that they are 
# insignificant for all practical purposes. They must be due to the priors
# in spAbundance, which currently are a little informative.
# (And in addition, the prior distributions are a little skewed, and this
# will add to some discrepancy.)




# What about covariates in the Nmix model ?
# -----------------------------------------

# Of course we can have covariates in the model !
# Remember this is simply a HM that is a combination of
# a Poisson GLM with a binomial GLM (a.k.a., logistic regresssion)

# N[i] ~ Poisson(lambda[i])
# log(lambda[i]) = beta0 + beta1 * rocks[i]

# C[i,t] ~ Binomial(N[i], p[i,t])
# logit(p[i,t]) = alpha0 + alpha1 * vegHt[i]


# Assume that we have a measure of how many rocks are in each site
# (Asp vipers love rocks)
rocks <- runif(n = nSites, min = -2, max = 2)
beta0 <- 1
beta1 <- 1
lambda <- exp(beta0 + beta1 * rocks)
N <- rpois(n = nSites, lambda = lambda)         # Spatial variability in abundance

# Look at data now
table(N)
mean(N)  ;  var(N)  ;  cat('Disp.index' = var(N) / mean(N))

# Simulate observation process
# Vegetation height negatively affects detection prob
vegHt <- runif(nSites, min = -2, max = 2)
alpha0 <- 0
alpha1 <- -1
p <- plogis(alpha0 + alpha1 * vegHt)
mean(p)            # Check average p over all sites (and visits)

# Simulate two repeat visits
C <- matrix(NA, nrow = nSites, ncol = 2)
for(t in 1:2){
  C[,t] <- rbinom(n = nSites, size = N, prob = p)
}

# Compare truth and measurements, and covariates
cbind('true' = N, 'observed' = C, 
  'rocks' = round(rocks,2), 'vegHt' = round(vegHt, 2))


# Fit the Nmix now with covariates
data.file <- unmarkedFramePCount(y = C, siteCovs = data.frame(rocks = rocks, vegHt = vegHt))
summary(data.file)
fm3 <- pcount(formula = ~ vegHt ~ rocks, data = data.file)
summary(fm3)
coef(fm3)                      # Just list the MLEs
confint(fm3, type = "state")   # 95% CIs
confint(fm3, type = "det")

# Compare with truth
cbind('beta0' = beta0, 'beta1' = beta1,        # coefs of abundance model 
  'alpha0' = alpha0, 'alpha1' = alpha1)        # coefs of detection model

# These look reasonably similar, and the true values of all four parameters
# will usually lie within the 95% CIs

# Estimate of abundance intercept on the 'normal' scale
exp(coef(fm3)[1])

# Estimate of detection intercept on the 'normal' scale
plogis(coef(fm3)[3])



# Fit the Nmix model with covariates with spAbundance
# ---------------------------------------------------

# Data set
abund.covs <- data.frame(rocks = rocks)
det.covs <- data.frame(vegHt = vegHt)
data.list <- list(y = C, abund.covs = abund.covs, det.covs = det.covs)

# Priors
prior.list <- list(beta.normal = list(mean = rep(0, 2), var = rep(100, 2)),
                   alpha.normal = list(mean = rep(0, 2), var = rep(2.72, 2))) 
# Starting values
inits.list <- list(alpha = 0, beta = 0, N = apply(C, 1, max, na.rm = TRUE))

# MCMC settings
n.batch <- 1000
batch.length <- 25
n.burn <- 500
n.thin <- 5
n.chains <- 3

out3 <- NMix(abund.formula = ~ rocks,
            det.formula = ~ vegHt, 
            data = data.list, 
            n.batch = n.batch, batch.length = batch.length, 
            inits = inits.list, 
            priors = prior.list, 
            accept.rate = 0.43, 
            n.omp.threads = 1, 
            verbose = TRUE, 
            n.report = 100,
            n.burn = n.burn, n.thin = n.thin, n.chains = n.chains) 

# Check convergence visually
plot(out3, param = 'alpha')
plot(out3, param = 'beta')

# Posterior summaries
summary(out3)

# Compare with MLEs
fm3


# And now we go on and add a site-specific random effect, w[i],
# to the model for abundance
# and impose a spatial covariance model on that ....

# N[i] ~ Poisson(lambda[i])
# log(lambda[i]) = beta0 + beta1 * rocks[i] + w[i]

# C[i,t] ~ Binomial(N[i], p[i,t])
# logit(p[i,t]) = alpha0 + alpha1 * vegHt[i]

