# 11e-summary-hbef-sfMsNMix.R: script to summarize the results from the spatial
#                              multi-species N-mixture model with four species
#                              in the Hubbard Brook Experimental Forest. There
#                              isn't much done in this script, it is more so
#                              just to show the workflow for fitting models with
#                              multiple chains in this way.
rm(list = ls())
library(spAbundance)

# Load the full model object with the three chains ------------------------
# Object loaded is called out.full
load('results/hbef-sfMsNmix-three-chains.rda')

# Has the model converged?
summary(out.full)
plot(out.full, 'beta', density = FALSE)
plot(out.full, 'alpha', density = FALSE)

# Summarize the results in whatever manner you like....
