# 11a-crimson-mantled-woodpecker-nmix.R: script to fit a single-species
#                                        N-mixture model to estimate how
#                                        abundance of the crimson-mantled woodpecker
#                                        varies across an elevation gradient in
#                                        the Andes mountain in western Bolivia.
#                                        These data are part of a very nice study
#                                        by Montaño‐Centellas et al. that looks at
#                                        drivers of avian community assembly across
#                                        a tropical elevation gradient.
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


# TODO 10a.2: fit the same model, but now include a linear effect of elevation.
#             Make sure to standardize elevation.


# TODO 10a.3: fit the same model, but now include a linear and quadratic effect of elevation.
#             Make sure to standardize elevation.


# TODO 10a.4: use the waicAbund() function to calculate the WAIC and compare the three models.
#             Which model performs the best?


# TODO 10a.5: use the ppcAbund() function to perform a posterior predictive check
#             for the top-performing model. Do the models adequately fit the data?


# TODO 10a.6: use the fitted() function to extract the fitted values for the
#             observed counts and generate a plot of the mean fitted values vs.
#             the observed count values. How does the plot look?


# TODO 10a.7: using the predict() function, generate a conditional effects plot
#             showing the relationship between elevation and abundance.

