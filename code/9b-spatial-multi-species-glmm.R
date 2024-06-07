# 9b-spatial-multi-species-glmm.R: in this script we will fit a multi-species
#                                  spatial GLMM to estimate relative abundance of 88 bird
#                                  species across the state of North Carolina, USA in 2019.
#                                  Using predictions of abundance across the state
#                                  for each species, we will then generate a map
#                                  of Shannon's diversity across the state as a
#                                  derived quantity from the multi-species Bayesian
#                                  model. These data come from the North American
#                                  Breeding Bird Survey.

# Data source citation:
#    Ziolkowski Jr., D.J., Lutmerding, M., English, W.B., Aponte, V.I.,
#    and Hudson, M-A.R., 2023, North American Breeding Bird Survey
#    Dataset 1966 - 2022: U.S. Geological Survey data release, https://doi.org/10.5066/P9GS9K64.
