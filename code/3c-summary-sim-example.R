# 3c-summary-sim-example.R: script to summarize the simulation results and
#                           generate a plot to compare the NNGP model results
#                           to the full GP model results.
rm(list = ls())
library(spOccupancy)
library(ggplot2)
library(sf)
library(patchwork)

# Load the model results --------------------------------------------------
# GP model results
load("results/spPGOcc-GP-sims.rda")
# NNGP model results
load("results/spPGOcc-NNGP-sims.rda")
# The true spatial random effects that we simulated
w.true <- c(dat$w)
# Plot the results --------------------------------------------------------
# Coordinates for plotting
coords.sf <- st_as_sf(data.frame(x = dat$coords[, 1], y = dat$coords[, 2]),
                      coords = c('x', 'y'))
# Data frame for plotting
plot.df <- data.frame(x = dat$coords[, 1],
                      y = dat$coords[, 2],
                      true = w.true,
                      nngp = w.means.nngp,
                      gp = w.means.gp)

true.plot <- ggplot(data = plot.df, aes(x = x, y = y, fill = true)) +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC',
                       na.value = NA) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_raster() +
  theme_light(base_size = 16) +
  labs(x = 'Easting', y = 'Northing', title = '(A) True w') +
  guides(fill = 'none') +
  theme(text = element_text(family="LM Roman 10"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 16),
        strip.text.x = element_text(color = 'black'),
        strip.text.y = element_text(color = 'black'))

gp.plot <- ggplot(data = plot.df, aes(x = x, y = y, fill = gp)) +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC',
                       na.value = NA) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_raster() +
  theme_light(base_size = 16) +
  labs(x = 'Easting', y = 'Northing', title = '(B) Full GP (424 min)') +
  guides(fill = 'none') +
  theme(text = element_text(family="LM Roman 10"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 16),
        strip.text.x = element_text(color = 'black'),
        strip.text.y = element_text(color = 'black'))

nngp.plot <- ggplot(data = plot.df, aes(x = x, y = y, fill = nngp)) +
  scale_fill_gradient2(midpoint = 0, low = '#B2182B', mid = 'white', high = '#2166AC',
                       na.value = NA) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_raster() +
  theme_light(base_size = 16) +
  labs(x = 'Easting', y = 'Northing', title = '(C) NNGP (7 min)') +
  guides(fill = 'none') +
  theme(text = element_text(family="LM Roman 10"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 16),
        strip.text.x = element_text(color = 'black'),
        strip.text.y = element_text(color = 'black'))

true.plot + gp.plot + nngp.plot
