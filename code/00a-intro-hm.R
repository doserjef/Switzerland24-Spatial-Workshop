
# Some code for HM in ecology intro on Monday

# Marc Kéry, 20 June 2024

# Load AHMbook library to get some demo functions for 
# (1) concepts of point pattern, abundance, and presence/absence
# and (2) of observation errors

library(AHMbook)

# Functions to demonstrate PP, N and z
?sim.fn
?simPPe

# Functions to demonstrate measurement errors: imperfect detection
?simNpC
?simPOP

# Fundamental thing is the point pattern
# N and z are derived quantities (aggregates of PP)
# Two things affect true N and true z:
# (1) Intensity of the process (i.e., density)
# (2) and spatial scale (i.e., size of spatial aggregation unit)

# Look at this

# Homogeneous point pattern
# --------------------------
# Increase density with constant spatial scale
set.seed(1)
sim.fn(quad.size = 10, cell.size = 1, intensity = 0.05)
sim.fn(quad.size = 10, cell.size = 1, intensity = 0.1)
sim.fn(quad.size = 10, cell.size = 1, intensity = 0.5)
sim.fn(quad.size = 10, cell.size = 1, intensity = 1)
sim.fn(quad.size = 10, cell.size = 1, intensity = 5)

# Increase scale (grain size) with constant density
set.seed(1)
sim.fn(quad.size = 20, cell.size = 1, intensity = 0.1)
sim.fn(quad.size = 20, cell.size = 2, intensity = 0.1)
sim.fn(quad.size = 20, cell.size = 2.5, intensity = 0.1)
sim.fn(quad.size = 20, cell.size = 5, intensity = 0.1)
sim.fn(quad.size = 20, cell.size = 10, intensity = 0.1)
sim.fn(quad.size = 20, cell.size = 20, intensity = 0.1)

# Can get occupancy equal to 1 for any rare occurring species when make
# cell size equal to study area size
sim.fn(quad.size = 20, cell.size = 20, intensity = 0.005)


# Inhomogeneous point pattern
# ---------------------------
# Increase density with constant spatial scale
set.seed(1)
str(dat <- simPPe(lscape.size = 150, buffer.width = 0, variance.X = 1, theta.X = 20,
       M = 10, beta = 1, quads.along.side = 6) )
str(dat <- simPPe(lscape.size = 150, buffer.width = 0, variance.X = 1, theta.X = 20,
       M = 50, beta = 1, quads.along.side = 6) )
str(dat <- simPPe(lscape.size = 150, buffer.width = 0, variance.X = 1, theta.X = 20,
       M = 100, beta = 1, quads.along.side = 6) )
str(dat <- simPPe(lscape.size = 150, buffer.width = 0, variance.X = 1, theta.X = 20,
       M = 200, beta = 1, quads.along.side = 6) )
str(dat <- simPPe(lscape.size = 150, buffer.width = 0, variance.X = 1, theta.X = 20,
       M = 300, beta = 1, quads.along.side = 6) )
str(dat <- simPPe(lscape.size = 150, buffer.width = 0, variance.X = 1, theta.X = 20,
       M = 500, beta = 1, quads.along.side = 6) )
str(dat <- simPPe(lscape.size = 150, buffer.width = 0, variance.X = 1, theta.X = 20,
       M = 1000, beta = 1, quads.along.side = 6) )
# Occupancy can be a proxy of abundance: perfect at very low abundance,
# useless once density so high that psi = 1


# Increase scale (grain size) with constant density
set.seed(1)
str(dat <- simPPe(lscape.size = 200, buffer.width = 0, variance.X = 1, theta.X = 20,
       M = 10, beta = 1, quads.along.side = 20) )
str(dat <- simPPe(lscape.size = 200, buffer.width = 0, variance.X = 1, theta.X = 20,
       M = 10, beta = 1, quads.along.side = 10) )
str(dat <- simPPe(lscape.size = 200, buffer.width = 0, variance.X = 1, theta.X = 20,
       M = 10, beta = 1, quads.along.side = 5) )
str(dat <- simPPe(lscape.size = 200, buffer.width = 0, variance.X = 1, theta.X = 20,
       M = 10, beta = 1, quads.along.side = 4) )
str(dat <- simPPe(lscape.size = 200, buffer.width = 0, variance.X = 1, theta.X = 20,
       M = 10, beta = 1, quads.along.side = 2) )
str(dat <- simPPe(lscape.size = 200, buffer.width = 0, variance.X = 1, theta.X = 20,
       M = 10, beta = 1, quads.along.side = 1) )











