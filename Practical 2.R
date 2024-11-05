### Point Referenced Spatial Data Modelling ###

## 1. Simulate a Random Field
library(geoR)

# Simulate n=100 random points using Exponential covariance functions with parameter 1 for partial sill
# and 0.25 for the range, and no nugget

sim.grf <- grf(n = 100,              # number of points (spatial locations) in each simulations.
               cov.model = "exponential", # exp cov function family
               cov.pars = c(1, .25), # a vector with 2 elements with values of the covariance 
               # parameters \sigma^2(partial sill) and \phi (range parameter)
               nugget = 0.0,         # the value of the nugget effect parameter
               xlims = c(0, 1), 
               ylims = c(0, 1) 
)

# Creating an sf object of these points
sim.grf.df <- data.frame(x = sim.grf$coords[,1],
                           y = sim.grf$coords[,2],
                           value = sim.grf$data)

library(sf)
sim.grf.sf <- st_as_sf(sim.grf.df, coords = c('x', 'y'))

# Plotting the points
library(ggplot2)
ggplot(data=sim.grf.sf) + 
  geom_sf(aes(color=value)) +
  theme_minimal() + 
  labs(y='y', x='x', title='GRF')


# Simulating GRF with 500 poiints, Gaussian covariance function partial sill 2 and range 1.5, and a
# nugget 0.1
library(geoR)
sim.grf <- grf(n = 500,              # number of points (spatial locations) in each simulations.
               cov.model = "gaussian", # Gaussian cov function family
               cov.pars = c(2, 1.5), # a vector with 2 elements with values of the covariance 
               # parameters \sigma^2(partial sill) and \phi (range parameter)
               nugget = 0.1,         # the value of the nugget effect parameter
               xlims = c(0, 1)
)

# Create an sf object
sim.grf.df <- data.frame(x = sim.grf$coords[,1],
                         y = sim.grf$coords[,2],
                         value = sim.grf$data)

library(sf)
sim.grf.sf <- st_as_sf( sim.grf.df, 
                        coords = c("x", "y"))

# Plot
library(ggplot2)
ggplot(data = sim.grf.sf) + 
  geom_sf(aes(color = value)) +
  theme_minimal() +
  labs(y= "y", x = "x", title ="GRF")


## 3.1 Semivariogram Cloud

# Computing and plotting the semivariogram cloud corresponding to delta(), a zero-mean intrinsic random field
library(sp)
library(sf)

# 1 Get the data as sf object
data(meuse)
meuse.df <- meuse
meuse.sf <- st_as_sf(meuse.df, 
                     coords = c("x", "y"),      # specify which columns are the coordinates
                     crs = 28992)  # CRS code

# 2 Compute the sample variogram cloud
library(gstat)
meuse.sv.cld <- variogram(object = log(zinc) ~ 1 + sqrt(dist), # formula defining the response vector and (possible) regressors
                          data = meuse.sf,                        # data frame
                          cloud = TRUE)                         # if TRUE, it calculates the semivariogram cloud

# ...with some investigation
names(meuse.sv.cld)
head(meuse.sv.cld)

# 3 Plot
plot(meuse.sv.cld,
     xlab='distance',
     ylab='gamma',
     main='semivariogram cloud')

plot(x = meuse.sv.cld$dist,
    y = meuse.sv.cld$gamma,
    cex = .5, 
    pch = 3, 
    col = 'blue', 
    xlab = 'distance', 
    ylab = 'gamma', 
    main = 'semivariogram cloud')

# Consider new object jura.pred from jura{gstat}. The quantity of interest is cadmium Cd in log-scale. Consider a 
# constant deterministic drift function mu = beta_0 unknown. Compute and plot the semivariogram cloud.

# First organising the data as an sf object
data(jura)
jura.pred.df <- as.data.frame(jura.pred)
jura.pred.sf <- st_as_sf(jura.pred.df,
                         coords = c('long', 'lat'),
                         crs=4326)

# Compute sample variogram cloud
jura.pred.sv.cld <- variogram(object = log(Cd) ~ 1,     # formula defining the response vector and regressors
                              data = jura.pred.sf,      # data frame
                              cloud = TRUE)             #  TRUE to calculate semivariogram cloud

# Investigation
names(jura.pred.sv.cld)
head(jura.pred.sv.cld)

# Plotting the cloud
plot(jura.pred.sv.cld,
     xlab = 'distance',
     ylab = 'gamme',
     main = 'semivariogram cloud')


## 3.2 Isotropic sample semivariogram (non-parametric estimate)

# Compute and plots non-parametric Matheron's semi-variogram estimate corresponding to delta()
rm(list = ls())
library(sp)
library(sf)

# 1 Get the data as sf object
data(meuse)
meuse.df <- meuse
meuse.sf <- st_as_sf(meuse.df, 
                     coords = c("x", "y"),      # specify which columns are the coordinates
                     crs = 28992  # CRS code
)

# Compute sample semivariogram cloud
svgm <- variogram(object = log(zinc) ~ 1 + sqrt(dist), # formula defining the response vector and (possible) regressors
                  data = meuse.sf,                        # data frame
                  cloud=FALSE
)

summary(svgm)

# Plot
plot(x = svgm$dist,
     y = svgm$gamma,
     cex = .5,
     pch = 3,
     col = 'blue',
     xlab = 'distance',
     ylab = 'gamma',
     main = 'sample variogram')

# Add labels
text(svgm$dist, 
     svgm$gamma, 
     svgm$np, 
     adj = c(0,2))


## 3.3 Anisotropic Sample Semivariogram (Non-Parametric Estimate)

# Function variogram computes non-parametric Matheron's estimate for specific angles by the argument 'alpha'
# Compute semivariogram estimate for angle 45 degrees, cutoff distance 1000, and width 50

rm(list = ls())
library(sp)
library(sf)

# 1 Get the data as sf object
data(meuse)
meuse.df <- meuse
meuse.sf <- st_as_sf(meuse.df, 
                     coords = c("x", "y"),      # specify which columns are the coordinates
                     crs = 28992)  # CRS code

# 2 Compute the semivariogram estimate
svgm.tmp <- variogram(object = log(zinc) ~ 1 + sqrt(dist), # formula defining the response vector and (possible) regressors
                      data = meuse.sf,               # data frame
                      alpha = c(45),              # direction in plane (x,y),
                      cutoff = 1000,              # spatial separation distance up to which point pairs are included in semivariance estimates
                      width = 50)                  # the width of subsequent distance intervals into which data point pairs are grouped for semivariance estimates

# Plot
plot( svgm.tmp )


# Consider jura.pred dataset again with quantity of interest Cd. Consider deterministic drift function mu = b_0
# unknown. Compute and plot sample semivariogram cloud angle 0,45,60,90 degrees, cutoff distance 1.5, width 0.1
rm(list=ls())
# Load the data
library(sp)
library(sf)
library(gstat)

data(jura)
jura.pred.df <- as.data.frame(jura.pred)
jura.pred.sf <- st_as_sf(jura.pred.df, 
                         coords = c("long", "lat"),      # specify which columns are the coordinates
                         crs = 4326  # CRS code
)

# 2 Compute the semivariogram estimate
svgm.tmp <- variogram( object = log(Cd) ~ 1 , # formula defining the response vector and (possible) regressors
                       data = jura.pred.sf,                # data frame
                       alpha = c(0,45,60,90),
                       cutoff = 1.5,
                       width = 0.1)

plot( svgm.tmp )

















