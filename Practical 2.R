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




