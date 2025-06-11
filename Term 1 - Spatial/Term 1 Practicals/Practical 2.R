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



## 4. Specification and Fit of the Parametric Semivariogram Model

# Bit tricky but the main steps are:
# 1 Compute a sample semivariogram given data as previously
# 2 Specify the parametric form of the semivariogram model
# 3 calibrate/estimate/fit the specified parameter semivariogram against the computed sample semivariogram


## 4.1. Step 2: Specification of the Parametric Semivariogram Model

# Specifies parametric semivariogram model resulting by adding several standard semivariogram models
# 1. Specify the semivariogram model
vgm.tmp <- vgm(psill = 0.5, model = 'Nug', range = 0)

vgm.tmp <- vgm(psill = 1, model = 'Sph', range = 300, add.to = vgm.tmp)

vgm.tmp <- vgm(psill = 0.8, model = 'Sph', range = 800, add.to = vgm.tmp)

print(vgm.tmp)


## 4.2 Step : Fit the Parametric Semivariogram Model Against the Sample Semivariogram

# Fit the semivariogram model against the sample variogram via weighted least squares
library(sp)
library(sf)
# 0 Get the data as sf object
data(meuse)
meuse.df <- meuse
meuse.sf <- st_as_sf(meuse.df, 
                     coords = c("x", "y"),      # specify which columns are the coordinates
                     crs = 28992  # CRS code
)
# 1 Specify the semi variogram model
smpl.vgm <- variogram(object = log(zinc) ~ 1 + sqrt(dist), # formula associated to the trend
                      data = meuse.sf,           # data frame
)
vgm.obj <- vgm(psill = 0.5,
               model = "Nug", 
               range = 0
)
vgm.obj <- vgm(psill = 1,                   # the value "=1" is just an initial value
               model = "Sph",
               range = 300,
               add.to = vgm.obj
)
# 2 Fit the semi variogram model against the sample variogram
est.vgm <- fit.variogram(object = smpl.vgm, # sample variogram, output of variogram
                         model = vgm.obj    # variogram model, output of vgm
)
# 3 Print the estimated values for the parameters of the semi variogram model   
print(est.vgm)


## Task: consider the dataset in object jura.pred where Cadmium is quantity of interest. Want to fit a semivariogram nugget + exponential. Assume isotropy and the
# sample variogram with cutoff distance equal to 1.5 and width 0.1
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

# Step 1
smpl.vgm <- variogram(object = log(Cd) ~ 1 , # formula defining the response vector and (possible) regressors
                      data = jura.pred.sf,            # data frame
                      cutoff = 1.5,              # spatial separation distance up to which point pairs are included in semivariance estimates
                      width = 0.1                  # the width of subsequent distance intervals into which data point pairs are grouped for semivariance estimates
)
# Step 2
vgm.model <- vgm(psill = 0.1,
                 model = "Nug", 
                 range = 0
)
vgm.model <- vgm(psill = 0.5,
                 model = "Exp",
                 range = 1.5,
                 add.to = vgm.model
)
print(vgm.model)

# Step 3
est.vgm.model <- fit.variogram(object = smpl.vgm, # sample variogram, output of variogram
                               model = vgm.model    # variogram model, output of vgm
)

print(est.vgm.model)
plot(est.vgm.model, 
     cutoff = 1.9)



## 5. Spatial Prediction by Kriging

# Main steps for universal and simple Kriging:
# 1 Specify semivariogram model and fit against sample semivariogram estimate
# 2 Train geostatistical model
# 3 Compute predictive mean and variance at unseen points


## 5.1. Train Geostatistical Model

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

# Set the formula of the trend  
frml <- log(zinc) ~ 1 + sqrt(dist)

# Estimate the variogram  
smpl.vgm <- variogram(object = frml, 
                      data = meuse.sf)

par.vgm.0 <- vgm(psill = 0.5,
                 model = "Nug", 
                 range = 0)

par.vgm.1 <- vgm(psill = 1,
                 model = "Sph",
                 range = 300,
                 add.to = par.vgm.0)

par.vgm.est <- fit.variogram(object = smpl.vgm, 
                             model = par.vgm.1)

# Compute the Kriging equations
uni.krige.est <- gstat( formula = frml,
                        data = meuse.sf, 
                        model = par.vgm.est)

# Summary
summary(uni.krige.est)


## 5.2. Compute Predictive Mean and Variance

# Following script computes predictive mean and variance at unseen locations
# 0 Get the data => from previous script
# 1 Train the geostatistical model => from the previous script
# 2.1 Get the unseen locations  
data(meuse.grid)
meuse.grid.df <- meuse.grid
meuse.grid.sf <- st_as_sf(meuse.grid.df, 
                          coords = c("x", "y"),      # specify which columns are the coordinates
                          crs = 28992)  # CRS code

# 2.2 Compute the predictions at the unseen locations
uni.krige.prd <- predict(object = uni.krige.est, # object of class gstat; output of gstat()
                         newdata = meuse.grid.sf)   # data frame with prediction locations

# 2.3 Print summaries
print(uni.krige.prd)
summary(uni.krige.prd)


# 2.4 Plot predictions
library(mapview)
mapview(uni.krige.prd, zcol='var1.pred', layer.name = 'predicted mean')
mapview(uni.krige.prd, zcol='var1.var', layer.name = 'predicted variance')

# For dataset jura.pred, quantity of interest is cadmium. Deterministic trend is 0.1 (simple Kriging). Semivariogram
# has a nugget plus exponential. Compute predictive mean and variance at unseen locations in jura.grid

# Getting the data and libraries
library(sp)
library(sf)
library(gstat)
data(jura)
jura.pred.df <- as.data.frame(jura.pred)
jura.pred.sf <- st_as_sf(jura.pred.df, coords = c('long', 'lat'), crs=4326)

jura.grid.df <- as.data.frame(jura.grid)
jura.grid.sf <- st_as_sf(jura.grid.df, coords = c('long', 'lat'), crs=4326)

# Setting formula of the trend
frml <- log(Cd) ~ 1

# Estimate semivariogram
smpl.vgm <- variogram(object = frml, data=jura.pred.sf)

par.vgm.model <- vgm(psill = 0.1, model = 'Nug', range=0)
par.vgm.model <- vgm(psill=0.5, model = 'Exp', range=1.5, add.to = par.vgm.model)

par.vgm.est <- fit.variogram(object=smpl.vgm, model = par.vgm.model)

# Compute Kriging equations
sim.krige.est <- gstat(formula = frml, data = jura.pred.sf, model = par.vgm.est, 
                       beta=c(0.1)) # vector for simple kriging with trend coeffs

krige.prd <- predict(object = sim.krige.est, newdata = jura.grid.sf)

# Plotting
mapview(krige.prd, zcol='var1.pred', layer.name = 'predicted mean')
mapview(krige.prd, zcol='var1.var', layer.name = 'predicted variance')



## 6. Cross Validation

# The following script:
# 1 Trains a geostatistical model and saves training output in krige.fit
# 2 Calculates z-scores of a 5-fold cv
# 3 Plots the z-scores

library(sp)
library(sf)
# 0 Get the data as sf object
data(meuse)
meuse.df <- meuse
meuse.sf <- st_as_sf(meuse.df, 
                     coords = c("x", "y"),      # specify which columns are the coordinates
                     crs = 28992)  # CRS code

# 1 train the geostatistical model 
# Set the formula of the trend  
frml <- log(zinc) ~ 1 + sqrt(dist)
# Estimate the variogram  
smpl.vgm <- variogram(object = frml, 
                      data = meuse.sf)

par.vgm.0 <- vgm(psill = 0.5,
                 model = "Nug", 
                 range = 0)

par.vgm.1 <- vgm(psill = 1,
                 model = "Sph",
                 range = 300,
                 add.to = par.vgm.0)

par.vgm.est <- fit.variogram(object = smpl.vgm, 
                             model = par.vgm.1)

# Train  the model
krige.fit <- gstat( formula = frml,
                    data = meuse.sf, 
                    model = par.vgm.est)

# 2 Perform n-fold Cross Validation
nf.cv = gstat.cv(krige.fit,    # object of class gstat
                 nfold = 5)

# 3 Plot z-scores
bubble(obj = nf.cv['zscore'])


# Now we perform a 2-fold cv for the jura.pred Cadmium data with same deterministic trend and semviariogram parameterisation
# Load the data
library(sp)
library(sf)
library(gstat)
data(jura)
jura.pred.df <- as.data.frame(jura.pred)
jura.pred.sf <- st_as_sf(jura.pred.df, 
                         coords = c("long", "lat"),      # specify which columns are the coordinates
                         crs = 4326)  # CRS code

jura.grid.df <- as.data.frame(jura.grid)
jura.grid.sf <- st_as_sf(jura.grid.df, 
                         coords = c("long", "lat"),      # specify which columns are the coordinates
                         crs = 4326)  # CRS code

# Set the formula of the trend  
frml <- log(Cd) ~ 1

# Estimate the variogram  
smpl.vgm <- variogram(object = frml, 
                      data = jura.pred.sf)

par.vgm.model <- vgm(psill = 0.1,
                     model = "Nug", 
                     range = 0.0)

par.vgm.model <- vgm(psill = 0.5,
                     model = "Exp", 
                     range = 1.5,
                     add.to = par.vgm.model)

par.vgm.est <- fit.variogram(object = smpl.vgm, 
                             model = par.vgm.model)

# Compute the Kriging equations
sim.krige.est <- gstat( formula = frml,
                        data = jura.pred.sf, 
                        model = par.vgm.est,
                        beta = c(0.1)) #for simple kriging, vector with the trend coefficients (including intercept)  

# 2 Perform n-fold Cross Validation
nf.cv = gstat.cv(sim.krige.est,    # object of class gstat
                 nfold = 2)

# 3 Plot the z-scores
bubble(obj = nf.cv["zscore"])

# 4. I observe that at some locations the z-scores are way too large, e.g. larger than +/-3. Perhaps the fitting is not 
# the best, and the model has to be reconsidered by using different parametric semivariogram or by using different  
# different trend or by including additional covariates.  



## 7. Change of Support (Block Kriging)

# The following script computes and prints block Kriging prediction for blocks size 50x50 given spatial locations meuse.grid
library(sp)
library(sf)

# 1 Get data as sf object
data(meuse)
meuse.df <- meuse
meuse.sf <- st_as_sf(meuse.df, coords = c('x', 'y'), crs=28992)

# Set formula of trend
frml <- log(zinc) ~ 1 + sqrt(dist)

# Estimate the variogram
smpl.vgm <- variogram(object = frml, 
                      data = meuse.sf)

par.vgm.0 <- vgm(psill = 0.5,
                 model = "Nug", 
                 range = 0)

par.vgm.1 <- vgm(psill = 1,
                 model = "Sph",
                 range = 300,
                 add.to = par.vgm.0)

par.vgm.est <- fit.variogram(object = smpl.vgm, 
                             model = par.vgm.1)

# Compute Kriging equations
krige.est <- gstat(formula = frml, locations = meuse.sf, model = par.vgm.est)

# New locations
data(meuse.grid)
meuse.grid.df <- meuse.grid
meuse.grid.sf <- st_as_sf(meuse.grid.df, 
                          coords = c("x", "y"),      # specify which columns are the coordinates
                          crs = 28992)  # CRS code

# Kriging predictions
krige.prd <- predict(object = krige.est, newdata = meuse.grid.sf, block = c(50,50))

# Plotting
mapview(krige.prd, zcol='var1.pred', layer.name='predicted mean')
mapview(krige.prd, zcol='var1.var', layer.name='predicted variance')


# If we instead work with blocks of circular shape radius 20, centered on the points of meuse.grid
xy <- expand.grid(x = seq(-20, 20, 4), y = seq(-20, 20,4))
xy <- xy[(xy$x^2 + xy$y^2) <= 20^2, ]
krige.prd <- predict(object = krige.est,      # object of class gstat; output of gstat()
                     newdata = meuse.grid.sf, # data frame with prediction locations
                     block = xy)

# Plotting
mapview(krige.prd, zcol='var1.pred', layer.name='predicted mean')
mapview(krige.prd, zcol='var1.var', layer.name='predicted variance')



## 8. Multivariate Geostatistics

# Multiple dependent spatial variables and analysed jointly
# Need to organise the available information which is done sequentially using the function gstat
library(sp)
library(sf)
# 0 Get the data as sf object
data(meuse)
meuse.df <- meuse
meuse.sf <- st_as_sf(meuse.df,
                     coords = c("x", "y"),      # specify which columns are the coordinates
                     crs = 28992)  # CRS code

# 1 organize the multivariate response
gstat.obj <- NULL
gstat.obj <- gstat(g = gstat.obj,
                   id = "logCA",
                   formula = log(cadmium) ~ 1,
                   data = meuse.sf)

gstat.obj <- gstat(g = gstat.obj,
                   id = "logCO",
                   formula = log(copper) ~ 1,
                   data = meuse.sf)

gstat.obj <- gstat(g = gstat.obj,
                   id = "logLE",
                   formula = log(lead) ~ 1,
                   data = meuse.sf)

gstat.obj <- gstat(g = gstat.obj,
                   id = "logZI",
                   formula = log(zinc) ~ 1,
                   data = meuse.sf)

gstat.obj

# Computing and plotting sample semivariogram
smpl.svg <- variogram(object = gstat.obj)
plot(smpl.svg)

# Function fit.lmc fits a Linear Model of Coregionalisation to a multivariable sample variogram
vm.fit <- gstat::fit.lmc(v = smpl.svg,          # multivariate sample variogram
                         g = gstat.obj,         # gstat object, output of gstat
                         model = vgm(1, 'Sph', 800, 1))   # variogram model, output of vgm

# Plot sample and parametric cross semivariogram in same plot
plot(smpl.svg, vm.fit)

# Computing cokriging predictions of all the response variables at locations in meuse.grid and plots
data(meuse.grid)
meuse.grid.df <- meuse.grid
meuse.grid.sf <- st_as_sf(meuse.grid.df, coords = c('x', 'y'), crs = 28992)

cok.maps <- predict(vm.fit, meuse.grid.sf)

library(ggplot2)
library(viridisLite)

ggplot() +
  geom_sf(data = cok.maps, aes(color = logZI.pred)) +
  scale_color_viridis_c(name = 'log(ZI)') + theme_bw()

ggplot() + 
  geom_sf(data = cok.maps, aes(color = logZI.var)) +
  scale_color_viridis_c(name = "log(ZI) variance") + theme_bw()

ggplot() + 
  geom_sf(data = cok.maps, aes(color = logLE.pred)) +
  scale_color_viridis_c(name = "log(LE)") + theme_bw()

ggplot() + 
  geom_sf(data = cok.maps, aes(color = logLE.var)) +
  scale_color_viridis_c(name = "log(LE) variance") + theme_bw()

mapview(cok.maps,             # Object of class SpatialPixelsDataFrame
        zcol = "logZI.pred")      # vector of characters with the attribute name(s)
        
mapview(cok.maps,             # Object of class SpatialPixelsDataFrame
        zcol = "logLE.pred")      # vector of characters with the attribute name(s)
        
mapview(cok.maps,             # Object of class SpatialPixelsDataFrame
        zcol = "logZI.var")      # vector of characters with the attribute name(s)
        
mapview(cok.maps,             # Object of class SpatialPixelsDataFrame
        zcol = "logLE.var")      # vector of characters with the attribute name(s)



# From jura dataset coords long and lat, consider responses with order Ni, Pb, Zn and do
# 1 Compute multivariable sample variogram and plot
# 2 Fit a LMC to multivariable sample variogram and plot choosing any parametric semivariogram
# 3 Compute cokriging predictions of all the response variables at locations in jura.grid
library(sp)
library(sf)
library(gstat)
# 0 Get the data as sf object
data(jura)
jura.pred.df <- as.data.frame(jura.pred)
jura.pred.sf <- st_as_sf(jura.pred.df, 
                         coords = c("long", "lat"), 
                         crs = 4326)

# 0 organize the multivariate response
gstat.obj <- NULL
gstat.obj <- gstat(g = gstat.obj,
                   id = "logNi",
                   formula = log(Ni) ~ 1,
                   data = jura.pred.sf)

gstat.obj <- gstat(g = gstat.obj,
                   id = "logPb",
                   formula = log(Pb) ~ 1,
                   data = jura.pred.sf)

gstat.obj <- gstat(g = gstat.obj,
                   id = "logZn",
                   formula = log(Zn) ~ 1,
                   data = jura.pred.sf)

gstat.obj

# 1
smpl.svg <- variogram(object = gstat.obj)
plot(smpl.svg)

# 2
vm.fit <- gstat::fit.lmc(v = smpl.svg, 
                         g = gstat.obj, 
                         model = vgm(1, "Sph", 300, 1) )

plot(smpl.svg, vm.fit)

# 3
jura.grid.df <- jura.grid
jura.grid.sf <- st_as_sf(jura.grid.df,
                         coords = c("long", "lat"),  
                         crs = 4326 )

cok.maps <- predict( vm.fit, 
                     jura.grid.sf)

library(ggplot2)
library(viridisLite)
ggplot() + 
  geom_sf(data = cok.maps, aes(color = logNi.pred)) +
  scale_color_viridis_c(name = "log(Ni)") + theme_bw()

ggplot() + 
  geom_sf(data = cok.maps, aes(color = logNi.var)) +
  scale_color_viridis_c(name = "log(Ni) variance") + theme_bw()

ggplot() + 
  geom_sf(data = cok.maps, aes(color = logPb.pred)) +
  scale_color_viridis_c(name = "log(Pb)") + theme_bw()

ggplot() + 
  geom_sf(data = cok.maps, aes(color = logPb.var)) +
  scale_color_viridis_c(name = "log(Pb) variance") + theme_bw()

mapview(cok.maps,  
      zcol = "logNi.pred")

mapview(cok.maps,  
      zcol = "logPb.pred")  

mapview(cok.maps,  
      zcol = "logNi.var")  

mapview(cok.maps,  
      zcol = "logPb.var")  

