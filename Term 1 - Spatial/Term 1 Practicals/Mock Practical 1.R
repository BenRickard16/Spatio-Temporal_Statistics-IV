## Spatial Data

# Question 2
# Loading necessary libraries
library(sf)

# Create the data frame
data.df <- data.frame(
  ID = 1:5,
  Longitude = c(-73.935242, -74.006015, -73.994484, -74.009785, -73.985130),
  Latitude = c(40.730610, 40.712776, 40.740129, 40.719296, 40.758896),
  Value = c(10, 20, 15, 25, 30))

# Convert to a sf object
data.sf <- st_as_sf(data.df, coords = c('Longitude', 'Latitude'), crs = 4326)

print(data.sf)


# Plotting the 'Value' on a map
data.map <- data.sf

library(tmap)

# Create and print map
tmap_mode('view')
map <- tm_shape(data.map) + 
  tm_dots(col = "Value", 
          palette = "Set1", 
          size = 0.5, 
          title = "Location ID") + 
  tm_layout(title = "Map of sf Object with CRS 4326")


# Question 3
# Create sf object with 2 points and their attribute values, and a polygon centered on specified point, plotted with ggplot2
# Single point (point as a vector)
p1_sfg <- st_point(c(1,4))     # a point
p2_sfg <- st_point(c(2,4))   # a point

# Polygon
p1 <- rbind(c(0, 0), c(0, 2), c(2, 2), c(2,0), c(0,0) )  
pol_sfg <- st_polygon(list(p1))

# Create sf object
p_sfc <- st_sfc(p1_sfg, 
                p2_sfg, 
                pol_sfg)

df <- data.frame( v1 = c("tree 1", "tree 2", "park") ) # the attributes  
p_sf <- st_sf(df, 
              geometry = p_sfc)

# Question 4
# Plot them (to be introduced later in the computer practical)
library(ggplot2)
ggplot(p_sf) + 
  geom_sf(aes( col = v1 ), size = 3) + 
  theme_bw()

# Calculate number of points in each unit of the map
library(spData)
library(sf)

gk.map <- st_read(
  system.file("shapes/columbus.shp", package="spData"), 
  quiet=TRUE)

gk.map <- gk.map[1:4,]

gk.map$Division <- c("A","B","C","D")

gk.map <- gk.map[,c("Division")]

set.seed(2024)

gk.points <- st_sample(gk.map,
                       size = 30)

# Intersection (first argument map, then points)
gk.inter <- st_intersects(gk.map, gk.points )
# Add point count to each polygon
gk.map$count <- lengths(gk.inter)
print(gk.map$count)



## Geostatistics

# Importing data
rm(list=ls())
scallop.mod <- read.csv('https://raw.githubusercontent.com/georgios-stats/Spatio-Temporal_Statistics_Michaelmas_2024/refs/heads/main/Misc/scallop.mod.csv',sep=',')
print(scallop.mod)
scallop.mod.df <- as.data.frame(scallop.mod)
library(sp)
library(sf)

# Question 7
scallop.mod.map <- st_as_sf(scallop.mod.df, 
                            coords = c("longitude", "latitude"),
                            crs = 4326)

scallop.mod.map$longitude <-  scallop.mod.df$longitude
scallop.mod.map$latitude <-  scallop.mod.df$latitude

# Compute the sample variogram cloud
library(gstat)
frml <- log.1.p.tot.catch ~ 1 
scallop.mod.map.cld <- variogram(object = frml,
                                 data = scallop.mod.map,
                                 cloud = TRUE)

# Plot cloud semivariogram
plot(scallop.mod.map.cld,
     xlab = 'distance', 
     ylab = 'gamma', 
     main = 'semivariogram cloud')

# Question 9
# Fit parametric semivariogram spherical with nugget and plot with cutoff 50
# 1 Specify the semi variogram model
library(gstat)
frml <- log.1.p.tot.catch ~ 1 
smpl.vgm <- variogram(object = frml,
                      data = scallop.mod.map)

vgm.obj <- vgm(psill = 0.5,
               model = "Nug", 
               range = 0)

vgm.obj <- vgm(psill = 1,
               model = "Sph",
               range = 100,
               add.to = vgm.obj)

# 2 Fit the semi variogram model against the sample variogram
est.vgm <- fit.variogram(object = smpl.vgm,
                         model = vgm.obj)

# 3 Print the estimated values for the parameters of the semi variogram model   
print(est.vgm)

# 4 plot it 
plot(est.vgm, 
     cutoff = 50)

# Question 14
# Use Kriging to estimate value and variance at unseen location (-73.7, 38.6)
krige.fit <- gstat( formula = frml,
                    data = scallop.mod.map, 
                    model = est.vgm)

scallop.mod.new.df <- data.frame(longitude = c(-73.7), 
                                latitude = c(38.6) )

scallop.mod.new.sf <- st_as_sf(scallop.mod.new.df, 
                               coords = c("longitude", "latitude"),
                               crs = 4326)

krige.prd <- predict(object = krige.fit,
                     newdata = scallop.mod.new.sf)

krige.prd



## Aerial Data Unit

# Load in SIDs data
library(sf)

nc.data <- st_read(system.file("shapes/sids.gpkg", package="spData")[1], quiet=TRUE)
nc.data.sf <- nc.data

# Question 19
# Hypothesis test to infer spatial autocorrelation for number of SIDs 5% standardised row proximity matrix and neighbours based on contiguous regions

# Build neighbours list
library(spdep)

nc.nb <- poly2nb(pl = nc.data.sf )

# Supplements a neighbors list with spatial weights for the chosen coding scheme (use style="B" for binary)
nc.listw <- nb2listw(neighbours = nc.nb, 
                     style = "W")

# Perform Global Moran’s I hypothesis test 
nc.data.gmoran.test <- moran.test(nc.data.sf$SID74,
                                  nc.listw, 
                                  alternative = "two.side")

# Print
print(nc.data.gmoran.test)


# Question 22
# Fit a SAR model
nc.data.sf$logSID79 <- log(nc.data.sf$SID79+1.0)
nc.data.sf$logBIR79 <- log(nc.data.sf$BIR79)

# Supplements a neighbors list with spatial weights for the chosen coding scheme (use style="B" for binary)
nc.listw <- nb2listw(neighbours = nc.nb, 
                     style = "W") 

# Fit the SAR model
library(spatialreg)

nc.sar.fit <- spautolm(formula = logSID79 ~ 1 + logBIR79,
                       data = nc.data.sf,
                       listw = nc.listw, 
                       family = "SAR") 

# Print the summary
summary(nc.sar.fit)

# Question 26
# Plot of estimated signal trend
nc.data$signal_trend <- nc.sar.fit$fit$signal_trend

library(tmap)
tmap_mode("plot")
tm_shape(nc.data) + 
  tm_polygons("signal_trend") 

# Insert the fitted value to site CNTY_ID 1825, NAME Ashe to the data set
nc.sar.fit$fit$fitted.values[1]

