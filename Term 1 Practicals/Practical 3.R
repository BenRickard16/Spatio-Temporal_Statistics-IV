### Aerial Unit Data Modelling

## 1. Loading a DataSet

# Leukemia in New York dataset
library(spData)
library(sf)

# Read the data
NY8.data <- st_read(
  system.file("shapes/NY8_utm18.shp", package="spData"), 
  quiet=TRUE)

# Plot map of transformed proportions store in variable Z
library(tmap)
NY8.map <- NY8.data
tm_shape(NY8.map) + 
  tm_polygons("Z") 


# Columbus data plotting Crime response variable
library(spData)
library(sf)
# Read the data
columbus.data <- st_read(
  system.file("shapes/columbus.shp", package="spData"), 
  quiet=TRUE)

# Plot 
library(tmap)
columbus.map <- columbus.data
tm_shape(columbus.map) + 
  tm_polygons("CRIME") 



## 2. Spatial Neighbourhoods and Proximity Matrics

## 2.1. Neighbours List Based on Contiguity

# Neighbours based on contiguity assume neighbours share a common boundary
# The following script:
# 0 Loads the data
# 1 Creates a workinng copy of the dataset
# 2 Creates a neighbours list based on contiguity
# 3 Prints neighbours of first 3 locations in list of polygons
# 4 Plots borders and neighbour list as a network

# 0 load the data
library(spData)
library(sf)
nydata.data <- st_read(
  system.file("shapes/NY8_utm18.shp", package="spData"), 
  quiet=TRUE)

# 1 Create a working copy of the spatial data set
nydata.map <- nydata.data

# 2 Compute the neighbors list 
library(spdep)
nydata.nb <- spdep::poly2nb(pl = nydata.map, # list of polygons 
                            queen = TRUE)  # bounds

# 3 Print the neighbors of the first 3 locations in the list of polygones   
nydata.nb[1:3]  

# 4 plots the boarders and the neighbor list as a network 
# Plot the borders 
plot(st_geometry(nydata.map),        # object of class sf
     border = "black")             # color of polygon border(s); using NA hides them

#  Plot the neighbors (superimpose them to the previous plot) 
plot.nb(x = nydata.nb,                    # an object of class nb
        coords = st_geometry(nydata.map), # matrix of region point coordinates, a Spatial object (points or polygons),
        col = "red",
        add = TRUE)                     # superimpose them to the previous plot


# Now with the Columbus data we:
# 1 Create a list of neighbours by contiguity
# 2 Print the neighbours of first 3 locations in list of polygons
# 3 Plot the borders and neighbour list as a network

# 0 load the data
library(spData)
library(sf)
columbus.data <- st_read(
  system.file("shapes/columbus.shp", package="spData"), 
  quiet=TRUE)

# 1 Create a working copy of the spatial data set
columbus.map <- columbus.data

# 2 Compute the neighbors list 
library(spdep)
columbus.nb <- spdep::poly2nb(pl = columbus.map, # list of polygons 
                              queen = TRUE)  # bounds

# 3 Print the neighbors of the first 3 locations in the list of polygones   
columbus.nb[1:3] 

# 4 plots the boarders and the neighbor list as a network 
# Plot the borders 
plot(st_geometry(columbus.map),        # object of class sf
     border = "black")             # color of polygon border(s); using NA hides them

# Plot the neighbors (superimpose them to the previous plot) 
plot.nb(x = columbus.nb,                    # an object of class nb
        coords = st_geometry(columbus.map), # matrix of region point coordinates, a Spatial object (points or polygons),
        col = "red",
        add = TRUE)                     # superimpose them to the previous plot


## 2.2. Neighbours List Based on k Nearest Neighbours

# Constructed by assuming neighbours of a given are are the areas with smallest k centroid distances

# The following script:
# 1 creates a 4 nearest neighbors list
# 2 prints the neighbors of the first 3 locations in the list of polygons.
# 3 plots the boarders and the neighbor list as a network

# load the data
library(spData)
library(sf)
nydata.data <- st_read(
  system.file("shapes/NY8_utm18.shp", package="spData"), 
  quiet=TRUE)

# 1 creates a $k$ nearest neighbors list
nydata.map <- nydata.data

library(spdep)
# Get the centroids  
coo <- st_centroid(x = nydata.map)   # object of class sfg, sfc or sf

# Compute the neighbors list
nydata.nb <- knn2nb(
  knearneigh(coo,  #  centroids
             k = 4)) #  number nearest neighbors

# 2 print the neighbors of the first 3 locations
nydata.nb[1:3]

# 3 Plot the boarders and the neighbor list as a network
plot(st_geometry(nydata.map),        # object of class sf
     border = "black")             # color of polygon border(s); using NA hides them

plot.nb(x = nydata.nb,                    # an object of class nb
        coords = st_geometry(nydata.map), # matrix of region point coordinates, a Spatial object (points or polygons),
        col = "red",
        add = TRUE)                     # superimpose them to the previous plot


# Now using Columbus data
# 1 Create a k = 2 nearest neighbour list of neighbours
# 2 Print the neighbours of the first 3 locations in the list of polygons
# 3 Plot the borders and the neighbour list as a network

# Load the data
library(spData)
library(sf)
columbus.data <- st_read(
  system.file("shapes/columbus.shp", package="spData"), 
  quiet=TRUE)

# Create a working copy of the spatial data set
columbus.map <- columbus.data

# 1 Create a $k=2$ nearest neighbors list of the neighbors.
columbus.nb <- knn2nb(
  knearneigh(x = st_centroid(columbus.map),  #  centroids
             k = 2))                      #  number nearest neighbors

# 2. print the neighbors of the first $3$ locations in the list of polygons.
columbus.nb[1:3]

# 3. plot the boarders and the neighbor list as a network  
# Plot the borders 
plot(st_geometry(columbus.map),             # object of class sf
     border = "black")             # color of polygon border(s); using NA hides them

# Plot the neighbors (superimpose them to the previous plot) 
plot.nb(x = columbus.nb,                    # an object of class nb
        coords = st_geometry(columbus.map), # matrix of region point coordinates, a Spatial object (points or polygons),
        col = "red",
        add = TRUE)                     # superimpose them to the previous plot



## 2.3. Proximity Matrix / Spatial Weights for Neighbours Lists

# The following code:
# 1 Creates a neighbourhood list based on 4 nearest neighbours
# 2 Creates the proximity matrix as a basic binary coding

# load the data
library(spData)
library(sf)
nydata.data <- st_read(
  system.file("shapes/NY8_utm18.shp", package="spData"), 
  quiet=TRUE)

# Create a working copy of the spatial data set
nydata.map <- nydata.data

# 1 Create the list of neighbors based on 4 nearest neighbors  
library(spdep)
nydata.nb <- knn2nb(
  knearneigh(x = st_centroid(nydata.map),  #  centroids
             k = 4))                         #  number nearest neighbors
 
# 2 Create the proximity matrix  
nydata.nbw <- nb2listw(neighbours = nydata.nb, # list of neighbors 
                       style = "B")          # coding style    

# 3 Printinfs 
nydata.nbw$weights[1:3]


# Now using Columbus data we:
# 1 Create a neighbourhood structure based on 2 nearest neighbour
# 2 Create a proximity matrix as a row standardised coding
# 3 Print the proximity matrix

library(spData)
library(sf)
columbus.data <- st_read(
  system.file("shapes/columbus.shp", package="spData"), 
  quiet=TRUE)


# Create a working copy of the spatial data set
columbus.map <- columbus.data

# 1 Create the list of neighbors based on 2 nearest neighbors  
library(spdep)
columbus.nb <- knn2nb(
  knearneigh(x = st_centroid(columbus.map),  #  centroids
             k = 2))                         #  number nearest neighbors

# 2 Create the spatial neighborhood matrix  
columbus.nbw <- nb2listw(neighbours = columbus.nb, # list of neighbors 
                         style = "W")          # coding style    

# 3 Print 
columbus.nbw$weights[1:3]



## 3. Spatial Autocorellation Hypothesis Tests

## 3.1. Moran's Scatterplot

# Moran's I scatterplot visualises spatial autocorrelation in data by displaying observations of each are
# against its spatially lagged values
# Spatially lagged values for a given area is calculated as a weighted average of the neighbouring values for that area

# The following script loads the data, creates a neighborhood structure and creates a proximity matrix.
# Libraries
library(spdep)

# Get the data
NY8.data <- st_read(
  system.file("shapes/NY8_utm18.shp", package="spData")[1], 
  quiet=TRUE)

# Make a fresh copy of the data in a working object
NY8.map <- NY8.data

# Build a neighbours list based on regions with contiguous boundaries
NY8.nb <- poly2nb(pl = NY8.data) # list of polygons

# Supplements a neighbours list with spatial weights for the chosen coding scheme 
NY8.listw <- nb2listw(neighbours = NY8.nb, # an object of class nb 
                      style="W")            # style can take values “W”, “B”, “C”, “U”, “minmax” and “S”


# Following script computes the spatially lagged values for each of locations
# Compute the spatially lagged value's for each of the locations
lagCases = lag.listw(x = NY8.listw,       # a listw object created for example by nb2listw
                     var = NY8.data$Cases)  # a numeric vector the same length as the neighbours list in listw

NY8.map$lagCases <- lagCases

# Plot 
plot(x = NY8.map$lagCases,   # standardized values of lagged Cases
     y = NY8.map$Cases,      # standardized values of Cases
     main = "Moran’s I scatterplot",
     xlab = "variable",
     ylab = "spatialy lagged variable")

cor(NY8.map$lagCases, NY8.map$Cases)


# Use Columbus data and variable Crime to:
# 1 Produce Moran's I scatter plot
# 2 Is there evidence of spatial autocorrelation?

# load the data
library(spData)
library(sf)
columbus.data <- st_read(
  system.file("shapes/columbus.shp", package="spData"), 
  quiet=TRUE)

# Create a working copy of the spatial data set
columbus.map <- columbus.data

# Build a neighbours list based on regions with contiguous boundaries
columbus.nb <- poly2nb(pl = columbus.map) # list of polygons

# Supplements a neighbours list with spatial weights for the chosen coding scheme 
columbus.listw <- nb2listw(neighbours = columbus.nb, # an object of class nb 
                           style="W")            # style can take values “W”, “B”, “C”, “U”, “minmax” and “S”

# Compute the spatially lagged value's for each of the locations
lagCRIME = lag.listw(x = columbus.listw,       # a listw object created for example by nb2listw
                     var =columbus.data$CRIME)  # a numeric vector the same length as the neighbours list in listw

columbus.map$lagCRIME <- lagCRIME

# Plot 
plot(x = columbus.map$lagCRIME,   # standardized values of lagged CRIME
     y = columbus.map$CRIME,      # standardized values of CRIME
     main = "Moran’s I scatterplot",
     xlab = "variable",
     ylab = "spatialy lagged variable")

cor(columbus.map$lagCRIME, columbus.map$CRIME)

# Yes there is some evidence of spatial autocorrelation. I see a slope



## 3.2. Global Moran's I

# The following script loads the data and creates a proximity matrix

# Libraries
library(spdep)

# Get the data
NY8.data <- st_read(
  system.file("shapes/NY8_utm18.shp", package="spData")[1], 
  quiet=TRUE)

NY8.map <- NY8.data

# Build a neighbours list based on regions with contiguous boundaries
NY8.nb <- poly2nb(pl = NY8.data) # list of polygons

# Supplements a neighbors list with spatial weights for the chosen coding scheme (use style="B" for binary)
NY8.listw <- nb2listw(neighbours = NY8.nb, # an object of class nb 
                      style = "B")            # style can take values “W”, “B”, “C”, “U”, “minmax” and “S”


# The following
# 1 Tests null hypothesis no spatial autocorrelation, against the alternative hypothesis there is spatial autocorrelation
# 2 Prints summary of the result 

# 1
NY8.data.gmoran.test <- moran.test(NY8.data$Cases, # the data as a numeric vector the same length as the neighbours list in listw 
                                   NY8.listw,      # a listw object created for example by nb2listw
                                   alternative = "two.side") # one of "greater", "smaller",  "two.side"

# 2 
print(NY8.data.gmoran.test)


# Consider Columbus data variable Crime and test spatial autocorrelation against positive spatial autocorrelation at 5%
library(spData)
library(sf)
columbus.data <- st_read(
  system.file("shapes/columbus.shp", package="spData"), 
  quiet=TRUE)

columbus.map <- columbus.data
columbus.nb <- poly2nb(pl = columbus.data )
columbus.listw <- nb2listw(neighbours = columbus.nb, 
                           style = "B") 

# 1
columbus.data.gmoran.test <- moran.test(columbus.data$CRIME, # the data as a numeric vector the same length as the neighbours list in listw 
                                        columbus.listw, 
                                        alternative = "greater") 

# 2 
print(columbus.data.gmoran.test)

# There is significant evidence to reject the null hypothesis at 5% significance



## 3.3, Local Moran's I_i

# Following script loads the dataset and generates proximity matrix
# Libraries
library(spdep)
library(mapview)

# Get the data
NY8.data <- st_read(
  system.file("shapes/NY8_utm18.shp", package="spData")[1], 
  quiet=TRUE)

NY8.map <- NY8.data

# Build a neighbours list based on regions with contiguous boundaries
NY8.nb <- poly2nb(pl = NY8.data) # list of polygons

# supplements a neighbours list with spatial weights for the chosen coding scheme (use style="B" for binary)
NY8.listw <- nb2listw(neighbours = NY8.nb, # an object of class nb 
                      style = "B")            # style can take values “W”, “B”, “C”, “U”, “minimax” and “S”

# Following script computes local Moran's I_i and associated values
local.moran <- localmoran(x = NY8.data$Cases,   # the values of the variable of interest as a numeric vector the same length as the neighbours list in listw
                          listw = NY8.listw,         # a listw object created for example by nb2listw
                          alternative = "two.sided")  # one of "greater", "less" or "two.sided"

print(local.moran[1:10,])                       # print the first 10 locations

# Following script visualises the outputs by generating maps plotting local Moran's index and z-scores against the location by using functions of package tmap
library(tmap)  # load the library  

tmap_mode(
  mode = "plot"   # this is for printing the standard mode plot (not interactive)  
  # mode = "view"  # this is for printing the interactive plot  
  )             # set the mode to be printable (not the interactive)

# Create a fresh copy of the spatial dataset  
NY8.map <- NY8.data  

# Variable of interest 
# tm_shape creates a tmap-element that specifies a spatial data object, which we refer to as shape.
p1 <- tm_shape(shp = NY8.map)    # shape object
                 
# tm_polygons fills the polygons and draws the polygon borders.
p1 <- p1 + tm_polygons(col = "Cases", # the name of a data variable 
                       title = "vble",         # title 
                       style = "hclust")        # method to process the color scale when col is a numeric variable, here by clustering. 
                  
# tm_layout specifies the map layout; e.g. controls title, margins, aspect ratio, colors, frame, legend, among many other things.
p1 <- p1 +  tm_layout(legend.outside = TRUE)  # determines whether the legend is plot outside of the map/facets. 

# Local Moran's I   

# Append the local Moran's I in the sf data.frame
NY8.map$lmI <- local.moran[, "Ii"] 

# tm_shape creates a tmap-element that specifies a spatial data object, which we refer to as shape.
p2 <- tm_shape(shp = NY8.map)    # shape object
                 
# tm_polygons fills the polygons and draws the polygon borders.
p2 <- p2 + tm_polygons(col = "lmI",       # the name of a data variable 
                       title = "Local Moran's I",  # title 
                       style = "hclust")            # method to process the color scale when col is a numeric variable, here by clustering. 
                  
# tm_layout specifies the map layout; e.g. controls title, margins, aspect ratio, colors, frame, legend, among many other things.
p2 <- p2 +  tm_layout(legend.outside = TRUE)  # determines whether the legend is plot outside of the map/facets. 
          
# Local Moran's I Z-score    

# Append the local Moran's I in the sf data.frame
NY8.map$lmZ <- local.moran[, "Z.Ii"] 

# tm_shape creates a tmap-element that specifies a spatial data object, which we refer to as shape.
p3 <- tm_shape(shp = NY8.map)    # shape object
                 
# tm_polygons fills the polygons and draws the polygon borders.
p3 <- p3 + tm_polygons(col = "lmZ",       # the name of a data variable 
                       title = "Z-score",          # title 
                       breaks = c(-Inf, 1.65, Inf)) # method to process the color scale when col is a numeric variable, here by setting the breaks manually 
                  
# tm_layout specifies the map layout; e.g. controls title, margins, aspect ratio, colors, frame, legend, among many other things.
p3 <- p3 +  tm_layout(legend.outside = TRUE)  # determines whether the legend is plot outside of the map/facets. 
          
# Arrange small multiples in a grid layout. 
tmap_arrange( p1, p2, p3)


# Consider variable crime from dataset Columbus. Compute local Moran's I_is and associated values, printing the output
# Then visualise the outputs by generating maps plotting local Moran's index and expectation against the location by using tmap
library(spData)
library(sf)
columbus.data <- st_read(
  system.file("shapes/columbus.shp", package="spData"), 
  quiet=TRUE)

columbus.map <- columbus.data
columbus.nb <- poly2nb(pl = columbus.data )
columbus.listw <- nb2listw(neighbours = columbus.nb, 
                           style = "B")  


local.moran <- localmoran(x = columbus.data$CRIME,   # the values of the variable of interest as a numeric vector the same length as the neighbours list in listw
                          listw = columbus.listw,         # a listw object created for example by nb2listw
                          alternative = "two.sided"  # one of "greater", "less" or "two.sided"
)
print(local.moran[1:10,])                       # print the first 10 locations

library(tmap)  # load the library  

tmap_mode(
  mode = "plot"   # this is for printing the standard mode plot (not interactive)  
  #  mode = "view"  # this is for printing the interactive plot  
)               # set the mode to be printable (not the interactive)

# Create a fresh copy of the spatial dataset  

columbus.map <- columbus.data  

# variable of interest ####################  

## tm_shape creates a tmap-element that specifies a spatial data object, which we refer to as shape.
p1 <- tm_shape(shp = columbus.map    # shape object
)                 
## tm_polygons fills the polygons and draws the polygon borders.
p1 <- p1 + tm_polygons(col = "CRIME", # the name of a data variable 
                       title = "vble",         # title 
                       style = "hclust"        # method to process the color scale when col is a numeric variable, here by clustering. 
)                  
## tm_layout specifies the map layout; e.g. controls title, margins, aspect ratio, colors, frame, legend, among many other things.
p1 <- p1 +  tm_layout(legend.outside = TRUE  # determines whether the legend is plot outside of the map/facets. 
)
# local Moran's I   ##########  

## append the local Moran's I in the sf data.frame
columbus.map$lmI <- local.moran[, "Ii"] 
## tm_shape creates a tmap-element that specifies a spatial data object, which we refer to as shape.
p2 <- tm_shape(shp = columbus.map    # shape object
)                 
## tm_polygons fills the polygons and draws the polygon borders.
p2 <- p2 + tm_polygons(col = "lmI",       # the name of a data variable 
                       title = "Local Moran's I",  # title 
                       style = "hclust"            # method to process the color scale when col is a numeric variable, here by clustering. 
)                  
## tm_layout specifies the map layout; e.g. controls title, margins, aspect ratio, colors, frame, legend, among many other things.
p2 <- p2 +  tm_layout(legend.outside = TRUE  # determines whether the legend is plot outside of the map/facets. 
)          
# local Moran's I Z-score   ##########  

## append the local Moran's I in the sf data.frame
columbus.map$lmZ <- local.moran[, "Z.Ii"] 
## tm_shape creates a tmap-element that specifies a spatial data object, which we refer to as shape.
p3 <- tm_shape(shp = columbus.map    # shape object
)                 
## tm_polygons fills the polygons and draws the polygon borders.
p3 <- p3 + tm_polygons(col = "lmZ",       # the name of a data variable 
                       title = "Z-score",          # title 
                       breaks = c(-Inf, 1.65, Inf) # method to process the color scale when col is a numeric variable, here by setting the breaks manually 
)                  
## tm_layout specifies the map layout; e.g. controls title, margins, aspect ratio, colors, frame, legend, among many other things.
p3 <- p3 +  tm_layout(legend.outside = TRUE  # determines whether the legend is plot outside of the map/facets. 
)          
## Arrange small multiples in a grid layout. 
png("./aa.png")
tmap_arrange( p1, p2, p3)
dev.off()



## 4. Fitting AutoRegressive Models

## 4.1. SAR

# Following script loads the data and computes proximity matrix
library(spData)

# Get the data
NY8.data <- st_read(
  system.file("shapes/NY8_utm18.shp", package="spData")[1], 
  quiet=TRUE)

# Create a working data set map  
NY8.map <- NY8.data

# Build a neighbours list based on regions with contiguous boundaries
NY8.nb <- poly2nb(pl = NY8.data)

# Supplements a neighbours list with spatial weights for the chosen coding scheme (use style="B" for binary)
NY8.listw <- nb2listw(neighbours = NY8.nb, # an object of class nb 
                      style = "B")           # style can take values “W”, “B”, “C”, “U”, “minmax” and “S”


# Following R script fits a SAR model and prints its summary
library(spatialreg)
# 1
NY.sar.fit <- spautolm(formula = Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, #a symbolic description of the model to be fit. The details of model specification are given for lm()  
                       data = NY8.data,   # an optional data frame containing the variables in the model.  
                       listw = NY8.listw, # a listw object created for example by nb2listw
                       family = "SAR")     # character string: either "SAR" or "CAR"

# 2
summary(NY.sar.fit)

# It seems that there is significant spatial correlation because the estimated value of λ is 0.0357302 and the p-value of the likelihood ratio test for λ=0
# is below 0.05

# The proximity to a TCE seems not to be significant, although its p-value is close to being significant at the 95 % level and it would be advisable not 
# to discard a possible association and to conduct further research on this. The other two covariates are significant, suggesting that census tracts with 
# larger percentages of older people and with lower percentages of house owners have higher transformed incidence rates.


# The following R script
# 1 Maps the estimated signal trends against location
# 2 Maps the estimated stochastic signal against the locations

# 1
NY8.data$signal_trend <- NY.sar.fit$fit$signal_trend
mapview(NY8.data, 
        zcol = "signal_trend")

# 2
NY8.data$signal_stochastic <- NY.sar.fit$fit$signal_stochastic
mapview(NY8.data, 
        zcol = "signal_stochastic")


# Now consider Columbus OH data, first loading the data and generating a proximity matrix
columbus.data <- sf::st_read(system.file("shapes/columbus.shp", package="spData")[1])

columbus.map <- columbus.data
library(spdep)
columbus.np <- spdep::poly2nb(pl = columbus.map, # list of polygons 
                              queen = TRUE)  # bounds  

columbus.listw <- nb2listw(columbus.np, 
                           style = "W")

# We now do the following
# 1 Fit a SAR model considering response to be crime and signal trend constant mu = 1
# 2 Is the spatial autocorrelation significant?
# 3 Map the estimated stochastic signal against the locations

# 1 
columbus.sar.fit <- spautolm(formula = CRIME ~ 1, #a symbolic description of the model to be fit. The details of model specification are given for lm()  
                             data = columbus.data,   # an optional data frame containing the variables in the model.  
                             listw = columbus.listw, # a listw object created for example by nb2listw
                             family = "SAR")     # character string: either "SAR" or "CAR"

# 2
summary(columbus.sar.fit)

# Yes the spatial autocorrelation is statistically significant at sig. level 5% because the p-value for lambda is 4.4485e-06 

# 3
library(mapview)
columbus.map$signal_stochastic <- columbus.sar.fit$fit$signal_stochastic
mapview(columbus.map, 
        zcol = "signal_stochastic")



## 4.2. CAR

# Get the data
NY8.data <- st_read(
  system.file("shapes/NY8_utm18.shp", package="spData")[1], 
  quiet=TRUE)

NY8.map <- NY8.data 

# Build a neighbours list based on regions with contiguous boundaries
NY8.nb <- poly2nb(pl = NY8.data)

# supplements a neighbours list with spatial weights for the chosen coding scheme (use style="B" for binary)
NY8.listw <- nb2listw(neighbours = NY8.nb, # an object of class nb 
                      style = "B")           # style can take values “W”, “B”, “C”, “U”, “minmax” and “S”

# We do the following:
# 1 Fit a CAR model with regressors in linear signal PEXPOSURE, PCTAGE65P and PCTOWNHOME spatial neighbourhood matrix in NY8.listw, and weights k_i
# are column POP8 and name output NY.car.fit
# 2 Print summary of fitting

# 1
NY.car.fit <- spautolm(formula = Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, #a symbolic description of the model to be fit. The details of model specification are given for lm()  
                       data = NY8.data,        # an optional data frame containing the variables in the model.  
                       listw = NY8.listw,   # a listw object created for example by nb2listw
                       family = "CAR",     # character string: either "SAR" or "CAR"
                       weights = NY8.data$POP8) # a vector of weights to be used in the fitting process

# 2. 
summary(NY.car.fit)





































