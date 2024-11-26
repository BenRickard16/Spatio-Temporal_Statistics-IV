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












