## 1. Introduction

# Installing required packages
install.packages(pkgs = c( "st", "sf", "mapview", "tmap", "ggplot2", "spData", "spdep", "SemiPar"))



## 2. Vector Data

# Check out the default R datasets availabke
path.to.shape.files <- system.file("shape", package='sf')
            
print(dir(path.to.shape.files))



## 4. Reading Vector Data Via SF

# Library sf
library(sf)

# st_read to read the spatial dataset nc{sf}
path.to.nc.shp <- system.file("shape/nc.shp", package = "sf")

map <- st_read(path.to.nc.shp, quiet=TRUE)

class(map)
head(map)


# Display the geometry of an sf object
rm( list = ls() )
library(sf)
path.to.nc.shp <- system.file("shape/nc.shp", 
                              package = "sf")
nc.map <- st_read(path.to.nc.shp, 
                  quiet = TRUE)
nc.geometry <- st_geometry(nc.map)
print(nc.geometry)

# TASK: change value of the variable NAME of 30th observation from Durham to New_Durham
rm( list = ls() )
library(sf)
path.to.nc.shp <- system.file("shape/nc.shp", 
                              package = "sf")
nc.map.1 <- st_read(path.to.nc.shp, 
                    quiet = TRUE)

nc.map.1[30,]$NAME <- 'New_Durham'
print(nc.map.1[30,'NAME'])



## 5. Making Maps with SF Objects

# First we use ggplot2
# Load spatial data
library(sf)
path.to.nc.shp <- system.file("shape/nc.shp", 
                              package = "sf")
nc.map <- st_read(path.to.nc.shp, 
                  quiet = TRUE)

# Load package
library(ggplot2)

# Create a map for the attribute SID74
ggplot(nc.map) +                   # sf object holding the spatial data to be plotted
  geom_sf(aes(fill = SID74)) +     # the name of the variable I wish to map as colour  
  ggtitle("North Carolina SIDS data") +
  xlab("x") +
  ylab("y")

# TASK: Map attribute lead from dataset meuse {sp}
rm( list = ls() )

# Load spatial data
library(sp)
data(meuse)
meuse <- st_as_sf(meuse, 
                  coords = c("x", "y"),      # specify which columns are the coordinates
                  crs = 28992)                # CRS code
meuse.map <- meuse
# Load package
library(ggplot2)

# Create the map for the attribute lead
ggplot(meuse.map) +                   # sf object holding the spatial data to be plotted
  geom_sf(aes(fill = lead)) +
  ggtitle("Meuse river data set") +
  xlab("x") +
  ylab("y")

# Now we use the mapview package and function
rm( list = ls() )
# Load spatial data
library(sf)
path.to.nc.shp <- system.file("shape/nc.shp", 
                              package = "sf")
nc.map <- st_read(path.to.nc.shp, 
                  quiet = TRUE)

# Load package
library(mapview)

# Create an Interactive Map for the attribute SID74
mapview(nc.map,        # sf object holding the spatial data to be plotted 
        zcol = "SID74", # attribute name(s) or column number(s) in attribute table of the column(s) to be rendered.
        map.types = "CartoDB.Positron") 

# TASK: Draw a basic interactuve map for attribute lead from meuse dataset
rm(list=ls())

# Load spatial data
library(sp)
data(meuse)
meuse <- st_as_sf(meuse, coords=c("x", "y")     # specify which columns are coordinates
                  , crs=28992)        # CRS code

# Load package
library(mapview)

# Creating the map for attribute lead
mapview(meuse,                        # sf object holding the spatial data
        zcol = 'lead',                # attribute name to be plotted
        map.types = 'OpenTopoMap')


# Now we do the same using library tmap
library(tmap)

# Get the spatial data
library(sf)
pathshp <- system.file('shape/nc.shp', package='sf')
nc.data <- st_read(pathshp, quiet=TRUE)
nc.map <- nc.data

# Set the mode, and create the thematic map
tmap_mode('plot')
tm_shape(nc.map) +
  tm_polygons('SID74')

# Do the same for interactive map for attribute lead in meuse
rm(list=ls())

# Accessing the libraries and data
library(sp)
library(sf)
library(gstat)

data(meuse)
meuse <- st_as_sf(meuse, coords = c('x', 'y'), crs = 28992)

# Creating the map
meuse.map <- meuse
library(tmap)

tmap_mode('view')

tm_shape(meuse.map) + 
  tm_dots('lead') +
  tm_basemap('CyclOSM')



## 6. Transforming Point Data to SF Object

# We wish to create a sf object from a data.frame containing the coords of locations and attributes
rm(list=ls())
library(mapview)
# Create the dataframe
df.data <- data.frame(
  place = c("London", "Paris", "Madrid", "Rome"),
  long = c(-0.118092, 2.349014, -3.703339, 12.496366),
  lat = c(51.509865, 48.864716, 40.416729, 41.902782),
  value = c(200, 300, 400, 500))
class(df.data)

# Use function st_as_sf() to turn data.frame into a sf object
sf.data <- st_as_sf(x = df.data,          #  enter data
                    coords = c("long", "lat"), # specify the coordinates
                    crs = 4326)

# Use function st_crs() to set the CRS to represent long and lat
st_crs(sf.data) <- 4326

# Plot
mapview(sf.data)


# Transform the following data.frame into sf object. Then produce a map of the variable logtotalCatch by using mapview.
library(HRW)  
data(scallop)
scallop$logtotalCatch <- log(scallop$totalCatch+1.0)

# Load the library
library(sf)
df.data <- scallop

# Transform data.frame into sf object 
sf.data <- st_as_sf( x = df.data,                     #  enter data
                     coords = c("longitude", "latitude"), # specify the coordinates
                     crs = 4326                           # set the CRS code
)

# Plot
mapview(sf.data, 
        zcol = "logtotalCatch")



## 7. Creating sf Objects (More Advanced)

# A sf object contains objects of class sf, sfc, and sfg
# sfg objects can be of type point, multipoing, or polygon

# General steps:
# 1. Create simple feature geometry objects sfg of type point, multipoint and polygon
# 2. Use st_sfc() to create a simple feature geometry list-column sfc with the sfg objects
# 3. Use st_sf() to put the data.frame with the attribute and the simple feature geometry list-column sfc together

# We create an sf object containing 2 single points, a set of points, and a polygon with one attribute, then maps by ggplot2

library(sf)

# Step 1
# Single point (as a vector)
p1_sfg <- st_point(c(2, 2))     # a point
p2_sfg <- st_point(c(2.5, 3))   # a point

# Set of points (as a matrix)
p <- rbind(c(6, 2), c(6.1, 2.6), c(6.8, 2.5),
           c(6.2, 1.5), c(6.8, 1.8))
mp_sfg <- st_multipoint(p) 

# Polygon - sequence of points that form a closed non-self intersecting ring
# The first ring denotes the exterior ring, zero or more subsequent rings denote holes in exterior ring
p1 <- rbind(c(10, 0), c(11, 0), c(13, 2),
            c(12, 4), c(11, 4), c(10, 0))           # exterior ring  
p2 <- rbind(c(11, 1), c(11, 2), c(12, 2), c(11, 1)) # interior ring  
pol_sfg <- st_polygon(list(p1, p2))

# Step 2
# Create sf object
p_sfc <- st_sfc(p1_sfg, 
                p2_sfg, 
                mp_sfg, 
                pol_sfg
)

# Step 3
df <- data.frame( v1 = c("A", "B", "C", "D") ) # the attributes  
p_sf <- st_sf(df,                              # the attributes
              geometry = p_sfc                 # the geometry 
)

# Plot them
library(ggplot2)
ggplot(p_sf) + 
  geom_sf( aes( col = v1 ), size = 3)


# Create a sf object containing two single points (1.5,1.5) and (2.5,1.5) with attribute values Peter and Bob,
# and a polygon ring where the exterior ring is a square side length 3, and interior ring a square whose side 
# has length 1 located in center of exterior ring. The attribute has value John, and centered at (1.5,1.5)
library(sf)

# STEP 1
# Single point (point as a vector)
p1_sfg <- st_point(c(1.5,1.5))     # a point
p2_sfg <- st_point(c(3.5,4.5))   # a point

# Polygon. Sequence of points that form a closed, non-self intersecting ring. The first ring denotes the exterior 
# ring, zero or more subsequent rings denote holes in the exterior ring
p1 <- rbind(c(0, 0), c(0, 3), c(3, 3), c(3,0), c(0,0) )   # exterior ring  
p2 <- rbind(c(1, 1), c(1, 2), c(2, 2), c(2, 1), c(1, 1)) # interior ring  
pol_sfg <- st_polygon(list(p1, p2))

# STEP 2
# Create sf object
p_sfc <- st_sfc(p1_sfg, 
                p2_sfg, 
                pol_sfg
)

# STEP 3
df <- data.frame( v1 = c("Peter", "Bob", "John") ) # the attributes  
p_sf <- st_sf(df, 
              geometry = p_sfc)

# Plot them (to be introduced later in the computer practical)
library(ggplot2)
ggplot(p_sf) + 
  geom_sf(aes( col = v1 ), size = 3) + 
  theme_bw()



## 8. Deleting and Combining Polygons

# Following script creates nc.map.1 from nc.map by removing polygons with FIPS values '37125', '37051'
library(sf)
path.to.nc.shp <- system.file('shape/nc.shp', package='sf')
nc.map <- st_read(path.to.nc.shp, quiet = TRUE)

# Delete polygon
nc.map.1 <- nc.map[-which(nc.map$FIPS %in% c('37125', '37051')), ]

# Plot
plot(nc.map.1['AREA'])


# Following creates nc.map.2 as a union of all nc.map.1 polygons

# Union
nc.map.2 <- st_union(nc.map.1, by_feature = FALSE)
nc.map.2 <- st_sf(nc.map.2)

# Plot
ggplot(st_union(nc.map.2, by_feature = FALSE)) +
  geom_sf()


# Create a sf object from nc.map where counties with NAME values Durham, Northampton, and Stokes
# are removed and plot
path.to.nc.shp <- system.file('shape/nc.shp', package='sf')
nc.map <- st_read(path.to.nc.shp, quiet=TRUE)

# Delete polygons
nc.map.1 <- nc.map[-which(nc.map$NAME %in% c('Durham', 'Northampton', 'Stokes')),]

# Plot
plot(nc.map.1['AREA'])



## 9. Counting the Number of Points Within Polygons

# 1. Uses function st_sample() to sample random points at locations in the map
# 2. Uses function st_intersects() to count the number of points within the polygons of a sf object.
# Returned object is a list with feature ids intersected in each of the polygons
# 3. Uses function lengths() to calculate the number of points inside each feature

# Map with divisions (sf object)
nc.map <- read_sf(system.file('shape/nc.shp', package='sf'))

# 1. Points over map (simple feature geometry list-column sfc)
nc.points <- st_sample(nc.map, size=100)

# Map of points within polygons
ggplot() +
  geom_sf(data=nc.map) +
  geom_sf(data=nc.points)

# 2. Intersection (first argument map, then points)
nc.inter <- st_intersects(nc.map, nc.points)

# 3. Add point count to each polygon
nc.map$count <- lengths(nc.inter)

# Map of number of points within polygons
library(ggplot2)
ggplot(nc.map) +
  geom_sf(aes(fill=count))


# Now count number of points in North Carolina whose area is smalelr than 0.1 only

# Map with divisions (sf object)
nc.map <- read_sf(system.file('shape/nc.shp', package='sf'))

# 1. Points over map (simple feature geometry list-column sfc)
nc.points <- st_sample(nc.map, size=100)

nc.map <- nc.map[-which(nc.map$AREA>0.1),]

# Map of points within polygons
ggplot() +
  geom_sf(data=nc.map)+
  geom_sf(data=nc.points)

# 2. Intersection (first argument map, then points)
nc.inter <- st_intersects(nc.map, nc.points)

# 3. Add point count to each polygon
nc.map$count <- lengths(nc.inter)

# Map of number of points within polygons
library(ggplot2)
ggplot(nc.map) +
  geom_sf(aes(fill = count))
