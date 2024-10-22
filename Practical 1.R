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
