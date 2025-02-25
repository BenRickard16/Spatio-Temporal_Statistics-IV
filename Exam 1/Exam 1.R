library(sf)
library(spData)
library(ggplot2)
library(mapview)
library(tmap)
library(gstat)
library(spdep)
library(spatialreg)

###############

data.df <- data.frame(ID = c(1,2,3,4,5), 
                 Longitude = c(27389845, 26421228, 27050077, 28437658, 26681283),
                 Latitude = c(17287911, 17513182, 16815791, 17251164, 17217992),
                 Value = c(-1,-2,-4,-6,-9))

data.sf <- st_as_sf(data.df, coords = c('Longitude', 'Latitude'), crs = 6576)

data.map <- data.sf
tmap_mode('view')

map <- tm_shape(data.map) + 
  tm_dots(col = "Value", 
          palette = "Set1", 
          size = 0.5, 
          title = "Value") + 
  tm_basemap('CyclOSM')


##########
p1_sfg <- st_point(c(-80.1918, 25.7617))     # a point
p2_sfg <- st_point(c(-66.1057, 18.4655))
p3_sfg <- st_point(c(-64.7507, 32.3078))
p1 <- rbind(c(-80.1918, 25.7617), c(-66.1057, 18.4655), c(-64.7507, 32.3078), c(-80.1918, 25.7617) )  
pol_sfg <- st_polygon(list(p1))
# Create sf object
p_sfc <- st_sfc(p1_sfg, 
                p2_sfg,
                p3_sfg,
                pol_sfg)

df <- data.frame( v1 = c('Miami', 'San Juan', 'Bermuda', 'Bermuda Triangle') ) 
p_sf <- st_sf(df, 
              geometry = p_sfc)


library(ggplot2)
ggplot(p_sf) + 
  geom_sf(aes( col = v1 ), size = 3) + 
  theme_bw()

#################


cyprus.districts <- sf::st_read("./fcbt-2024-sdo-cyprus.districts.shp", quiet=TRUE) 

cyprus.sites <- sf::st_read("./fcbt-2024-sdo-cyprus.sites.shp", quiet=TRUE)

# Intersection (first argument map, then points)
gk.inter <- st_intersects(cyprus.districts, cyprus.sites )
# Add point count to each polygon
cyprus.districts$count <- lengths(gk.inter)
print(cyprus.districts$count)



######################
download.file(url="https://www.maths.dur.ac.uk/users/georgios.karagiannis/temp1/fcbt-st4-2024/fcbt-2024-prd.csv", "./fcbt-2024-prd.csv") 

usa.pfd.df <- read.csv(file='./fcbt-2024-prd.csv',sep=',')

usa.pfd.sf <- st_as_sf(usa.pfd.df,
                       
                       coords = c("longitude", "latitude"),
                       
                       crs = 4326)

frml <- Z~1
usa.cloud <- variogram(object = frml, data = usa.pfd.sf, cloud=TRUE)
plot(usa.cloud,
     xlab = 'distance', 
     ylab = 'gamma', 
     main = 'semivariogram cloud')

####

frml <- Z ~ X
smpl.vgm <- variogram(object = frml,
                      data = usa.pfd.sf)

vgm.obj <- vgm(psill = 0.48,
               model = "Nug", 
               range = 0)

vgm.obj <- vgm(psill = 1.55,
               model = "Exp",
               range = 243,
               add.to = vgm.obj)

# 2 Fit the semi variogram model against the sample variogram
est.vgm <- fit.variogram(object = smpl.vgm,
                         model = vgm.obj)

# 3 Print the estimated values for the parameters of the semi variogram model   
print(est.vgm)

# 4 plot it 0
plot(est.vgm, 
     cutoff = 1500)

####

krige.fit <- gstat( formula = frml,
                    data = usa.pfd.sf, 
                    model = est.vgm)

new.df <- data.frame(X = 0,
                    longitude = c(-92.1735), 
                    latitude = c(38.5767) )

new.sf <- st_as_sf(new.df, 
                               coords = c("longitude", "latitude"),
                               crs = 4326)

krige.prd <- predict(object = krige.fit,
                     newdata = new.sf)

krige.prd


######################

usa.data <- sf::st_read("./fcbt-2024-aud.shp")

neighbours <- poly2nb(pl=usa.data, queen=TRUE)
prox.matrix <- nb2listw(neighbours, style='W')

global.moran <- moran.test(usa.data$Z, prox.matrix, alternative = 'greater')
####################

local.moran <- localmoran(usa.data$Z, prox.matrix, alternative = 'two.side')


library(tmap)  # load the library  

tmap_mode(
  mode = "plot"   # this is for printing the standard mode plot (not interactive)  
  #  mode = "view"  # this is for printing the interactive plot  
)               # set the mode to be printable (not the interactive)



# variable of interest ####################  


# local Moran's I   ##########  

## append the local Moran's I in the sf data.frame
usa.data$lmI <- local.moran[, "Ii"] 
## tm_shape creates a tmap-element that specifies a spatial data object, which we refer to as shape.
p2 <- tm_shape(shp = usa.data    # shape object
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
usa.data$lmZ <- local.moran[, "Z.Ii"] 
## tm_shape creates a tmap-element that specifies a spatial data object, which we refer to as shape.
p3 <- tm_shape(shp = usa.data    # shape object
)                 
## tm_polygons fills the polygons and draws the polygon borders.
p3 <- p3 + tm_polygons(col = "lmZ",       # the name of a data variable 
                       title = "Z-score",          # title 
                       breaks = c(-Inf, 1.96, Inf) # method to process the color scale when col is a numeric variable, here by setting the breaks manually 
)                  
## tm_layout specifies the map layout; e.g. controls title, margins, aspect ratio, colors, frame, legend, among many other things.
p3 <- p3 +  tm_layout(legend.outside = TRUE  # determines whether the legend is plot outside of the map/facets. 
)          
## Arrange small multiples in a grid layout. 
tmap_arrange( p2, p3)



########################

usa.data <- sf::st_read("./fcbt-2024-aud.shp")

neighbours <- poly2nb(pl=usa.data, queen=TRUE)
prox.matrix <- nb2listw(neighbours, style='B')

frml <- Y ~ X
car.model <- spautolm(formula = frml, data = usa.data, listw = prox.matrix, family='CAR')
summary(car.model)


usa.data$signal_stochastic <- car.model$fit$signal_stochastic

library(tmap)
tmap_mode("plot")
tm_shape(usa.data) + 
  tm_polygons("signal_stochastic") 
