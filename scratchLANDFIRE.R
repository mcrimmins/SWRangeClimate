# trying out LANDFIRE data API
# MAC 03/15/2025
# using https://github.com/bcknr/rlandfire

library(rlandfire)
library(sf)
library(terra)

#####
# set AOI
shp <- sf::st_read(dsn = "C:/Users/Crimmins/OneDrive - University of Arizona/RProjects/ClimateReports/Micro_Apps/Kaibab National Forest/Data/shapes/AdministrativeForest.gdb")
shp<- sf::st_transform(shp, crs = sf::st_crs("+proj=longlat +datum=WGS84"))

# subset to selected NF
shp<-subset(shp, shp$FORESTNAME=="Kaibab National Forest")

# get ext
#aoi<-terra::ext(shp)+1

# get AOI using landfire function
aoi <- getAOI(shp, extend = 0.5)
aoi

# set up query
# products? https://lfps.usgs.gov/lfps/helpdocs/productstable.html
products <- c("200BPS")
path <- tempfile(fileext = ".zip")
projection<-4326
#resolution<-60

resp <- landfireAPI(products = products,
                    aoi = aoi,
                    projection = projection, 
                    #resolution = resolution,
                    path = path,
                    verbose = TRUE)

lf_dir <- file.path(tempdir(), "lf")
utils::unzip(path, exdir = lf_dir)

lf <- terra::rast(list.files(lf_dir, pattern = ".tif$", 
                             full.names = TRUE, 
                             recursive = TRUE)[1])
plot(lf)

########
library(terra)
library(ggplot2)
library(RColorBrewer)

# Load raster (modify the path accordingly)
r<-lf

activeCat(r) <- 5

# Convert raster to dataframe for ggplot
r_df <- as.data.frame(r, xy = TRUE)

