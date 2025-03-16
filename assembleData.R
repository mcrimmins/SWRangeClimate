# assemble data for RPMS ML analysis
# MAC 01/30/25

library(terra)
library(climateR)
library(AOI)
library(sf)

terraOptions(progress = 1)

#####
# set AOI
shp <- sf::st_read(dsn = "C:/Users/Crimmins/OneDrive - University of Arizona/RProjects/ClimateReports/Micro_Apps/Kaibab National Forest/Data/shapes/AdministrativeForest.gdb")
shp<- sf::st_transform(shp, crs = sf::st_crs("+proj=longlat +datum=WGS84"))

# subset to selected NF
shp<-subset(shp, shp$FORESTNAME=="Kaibab National Forest")

# get ext
aoi<-ext(shp)+1
#####

#####
# process RPMS data

rpmsDir<-"C:/Users/Crimmins/OneDrive - University of Arizona/RProjects/RPMSprocessing/swRPMS_data"

# Load files into SpatRaster collection
file_list <- list.files(path=rpmsDir, pattern=".tif", full.names=TRUE, recursive=FALSE)

# get range of years
yrs <- as.numeric(gsub(".*_(\\d{4})\\.tif$", "\\1", file_list))

# Drop 2012 from list, bad data?
file_list <- file_list[-29]
# Leave out most recent year
rpmsStack <- rast(file_list[1:(length(file_list))])
#rpmsStack <- rast(file_list[1:(length(file_list))])

# crop to aoi
rpmsStack <- crop(rpmsStack, aoi)

# set names on prod_stack layers
names(rpmsStack) <- yrs[-29]

# write data to file
writeRaster(rpmsStack, filename = "KNF_RPMS_1984_2024.tif", overwrite=TRUE)
#####

#####
# get gridmet data

# Define the latitude and longitude bounds
# -112.902832,34.997416,-111.556091,36.870834
# xmin <- -112.902832  # Minimum longitude
# xmax <- -111.556091   # Maximum longitude
# ymin <- 34.997416    # Minimum latitude
# ymax <- 36.870834    # Maximum latitude

# Create a bounding box as an sf object
#bbox <- st_bbox(c(xmin = xmin, ymin = ymin, xmax = xmax, ymax = ymax), 
#                crs = st_crs(4326))  # EPSG:4326 is for WGS 84
bbox <- st_bbox(aoi, crs = st_crs(4326))  # EPSG:4326 is for WGS 84

# Convert to an sf geometry for visualization or further manipulation
bbox_sf <- st_as_sfc(bbox)

# Print the bounding box
print(bbox_sf)

aoi<-bbox_get(bbox_sf)
# Visualize (optional)

system.time({
  gridmet_pr = getGridMET(AOI = aoi,
                          varname = "pr",
                          startDate = "1984-01-01",
                          endDate  = "2024-12-31")
})

layer_months <- format(as.Date(names(gridmet_pr$precipitation_amount),"pr_%Y-%m-%d"), "%Y-%m")

monthly_totals <- terra::tapp(gridmet_pr$precipitation_amount, index = layer_months, fun = sum)

# write data to file
writeRaster(monthly_totals, filename = "./data/KNF_GridMet_monthly_pr_1984_2024.tif", overwrite=TRUE)

#####

#####
# process to common resolution

rpmsStack <- rast("./data/KNF_RPMS_1984_2024.tif")
monthly_totals<-rast("./data/KNF_GridMet_monthly_pr_1984_2024.tif")

# Set number of threads (adjust based on your system's cores)
terra::terraOptions()
terraOptions(memfrac = 0.8)

resmpRPMS <- resample(rpmsStack,monthly_totals, method="bilinear")

#names(resmpRPMS) <- names(rpmsStack)

writeRaster(resmpRPMS, filename = "./data/4K_KNF_RPMS_1984_2024.tif", overwrite=TRUE)
#####



