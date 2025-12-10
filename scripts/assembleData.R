# assemble data for RPMS ML analysis
# MAC 01/30/25

library(terra)
library(climateR)
library(AOI)
library(sf)

terraOptions(progress = 1)

#####
# set AOI
shp <- sf::st_read(dsn = "./data/shapes/AdministrativeForest.gdb")
shp<- sf::st_transform(shp, crs = sf::st_crs("+proj=longlat +datum=WGS84"))

# subset to selected NF
shp<-subset(shp, shp$FORESTNAME=="Kaibab National Forest")

# get ext
aoi<-ext(shp)+1
#####

#####
# process RPMS data

rpmsDir<-"./data/swRPMS_data"

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

# mask 2024 with NAs from previous years
# Apply mask (NA values in rpms_1984 will be applied to other_layer)
#rpmsStack[[40]] <- mask(rpmsStack[[40]], rpmsStack[[39]])

# set names on prod_stack layers
#names(rpmsStack) <- yrs[-29]

# apply mask of any NAs to all layers
# Step 1: Identify NA locations in each layer
  na_mask <- app(rpmsStack, fun = function(x) any(is.na(x)))
# Step 2: Convert logical mask to numeric (1 for valid, NA for missing)
  na_mask[na_mask == 1] <- NA
  na_mask[na_mask == 0] <- 1
# Step 3: Apply the NA mask to the original stack
test <- mask(rpmsStack, na_mask)

# set names on prod_stack layers
names(rpmsStack) <- yrs[-29]

# write data to file
writeRaster(rpmsStack, filename = "./data/processed/KNF_RPMS_1984_2024.tif", overwrite=TRUE)
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
writeRaster(monthly_totals, filename = "./data/processed/KNF_GridMet_monthly_pr_1984_2024.tif", overwrite=TRUE)

#####

#####
# process to common resolution

rpmsStack <- rast("./data/processed/KNF_RPMS_1984_2024.tif")
monthly_totals<-rast("./data/processed/KNF_GridMet_monthly_pr_1984_2024.tif")

# Set number of threads (adjust based on your system's cores)
terra::terraOptions()
terraOptions(memfrac = 0.8)

resmpRPMS <- resample(rpmsStack,monthly_totals, method="bilinear")

#names(resmpRPMS) <- names(rpmsStack)

writeRaster(resmpRPMS, filename = "./data/processed/4K_KNF_RPMS_1984_2024.tif", overwrite=TRUE)
#####

##### get Landfire landcover data
library(rlandfire)

# get AOI using landfire function
aoiLF <- getAOI(shp, extend = 1)
aoiLF

# set up query
# products? https://lfps.usgs.gov/lfps/helpdocs/productstable.html
products <- c("200BPS")
path <- tempfile(fileext = ".zip")
projection<-4326
#resolution<-60

resp <- landfireAPI(products = products,
                    aoi = aoiLF,
                    projection = projection, 
                    #resolution = resolution,
                    path = path,
                    verbose = TRUE)

lf_dir <- file.path(tempdir(), "lf")
utils::unzip(path, exdir = lf_dir)
# load as spatRast
lf <- terra::rast(list.files(lf_dir, pattern = ".tif$", 
                             full.names = TRUE, 
                             recursive = TRUE)[1])
# check by plotting
plot(lf)
activeCat(lf) <- 5 # set active layer to GROUPVEG
plot(lf)
plot(shp, add = TRUE)
#write out to save
writeRaster(lf, filename = "./data/processed/KNF_Landfire_200BPS.tif", overwrite=TRUE)

# load to resample
lf <- rast("./data/processed/KNF_Landfire_200BPS.tif")
activeCat(lf) <- 5 # set to GROUP VEG
monthly_totals<-rast("./data/processed/KNF_GridMet_monthly_pr_1984_2024.tif")

# Set number of threads (adjust based on your system's cores)
terra::terraOptions()
terraOptions(memfrac = 0.8)

# resample to gridmet
resmpLF <- resample(lf,monthly_totals, method="near")

writeRaster(resmpLF, filename = "./data/processed/4K_KNF_200BPS_GROUPVEG.tif", overwrite=TRUE)
#####

##### RAP cover data processing
# Path to directory containing rasters
  raster_dir <- "./data/climEng"
# List all raster files (adjust pattern if needed, e.g., ".tif")
  raster_files <- list.files(raster_dir, pattern = "\\.tif$", full.names = TRUE)
# Read all rasters and combine into a SpatRaster
  rasters <- rast(raster_files)
# Create a categorical map of the dominant land cover type
  dominant_cover <- which.max(rasters)
# Set factor levels using layer names
  categories <- names(rasters)
  dominant_cover <- as.factor(dominant_cover)
  levels(dominant_cover) <- data.frame(ID = 1:length(categories), class = categories)

# resample to gridmet
  resmpRAP <- resample(dominant_cover,monthly_totals, method="near")

  plot(resmpRAP, main="Dominant Land Cover Type -- 2024 RAP")
  plot(shp, add=TRUE, col=NA, border="red")  
      
  writeRaster(resmpRAP, filename = "./data/processed/4K_KNF_2024_dominantRAPcover.tif", overwrite=TRUE)

# plot the dominant land cover type    
plot(dominant_cover, main="Dominant Land Cover Type -- 2024 RAP")
plot(shp, add=TRUE, col=NA, border="red")  

##### process TEU data
# load TEUs
teu<-sf::st_read(dsn = "./data/shapes/TEU.gdb")
teu<-subset(teu, teu$PROJECT=="KAIBAB")
teu<- sf::st_transform(teu, crs = sf::st_crs("+proj=longlat +datum=WGS84"))

# test plot
library(ggplot2)
ggplot(data = teu) +
  geom_sf(aes(fill = PNV_CLASS),color = "grey90", linewidth = 0.1) +
  scale_fill_viridis_d(option = "D") +
  theme_minimal()
#

# load gridmet for resolution sample
monthly_totals<-rast("./data/processed/KNF_GridMet_monthly_pr_1984_2024.tif")

# rasterize teus
teu$PNV_CLASS <- as.factor(teu$PNV_CLASS)
teu_vect <- vect(teu)
teu_rast <- rasterize(teu_vect, monthly_totals[[1]], field = "PNV_CLASS")

plot(teu_rast, main="TEU PNV_CLASS")

writeRaster(teu_rast, filename = "./data/processed/4K_KNF_TEU_PNV_CLASS.tif", overwrite=TRUE)
