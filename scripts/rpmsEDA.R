# RPMS and climate data EDA
# MAC 03/19/2025

library(terra)
library(sf)
library(rasterVis)
library(ggplot2)

# load grids
rpms<-rast("./data/processed/KNF_RPMS_1984_2024.tif")

# Load SpatRaster data
veg_prod <- rast("./data/processed/4K_KNF_RPMS_1984_2024.tif")  # Annual vegetation production
precip <- rast("./data/processed/KNF_GridMet_monthly_pr_1984_2024.tif")  # Monthly precipitation
veg_type <- rast("./data/processed/4K_KNF_200BPS_GROUPVEG.tif")  # Vegetation type
activeCat(veg_type)<-5

# set new extent to all spatRast
newExt<-ext(veg_prod)
newExt<-ext(newExt[1],newExt[2],newExt[3],newExt[4]-1)
veg_prod<-crop(veg_prod,newExt)
precip<-crop(precip,newExt)
veg_type<-crop(veg_type,newExt)

# load boundary
shp <- sf::st_read(dsn = "./data/shapes/AdministrativeForest.gdb")
shp<- sf::st_transform(shp, crs = sf::st_crs("+proj=longlat +datum=WGS84"))

# subset to selected NF
shp<-subset(shp, shp$FORESTNAME=="Kaibab National Forest")

#### make some maps
# plot rpms
plot(rpms[[39]], main="RPMS 2023")
# add shp with no fill
plot(shp, add=TRUE, col=NA, border="red")

# plot precip
plot(precip[[480]], main="Precipitation 12-2023")
# add shp with no fill
plot(shp, add=TRUE, col=NA, border="red")

# plot veg_prod
plot(veg_prod[[39]], main="RPMS 4km 2023")
# add shp with no fill
plot(shp, add=TRUE, col=NA, border="red")

# plot all veg prod
levelplot(veg_prod, main="RPMS 4km 1984-2024")

# plot time series hist
bwplot(veg_prod, main="RPMS 4km 1984-2024")

# percentile rank of veg_prod
veg_percentile <- app(veg_prod, fun = percentile_fun)
names(veg_percentile)<-names(veg_prod)

# hist plot of percentile rank
bwplot(veg_percentile, main="RPMS 4km 1984-2024 percentile rank")

# plot of yearly medians
  # Compute the median for each layer while removing NAs.
  medians <- global(veg_prod, fun = median, na.rm = TRUE)
  # Extract the layer names (assumed to be years) and convert them to numeric.
  years <- as.numeric(names(veg_prod))
  # Create a data frame with the years and corresponding medians.
  df <- data.frame(year = years, median = medians$global)
  # Plot the time series using ggplot2.
  ggplot(df, aes(x = year, y = median)) +
    geom_line() +
    geom_point() +
    theme_minimal() +
    labs(x = "Year", y = "Median RPMS", title = "Median Annual RPMS 4km 1984-2024")

# plot of vegetation types
plot(veg_type, main="Landfire BPS VegGroups -- 4km")
# add shp with no fill
plot(shp, add=TRUE, col=NA, border="black")
# with mask
plot(mask(veg_type, veg_prod[[40]]))






