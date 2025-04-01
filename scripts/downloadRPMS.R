# download and crop RPMS data for the SW US using terra package
# adapted from original downloadRPMS.R on virtual machine
# MAC 10/31/2024

library(terra)

# Set the download timeout to 600 seconds
options(timeout = 600)

# Start the clock!
ptm <- proc.time()

# Loop through the years
for (yr in 2024:2024) {
  # Construct the URL and download the file
  print("downloading file")
  download.file(
    paste0("https://fuelcast.net/rpms-download?date=", yr, "&file=rpms_", yr, ".tif"),
    destfile = "temp.tif", extra = "--no-verbose", mode = "wb"
  )
  
  # Load the raster into a SpatRaster object
  rpms <- rast("temp.tif")
  
  # Crop the raster to the specified extent
  print("cropping")
  extent_crop <- ext(-114.982910, -102.864990, 31.269161, 37.072710)
  rpms <- crop(rpms, extent_crop)
  
  # Save the cropped raster as a GeoTIFF
  print("saving cropped tif")
  writeRaster(rpms, filename = paste0("./data/swRPMS_data/swRPMS_", yr, ".tif"), filetype = "GTiff", overwrite = TRUE)
  
  print(yr)
  unlink("temp.tif")
}

# Stop the clock and print the elapsed time
proc.time() - ptm

