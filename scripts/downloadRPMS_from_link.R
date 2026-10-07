# download and crop RPMS data for the SW US using terra package
# adapted from original downloadRPMS.R on virtual machine
# MAC 10/31/2024
##### download from provided link -- adapted from /SWRangeClimate/scripts/downloadRPMS.R
# MAC 10/14/2025

library(terra)

# Set the download timeout to 600 seconds
options(timeout = 6000)

# Start the clock!
ptm <- proc.time()

# Download and process from link
  # Construct the URL and download the file
  print("downloading file")
  download.file(
    paste0("https://storage.googleapis.com/fuelcast-data/rpms-testing/ls-snc/merged/orig/rpms-max-ppa-ls-snc-2025-09-25.tif?X-Goog-Algorithm=GOOG4-RSA-SHA256&X-Goog-Credential=dask-raster-worker%40fuelcast.iam.gserviceaccount.com%2F20251014%2Fauto%2Fstorage%2Fgoog4_request&X-Goog-Date=20251014T182639Z&X-Goog-Expires=604800&X-Goog-SignedHeaders=host&X-Goog-Signature=11e855adf927ced63b9c1515b625f851856903bff0df32de55974c415557538e9819fc608858554cb815f0a05c6af90938b16b74fa7dcbe94143dbca517ddc6fcde08d08608ffbf5c137649b12bca364c40f8e8fa9b7abaad357e300df32b57f387fa4692d500fc09b20782a489fbfbdb029a22f4d2d03e077b16269a68e8ae7c4945cad44168f6d7fdbe7a4ea6d91547d445f7133cb2d0db380c6c09ad2880a53f0a91a5aa7f637c3ca20a161c15bec1de790e8679f9340ced1bc0df6fb497848528e471f4901b27f81104d7fafc1337f9ea5fe27f1d03c75eaba0f33b96a6b8193ddd33ca017b856f161cbc6f5f5a4ea7fd98a9d6a9d549604f90a2cb122c7"),
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
  yr=2025
  writeRaster(rpms, filename = paste0("./data/swRPMS_data/swRPMS_", yr, ".tif"), filetype = "GTiff", overwrite = TRUE)
  
  print(yr)
  unlink("temp.tif")


# Stop the clock and print the elapsed time
proc.time() - ptm
