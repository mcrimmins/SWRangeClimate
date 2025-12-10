# get RAP cover from Climate Engine
# code adapted from climEng_RAP.R
# MAC 02/16/2025
# links to info
# https://support.climateengine.org
# https://docs.climateengine.org/docs/build/html/overview.html
# https://api.climateengine.org/docs


library(raster)
library(httr) # HTTP API requests
library(httr2) # HTTP API requests
library(googleCloudStorageR)
library(tidyverse)

httr::set_config(config(timeout(120))) 

# load key, proj info
source('./scripts/climEngKey.R')

# google bucket info
bucket<-"clim-engine"
# data directory
dataDir<-"climEng"


#####
# set AOI bbox
shp <- sf::st_read(dsn = "./data/shapes/AdministrativeForest.gdb")
shp<- sf::st_transform(shp, crs = sf::st_crs("+proj=longlat +datum=WGS84"))

# subset to selected NF
shp<-subset(shp, shp$FORESTNAME=="Kaibab National Forest")

# get ext
aoi<-terra::ext(shp)+1
#####
bbox<-paste0("[",aoi[1],",",aoi[3],",",aoi[2],",",aoi[4],"]")
#bbox<-"[-111.809578,35.141283,-111.798720,35.149530]"
#centroid<-paste0("[[",(-111.809578+-111.798720)/2,",",(35.141283+35.149530)/2,"]]")

##### get time series -----
# get timeseries of RAP data
# following https://github.com/Google-Drought/SupportSiteTutorials/blob/a038ee1e69fff32008d8619c0acc6db082d5795d/Timeseries/Blends_Example.Rmd

# Define root url for Climate Engine API
# root_url <- 'https://api.climateengine.org/'
# # Define endpoint for initial data request
# endpoint <- "timeseries/native/coordinates"
# 
# # Define API arguments time-series endpoint to get long-term blend data 
# query <- list(dataset = 'RAP_PRODUCTION_16DAY',
#               variable = "herbaceousAGB",
#               start_date = '1986-01-01',
#               end_date = Sys.Date(),
#               #end_date = '1990-01-01',
#               #buffer = '1000',
#               #coordinates = paste0("[[",centroids@coords[1,1],",",centroids@coords[1,2],"]]"),
#               coordinates = centroid,
#               area_reducer = 'mean')
# 
# # Run GET request to get data
# getTS <- GET(paste0(root_url, endpoint), config = add_headers(Authorization = key), query = query)
# print(getTS)
#####

##### alt download ----
# library(httr2)
# req <- request(paste0(root_url, endpoint)) |>
#   req_headers(Authorization = key) |>
#   req_url_query(!!!query) |>
#   req_timeout(seconds = 120) |>
#   req_perform()
# print(resp_body_json(req))
#####

##### download raw RAP raster ----
# Define root url for Climate Engine API
root_url <- 'https://api.climateengine.org'

tempFile<-"rapLTRcover"
exportPath<-paste0(bucket,"/",tempFile)

print(paste0("Processing ", tempFile))

endpoint = '/raster/export/values'


# RAP cover variables https://docs.climateengine.org/docs/build/html/variables.html#rst-rap-cover-30m-yearly

##### RAP cover -----
query <- list(dataset = 'RAP_COVER',
              variable = "LTR", # AFG, PFG, SHR, TRE, BGR, LTR 
              temporal_statistic = "mean",
              bounding_box = bbox,
              export_path = exportPath,
              export_resolution = 4000,
              start_date = '2024-01-01',
              end_date = '2024-12-31'
)

# Run GET request to get data
get_raster <- GET(paste0(root_url, endpoint), config = add_headers(Authorization = key), query = query)
print(get_raster)

# alt download
req <- request(paste0(root_url, endpoint)) |>
  req_headers(Authorization = key) |>
  req_url_query(!!!query) |>
  req_timeout(seconds = 120) |>
  req_perform()
print(resp_body_json(req))


# download file from google cloud
objG<-gcs_list_objects(bucket)

# wait for processing
ptm <- proc.time()
while(length(objG)==0){
  print("waiting for raster to process")
  objG<-gcs_list_objects(bucket)
  Sys.sleep(10)
}
rawTime<-proc.time() - ptm

# download from bucket
gcs_get_object(objG$name[[which(objG$name==paste0(tempFile,".tif"))]], saveToDisk = paste0("./data/",dataDir,"/",tempFile,".tif"), bucket = bucket, overwrite = TRUE)
#####

# delete file
gcs_delete_object(objG$name[[which(objG$name==paste0(tempFile,".tif"))]])

#json_string <- jsonlite::toJSON(query, pretty = TRUE)  # Use pretty=TRUE for better readability
#print(json_string)

# load tif as a test
r<-terra::rast(paste0("./data/",dataDir,"/",tempFile,".tif"))
terra::plot(r)