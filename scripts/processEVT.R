# processing and analyzing EVT data
# MAC 05/21/25

library(terra)
library(sf)
library(RColorBrewer)

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
# # load in Landfire EVT https://www.landfire.gov/vegetation/evt
# evt_raster <- rast("./data/landfire/EVT/LF2024_EVT_250_CONUS/LC24_EVT_250.tif")
# # reproject
# evt_raster <- project(evt_raster, shp)
# # crop to AOI
# evt_raster <- crop(evt_raster, aoi)
# # save to file
# writeRaster(evt_raster, filename = "./data/landfire/EVT/KNF_LF_EVT_2024.tif", overwrite=TRUE)

# set activeCat for plotting
evt_raster <- rast("./data/landfire/EVT/KNF_LF_EVT_2024.tif")
# categories
cats(evt_raster)
activeCat(evt_raster) <- "EVT_NAME" # set active layer to EVT_NAME, EVT_ORDER
plot(evt_raster)
plot(shp, add=TRUE, col=NA, border="red")  
# create single cat rast of EVT_NAME
evt_cats <- cats(evt_raster)[[1]]  # This must include 'VALUE' and 'EVT_NAME'
# Build a new category table using EVT_ORDER
evt_order_labels <- evt_cats[, c("VALUE", "EVT_NAME")]
colnames(evt_order_labels) <- c("ID", "label")  # terra expects 'ID' and 'label'
#Assign this label table to the raster
levels(evt_raster) <- evt_order_labels
levels(evt_raster)[[1]]$label


##### get EVT ORDER for mask
# Extract the full category table
evt_order <- rast("./data/landfire/EVT/KNF_LF_EVT_2024.tif")
evt_cats <- cats(evt_order)[[1]]  # This must include 'VALUE' and 'EVT_ORDER'

# Build a new category table using EVT_ORDER
evt_order_labels <- evt_cats[, c("VALUE", "EVT_ORDER")]
colnames(evt_order_labels) <- c("ID", "label")  # terra expects 'ID' and 'label'
#Assign this label table to the raster
levels(evt_order) <- evt_order_labels
levels(evt_order)[[1]]$label
# Create logical mask where TRUE = Tree-dominated
#tree_mask <- evt_order == "Tree-dominated"

# 1. Get the current level table
lvl_tbl <- levels(evt_order)[[1]]
# 2. Identify which level IDs correspond to "Tree-dominated"
tree_ids <- lvl_tbl$ID[lvl_tbl$label == "Tree-dominated"]
# 3. Create a logical mask based on underlying integer values
tree_mask <- values(evt_order)[,1] %in% tree_ids
# 4. Turn logical vector into a raster
tree_mask_raster <- setValues(rast(evt_order), tree_mask)
# 5. Plot to verify
plot(tree_mask_raster, main = "Tree-dominated Mask (by label ID match)")
# 6. Save the tree mask raster
writeRaster(tree_mask_raster, filename = "./data/landfire/EVT/KNF_LF_EVT_treeMask_2024.tif", overwrite=TRUE)

# Promote evt_raster to integer type that supports NA
evt_raster_safe <- classify(evt_raster, cbind(NA, NA), datatype = "INT2S")
# apply tree_mask to EVT NAME raster
non_tree_raster <- mask(evt_raster_safe, tree_mask_raster, maskvalues = TRUE)

#Extract the full label table from the original
evt_labels <- levels(evt_raster)[[1]]  # or evt_raster_safe
# Reattach it to the masked raster
levels(non_tree_raster) <- evt_labels
# Plot with labels
plot(non_tree_raster, type = "classes", main = "Non-Tree-Dominated EVT (Labeled)")

#### GET TOP 10 LC types
code_values <- values(non_tree_raster)

# Get original levels table (with EVT_NAME and code)
orig_levels <- levels(non_tree_raster)[[1]]

# Count frequency of codes
#code_freq <- sort(table(code_values), decreasing = TRUE)
code_freq <- freq(non_tree_raster)
top_codes <- code_freq$value[order(-code_freq$count)][1:5]

# Extract the levels table
levels_table <- levels(non_tree_raster)[[1]]

# Match the top10_labels to their IDs
top_ids <- levels_table$VALUE[levels_table$EVT_NAME %in% top_codes]

# Create a mask: TRUE for cells with ID in top10_ids
top_mask <- non_tree_raster %in% top_ids

# Mask the original raster to keep only top 10 classes
filtered_raster <- mask(non_tree_raster, top_mask, maskvalues = FALSE)

# Plot the filtered raster
plot(filtered_raster, type = "classes", main = "Top 10 Vegetation Types")

# save to file
writeRaster(filtered_raster, filename = "./data/landfire/EVT/KNF_LF_EVT_treeMask_top5_2024.tif", overwrite=TRUE)


