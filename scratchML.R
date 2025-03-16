# scratch code to develop ML model on RPMS with climate
# MAC 01/24/25

library(terra)
library(dplyr)
library(tidyr)
library(caret)
library(ggplot2)

# Parameters for synthetic data
nrows <- 100
ncols <- 100
xmin <- 0
ymax <- 100
res <- 1
years <- 1981:2023
months <- 1:12

# Create a raster template
ext <- ext(xmin, xmin + ncols * res, ymax - nrows * res, ymax)
template <- rast(ext, nrows = nrows, ncols = ncols, res = res, crs = "EPSG:4326")

# Generate Monthly Precipitation Layers
precip_layers <- list()
for (year in years) {
  for (month in months) {
    precip <- template
    values(precip) <- runif(nrows * ncols, min = 0, max = 200)
    layer_name <- sprintf("precip_%d_%02d", year, month)
    names(precip) <- layer_name
    precip_layers[[layer_name]] <- precip
  }
}
precip_stack <- rast(precip_layers)

# Generate Annual Production Layers
prod_layers <- list()
for (year in years) {
  prod <- template
  values(prod) <- runif(nrows * ncols, min = 100, max = 1000)
  layer_name <- sprintf("prod_%d", year)
  names(prod) <- layer_name
  prod_layers[[layer_name]] <- prod
}
prod_stack <- rast(prod_layers)

# Generate Vegetation Type Raster
veg_raster <- template
values(veg_raster) <- sample(1:3, nrows * ncols, replace = TRUE, prob = c(0.5, 0.3, 0.2))
names(veg_raster) <- "veg_type"

# Align Vegetation Raster with Precipitation Stack
veg_raster <- resample(veg_raster, precip_stack)

# Define the aggregation factor (e.g., 2x2 blocks)
agg_factor <- 4  # Each cell will represent a 2x2 block of original cells

# Reduce resolution of raster layers
precip_stack <- aggregate(precip_stack, fact = agg_factor, fun = mean, na.rm = TRUE)
prod_stack <- aggregate(prod_stack, fact = agg_factor, fun = mean, na.rm = TRUE)

# Custom modal function for categorical raster aggregation
modal_function <- function(x) {
  x <- na.omit(x)  # Remove NA values manually
  if (length(x) == 0) return(NA)  # If all values are NA, return NA
  ux <- unique(x)  # Get unique values
  ux[which.max(tabulate(match(x, ux)))]  # Return most frequent value (mode)
}

# Apply the aggregation with the custom modal function
veg_raster <- aggregate(veg_raster, fact = agg_factor, fun = modal_function)


# Combine Precipitation Data
precip_df <- as.data.frame(precip_stack, xy = TRUE, na.rm = TRUE) %>%
  pivot_longer(
    cols = starts_with("precip"),
    names_to = "layer",
    values_to = "precip"
  ) %>%
  separate(layer, into = c("prefix", "year", "month"), sep = "_") %>%
  mutate(year = as.integer(year), month = as.integer(month)) %>%
  select(-prefix)

# Combine Production Data
prod_df <- as.data.frame(prod_stack, xy = TRUE, na.rm = TRUE) %>%
  pivot_longer(
    cols = starts_with("prod"),
    names_to = "layer",
    values_to = "prod"
  ) %>%
  separate(layer, into = c("prefix", "year"), sep = "_") %>%
  mutate(year = as.integer(year)) %>%
  select(-prefix)

# Combine Vegetation Data
veg_df <- as.data.frame(veg_raster, xy = TRUE, na.rm = TRUE)
names(veg_df) <- c("x", "y", "veg_type")

# Merge All Data with Monthly Resolution
merged_df <- precip_df %>%
  left_join(prod_df, by = c("x", "y", "year")) %>%
  left_join(veg_df, by = c("x", "y")) %>%
  group_by(x, y, year) %>%
  mutate(
    precip_total = sum(precip, na.rm = TRUE),  # Total annual precipitation
    precip_seasonal = sum(precip[month %in% c(3:6)], na.rm = TRUE),  # Spring precipitation
    cum_precip = cumsum(precip),  # Cumulative precipitation
    lag_precip = lag(precip, n = 1)  # Lagged precipitation (previous month)
  ) %>%
  ungroup() %>%
  na.omit()

# Add Month as a Factor
merged_df <- merged_df %>%
  mutate(month = factor(month))

# Split Data for Training and Testing
# set.seed(123)
# train_indices <- createDataPartition(merged_df$prod, p = 0.8, list = FALSE)
# train_data <- merged_df[train_indices, ]
# test_data <- merged_df[-train_indices, ]
set.seed(123)
train_data<- merged_df %>% sample_frac(0.1)  # Use 10% of the training data
# Use anti_join to get the remaining 90% as testing data
test_data <- anti_join(merged_df, train_data, by = c("x", "y", "year", "month"))

# Train a Machine Learning Model
model_formula <- prod ~ precip + precip_total + precip_seasonal + cum_precip + lag_precip + veg_type + month
rf_model <- train(
  model_formula,
  data = train_data,
  method = "ranger",
  trControl = trainControl(method = "cv", number = 5),
  tuneLength = 3,
  num.trees = 10  # Reduce number of trees
)

# Evaluate Model
predictions <- predict(rf_model, newdata = test_data)
results <- data.frame(
  Observed = test_data$prod,
  Predicted = predictions
)
correlation <- cor(results$Observed, results$Predicted)
print(correlation)

# Visualization of Monthly Effects
ggplot(merged_df, aes(x = month, y = prod, color = veg_type)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "loess") +
  facet_wrap(~ veg_type) +
  labs(title = "Monthly Precipitation Effects on Annual Production", x = "Month", y = "Production")
