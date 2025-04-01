
library(terra)

# Assume 'r' is your SpatRaster object.
# Set a threshold for veg_prod; adjust this value as needed.
veg_threshold <- 0  # example threshold

# Define a function that converts a pixel's values (across layers) to percentile ranks
percentile_fun <- function(x) {
  # Filter out NAs and values not exceeding the veg_prod threshold
  valid <- x[!is.na(x) & x > veg_threshold]
  
  # If no valid values, return NA for each element in x
  if (length(valid) == 0) {
    return(rep(NA, length(x)))
  }
  
  # Create an ECDF function using only the valid values
  ec <- ecdf(valid)
  
  # For each value in x:
  # - Return NA if the value is NA or not above the threshold.
  # - Otherwise, compute its percentile rank.
  sapply(x, function(val) {
    if (is.na(val) || val <= veg_threshold) {
      NA
    } else {
      ec(val) * 100  # Multiply by 100 to convert to percentile
    }
  })
}

# Apply the function to each pixel in the raster.
r_percentile <- app(r, fun = percentile_fun)

# 'r_percentile' now contains the percentile ranks for each pixel's values across layers.




