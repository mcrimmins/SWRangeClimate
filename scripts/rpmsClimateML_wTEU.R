# experimental ML code on monthly precipitation and vegetation production data
# USING TEUs, adapted from rpmsClimateML_wBPS.T
# MAC 5/8/25

library(terra)   # For spatial data
library(dplyr)   # For data manipulation
library(tidyr)   # For reshaping data
#library(ranger)  # Random forest modeling
library(ggplot2) # Visualization
# library(viridis) # For color scales
# library(xgboost) # For XGBoost modeling
# library(mgcv)    # For Generalized Additive Models (GAMs)
# library(GWmodel) # For Geographically Weighted Regression (GWR)
# library(stringr) # For string manipulation
#library(caret)
#library(SPEI)
library(mgcv)

terraOptions(progress=1)

# Load SpatRaster data
#veg_prod <- rast("./data/processed/4K_KNF_RPMS_1984_2024.tif")  # Annual vegetation production
#precip <- rast("./data/processed/KNF_GridMet_monthly_pr_1984_2024.tif")  # Monthly precipitation
#veg_type <- rast("./data/processed/4K_KNF_TEU_PNV_CLASS.tif")  # Vegetation type

# Load SpatRaster data
#veg_prod <- rast("./data/processed/800m_KNF_RPMS_1984_2024.tif")  # Annual vegetation production
veg_prod <- rast("./data/processed/800m_KNF_RPMS_1984_2024_EVTmask.tif")  # Annual vegetation production with EVT mask
precip <- rast("./data/processed/KNF_PRISM_800m_monthly_pr_1984_2024.tif")  # Monthly precipitation
#veg_type <- rast("./data/processed/800m_KNF_TEU_PNV_CLASS.tif")  # Vegetation type
veg_type <- rast("./data/processed/800m_KNF_LF_EVT_treeMask_top5_2024.tif")  # Vegetation type


# KNF boundary for plotting
shp <- sf::st_read(dsn = "./data/shapes/AdministrativeForest.gdb")
shp<- sf::st_transform(shp, crs = sf::st_crs("+proj=longlat +datum=WGS84"))
# subset to selected NF
shp<-subset(shp, shp$FORESTNAME=="Kaibab National Forest")

# mask out veg_type by veg_prod
veg_type<-mask(veg_type, veg_prod[[40]])

# make map of veg type
veg_df <- as.data.frame(veg_type, xy = TRUE, na.rm = TRUE)
ggplot(veg_df, aes(x = x, y = y, fill = droplevels(as.factor(EVT_NAME)))) +
  geom_sf(data = shp, fill = NA, color = "black", inherit.aes = FALSE) +
  geom_raster() +
  coord_sf() +
  theme_minimal() +
  labs(title = "Vegetation Type Map", x = "Longitude", y = "Latitude")


# apply veg_type as mask to precip and to veg_prod
precip <- mask(precip, veg_type)
veg_prod <- mask(veg_prod, veg_type)

##### prepare data
# Get year labels from names
months <- names(precip)
years <- unique(substr(months, 2, 5))

# Split into yearly stacks
precip_yearly <- lapply(years, function(yr) {
  lyr_idx <- grep(paste0("", yr), names(precip)) # paste X in if needed
  precip[[lyr_idx]]
})
names(precip_yearly) <- years

# Compute cumulative monthly precipitation for each year
#precip_cumsum_yearly <- lapply(precip_yearly, function(p) app(p, fun = cumsum))

# Also compute total annual precipitation
precip_total_annual <- sds(lapply(precip_yearly, function(p) sum(p)))

# Stack annual production and vegetation into a single list
data_stack <- lapply(1:nlyr(veg_prod), function(i) {
  c(veg_prod[[i]], precip_total_annual[[i]], veg_type)
})

# Get year labels from raster names (remove year 2012 from precip set to match veg_prod)
#years_precip <- unique(substr(names(precip), 2, 5))
years_precip <- unique(substr(names(precip), 1, 4))
years_prod <- names(veg_prod)
valid_years <- years_prod  # assume 2012 is already excluded from veg_prod

# Align annual precip stack to veg_prod years (excludes 2012)
precip_total_matched <- precip_total_annual[match(valid_years, years_precip)]

# Extract veg_type category labels
veg_lookup <- levels(veg_type)[[1]]  # should contain ID and label (e.g., PNV_CLASS)
colnames(veg_lookup) <- c("value", "label")  # rename for clarity

# function to get month indices for early and late growing seasons
get_month_indices <- function(year, start_yr = 1984) {
  offset <- (as.numeric(year) - start_yr) * 12
  list(
    early = offset + 1:6,   # Jan to Jun
    late  = offset + 7:10   # Jul to Oct
  )
}

# Construct all_data with mapping and cleaning in one go
all_data <- lapply(seq_along(valid_years), function(i) {
  # Indexes for monthly layers
  idxs <- get_month_indices(valid_years[i])
  
  # Extract raster values (as vectors)
  prod_vals <- as.vector(values(veg_prod[[i]]))
  prec_vals <- as.vector(values(precip_total_matched[[i]]))
  veg_vals  <- as.vector(values(veg_type))
  
  early_vals <- as.vector(values(sum(precip[[idxs$early]])))
  late_vals  <- as.vector(values(sum(precip[[idxs$late]])))
  
  tibble(
    year = as.numeric(valid_years[i]),
    veg_prod = prod_vals,
    annual_precip = prec_vals,
    veg_type_code = as.integer(veg_vals),
    early_precip = early_vals,
    late_precip = late_vals
  )
}) %>%
  bind_rows() %>%
  left_join(veg_lookup, by = c("veg_type_code" = "value")) %>%
  #rename(veg_type = PNV_CLASS) %>%
  #rename(veg_type = EVT_NAME) %>%
  mutate(
    #veg_type = sub("\\s.*", "", veg_type),
    #veg_type = factor(veg_type),
    early_frac = early_precip / annual_precip,
    late_frac  = late_precip / annual_precip
  ) %>%
  drop_na(veg_prod, annual_precip, veg_type_code, early_precip, late_precip)
# Convert veg_type_code to factor with labels
all_data$veg_type <- factor(all_data$veg_type_code, levels = veg_lookup$value, labels = veg_lookup$label)
#####


##### EDA ----
ggplot(all_data, aes(x = veg_prod)) +
  geom_histogram(bins = 50, fill = "forestgreen", color = "grey") +
  #scale_x_log10() +  # optional if right-skewed
  facet_wrap(~veg_type, scales="free_y") +
  labs(title = "Distribution of Forage Production", x = "lbs/ac", y = "Count")+
  theme_bw()

all_data_long <- all_data %>%
  pivot_longer(cols = c(early_precip, late_precip), names_to = "season", values_to = "precip")

ggplot(all_data_long, aes(x = precip, fill = season)) +
  geom_histogram(bins = 50, position = "identity", alpha = 0.6) +
  facet_wrap(~veg_type) +
  labs(title = "Seasonal Precipitation by Vegetation Type", x = "Precip (mm)", y = "Count")

ggplot(all_data, aes(x = early_precip, y = veg_prod)) +
  #geom_point(alpha = 0.3) +
  geom_hex(bins = 50) +  # adjust bins for resolution
  #geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~veg_type, scales = "free") +
  labs(title = "Forage Production vs. Jan-Jun Precip", x = "Jan-Jun Precip (mm)", y = "Forage (lbs/ac)")

#####
# library(ggpointdensity)
# ggplot(all_data, aes(x = early_precip, y = veg_prod)) +
#   geom_pointdensity(adjust = 0.5) +
#   facet_wrap(~veg_type, scales = "free") +
#   scale_color_viridis_c() +
#   labs(title = "Forage Production vs. Jan-Jun Precip",
#        x = "Jan-Jun Precip (mm)", y = "Forage (lbs/ac)")
#####

ggplot(all_data, aes(x = late_precip, y = veg_prod)) +
  #geom_point(alpha = 0.3) +
  geom_hex(bins = 50) +  # adjust bins for resolution
  #geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~veg_type, scales = "free") +
  labs(title = "Forage Production vs. Jul-Oct Precip", x = "Jul-Oct Precip (mm)", y = "Forage (lbs/ac)")

ggplot(all_data, aes(x = annual_precip, y = veg_prod)) +
  #geom_point(alpha = 0.3) +
  geom_hex(bins = 50) +  # adjust bins for resolution
  #geom_smooth(method = "loess", se = FALSE) +
  facet_wrap(~veg_type, scales = "free") +
  labs(title = "Forage Production vs. Annual Precip", x = "Annual Precip (mm)", y = "Forage (lbs/ac)")

ggplot(all_data, aes(x = early_precip, y = late_precip)) +
  geom_density_2d_filled(alpha = 0.7) +
  facet_wrap(~veg_type) +
  labs(title = "Density of Seasonal Precip Combinations", x = "Early Precip", y = "Late Precip")

corrs<-all_data %>%
  group_by(veg_type) %>%
  summarise(
    cor_early = cor(veg_prod, early_precip, use = "complete.obs", method="spearman"),
    cor_late = cor(veg_prod, late_precip, use = "complete.obs", method="spearman"),
    mean_prod = mean(veg_prod),
    n = n()
  )
#####


##### MODELING
# modeling with ANNUAL PRECIP ----
mod_lm <- lm(veg_prod ~ annual_precip * veg_type, data = all_data)
summary(mod_lm)

ggplot(all_data, aes(x = annual_precip, y = veg_prod)) +
  #geom_point(alpha = 0.3) +
  geom_hex(bins = 50) +  # adjust bins for resolution
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~veg_type) +
  theme_minimal() +
  labs(x = "Annual Precipitation (mm)", y = "Forage Production (kg/ha)")

## GAM with Annual Precipitation
mod_gam <- gam(veg_prod ~ s(annual_precip, by = veg_type) + veg_type,
               data = all_data, method = "REML")
summary(mod_gam)
plot(mod_gam, pages = 1, rug = TRUE)

AIC(mod_lm, mod_gam)

# Create new data for prediction
new_data <- expand.grid(
  annual_precip = seq(min(all_data$annual_precip), max(all_data$annual_precip), length.out = 100),
  veg_type = levels(all_data$veg_type)
)

# Predict from GAM
new_data$predicted_prod <- predict(mod_gam, newdata = new_data)

# Plot
ggplot(new_data, aes(x = annual_precip, y = predicted_prod, color = veg_type)) +
  geom_line(size = 1) +
  theme_minimal() +
  labs(x = "Annual Precipitation (mm)", y = "Predicted Forage Production (kg/ha)")

# Calculate thresholds for 80% of max production
thresholds <- new_data %>%
  group_by(veg_type) %>%
  mutate(pct_max = predicted_prod / max(predicted_prod)) %>%
  filter(pct_max >= 0.8) %>%
  summarise(min_precip_for_80pct = min(annual_precip))

# Function to calculate thresholds for % of max with bootstrapping
library(mgcv)
library(purrr)
library(tibble)
library(pbapply)
bootstrap_threshold <- function(data, veg, n = 500, cutoff = 0.8) {
  thresholds <- pbapply::pbreplicate(n, {
    samp <- data %>% 
      filter(veg_type == veg) %>% 
      sample_frac(replace = TRUE)
    
    mod <- gam(veg_prod ~ s(annual_precip), data = samp)
    
    new_precip <- seq(min(samp$annual_precip), max(samp$annual_precip), length.out = 200)
    preds <- predict(mod, newdata = data.frame(annual_precip = new_precip))
    
    max_prod <- max(preds, na.rm = TRUE)
    which_val <- which(preds >= max_prod * cutoff)[1]
    
    if (!is.na(which_val)) {
      new_precip[which_val]
    } else {
      NA_real_
    }
  })
  
  tibble(
    veg_type = veg,
    mean_thresh = mean(thresholds, na.rm = TRUE),
    ci_lower = quantile(thresholds, 0.025, na.rm = TRUE),
    ci_upper = quantile(thresholds, 0.975, na.rm = TRUE)
  )
}

# Apply the bootstrap threshold function to each vegetation type
veg_types <- unique(all_data$veg_type)

threshold_results <- map_dfr(veg_types, ~bootstrap_threshold(all_data, .x, n = 100))

print(threshold_results)

# Optional: plot thresholds with uncertainty
ggplot(threshold_results, aes(x = veg_type, y = mean_thresh)) +
  geom_point() +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.2) +
  theme_minimal() +
  labs(x = "Vegetation Type", y = "Precip. Threshold (80% max prod)", 
       title = "Estimated Precipitation Thresholds (80%) with 95% CI")

##### EARLY VS LATE precip analysis
library(mgcv)
library(pbapply)
library(dplyr)

bootstrap_joint_thresholds <- function(data, veg, n = 200, cutoff = 0.75) {
  # Fixed prediction grid
  early_seq <- seq(min(data$early_precip, na.rm = TRUE),
                   max(data$early_precip, na.rm = TRUE), length.out = 50)
  late_seq <- seq(min(data$late_precip, na.rm = TRUE),
                  max(data$late_precip, na.rm = TRUE), length.out = 50)
  base_grid <- expand.grid(early_precip = early_seq,
                           late_precip = late_seq)
  
  # Resampling function
  run_one <- function(i) {
    samp <- data %>%
      filter(veg_type == veg) %>%
      sample_frac(replace = TRUE)
    
    mod <- gam(veg_prod ~ te(early_precip, late_precip),
               data = samp, family = Gamma(link = "log"), method = "REML")
    
    if (inherits(mod, "try-error")) return(NULL)
    
    grid <- base_grid
    grid$pred <- predict(mod, newdata = grid, type = "response")
    maxval <- max(grid$pred, na.rm = TRUE)
    
    grid %>%
      mutate(
        pct_max = pred / maxval,
        boot_id = i
      ) %>%
      filter(pct_max >= cutoff)
  }
  
  # Apply with progress bar
  boot_results <- pbapply::pblapply(1:n, run_one)
  
  bind_rows(boot_results) %>% mutate(veg_type = veg)
}

veg_types <- unique(all_data$veg_type)

boot_thresholds_all <- map_dfr(veg_types, ~bootstrap_joint_thresholds(all_data, .x, n = 100))

boot_summary_table <- boot_thresholds_all %>%
  group_by(veg_type) %>%
  summarise(
    early_mean  = mean(early_precip, na.rm = TRUE),
    early_ci_lo = quantile(early_precip, 0.025, na.rm = TRUE),
    early_ci_hi = quantile(early_precip, 0.975, na.rm = TRUE),
    
    late_mean   = mean(late_precip, na.rm = TRUE),
    late_ci_lo  = quantile(late_precip, 0.025, na.rm = TRUE),
    late_ci_hi  = quantile(late_precip, 0.975, na.rm = TRUE),
    
    n_points    = n()
  )


ggplot(boot_thresholds_all, aes(x = early_precip, y = late_precip)) +
  geom_point(alpha = 0.5, size = 1, shape=20) +
  stat_density2d(aes(fill = ..level..), geom = "polygon", bins = 10, alpha = 0.4) +
  facet_wrap(~veg_type) +
  labs(
    title = "Bootstrapped Joint Thresholds: Early vs Late Precipitation",
    subtitle = "Region = combinations yielding ≥ 75% of max forage production",
    x = "Early Precipitation (Jan–Jun, mm)",
    y = "Late Precipitation (Jul–Oct, mm)"
  ) +
  theme_minimal()

ggplot(boot_thresholds_all, aes(x = early_precip, y = late_precip)) +
  stat_bin_2d(bins = 49, aes(fill = after_stat(count))) +
  scale_fill_viridis_c(name = "Point Count", option = "C") +
  facet_wrap(~veg_type) +
  labs(
    title = "Bootstrapped Joint Thresholds: Early vs Late Precipitation",
    subtitle = "Region = combinations yielding ≥ 75% of max forage production",
    x = "Early Precipitation (Jan–Jun, mm)",
    y = "Late Precipitation (Jul–Oct, mm)"
  ) +
  theme_minimal()


joint_threshold_summary <- boot_thresholds_all %>%
  group_by(veg_type) %>%
  summarise(
    early_mean  = mean(early_precip, na.rm = TRUE),
    early_lo    = quantile(early_precip, 0.025, na.rm = TRUE),
    early_hi    = quantile(early_precip, 0.975, na.rm = TRUE),
    
    late_mean   = mean(late_precip, na.rm = TRUE),
    late_lo     = quantile(late_precip, 0.025, na.rm = TRUE),
    late_hi     = quantile(late_precip, 0.975, na.rm = TRUE),
    
    n_points    = n()
  ) %>%
  arrange(veg_type)

joint_threshold_summary <- joint_threshold_summary %>%
  mutate(
    season_bias = early_mean / (early_mean + late_mean)
  )

joint_threshold_summary <- joint_threshold_summary %>%
  mutate(
    seasonal_strategy = case_when(
      season_bias >= 0.7 ~ "early-dominated",
      season_bias <= 0.3 ~ "late-dominated",
      TRUE ~ "mixed"
    )
  )

ggplot(joint_threshold_summary, aes(x = early_mean, y = late_mean, label = veg_type)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = early_lo, xmax = early_hi), height = 0.1) +
  geom_errorbar(aes(ymin = late_lo, ymax = late_hi), width = 0.1) +
  geom_text(nudge_y = 10, size = 3) +
  theme_minimal() +
  labs(x = "Early Precip Threshold (mm)", y = "Late Precip Threshold (mm)",
       title = "Joint Precipitation Thresholds by Vegetation Type",
       subtitle = "75% Production Threshold with 95% Confidence Ranges")

# Visualize the GAM Smooths
mod_joint_bytype <- gam(veg_prod ~ te(early_precip, late_precip, by = veg_type) + veg_type,
                        data = all_data, family = Gamma(link = "log"), method = "REML")

veg_levels <- unique(all_data$veg_type)

grid_all <- do.call(rbind, lapply(veg_levels, function(vt) {
  grid <- expand.grid(
    early_precip = seq(min(all_data$early_precip), max(all_data$early_precip), length.out = 50),
    late_precip = seq(min(all_data$late_precip), max(all_data$late_precip), length.out = 50),
    veg_type = vt
  )
  grid$predicted <- predict(mod_joint_bytype, newdata = grid, type = "response")
  return(grid)
}))

grid_all <- grid_all %>%
  group_by(veg_type) %>%
  mutate(
    max_prod = max(predicted, na.rm = TRUE),
    pct_max = predicted / max_prod
  ) %>%
  ungroup()

ggplot(grid_all, aes(x = early_precip, y = late_precip, fill = predicted)) +
  geom_tile() +
  geom_contour(aes(z = pct_max),
               breaks = 0.75,   # your threshold (80%)
               color = "white",
               linewidth = 0.6) +
  # geom_point(data = boot_thresholds_all,
  #            aes(x = early_precip, y = late_precip),
  #            inherit.aes = FALSE,
  #            color = "black", alpha = 0.1, size = 0.3) +
  facet_wrap(~veg_type) +
  scale_fill_viridis_c(
    limits = c(0, 3000),
    oob = scales::squish,
    name = "lbs/ac"
  ) +
  labs(
    title = "GAM Predicted Forage Production by Vegetation Type",
    subtitle = "White contour = 75% of max production",
    x = "Early Precip (mm)", y = "Late Precip (mm)"
  ) +
  theme_minimal()


pred <- predict(mod_joint_bytype, type = "response")
ggplot(all_data, aes(x = pred, y = veg_prod)) +
  geom_point(alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  coord_equal() +
  labs(title = "Predicted vs. Observed Forage Production",
       x = "Predicted", y = "Observed")


##### OPTIMAL precip for max production
optima_table <- grid_all %>%
  group_by(veg_type) %>%
  slice_max(predicted, n = 1, with_ties = FALSE) %>%
  ungroup()

ggplot(grid_all, aes(x = early_precip, y = late_precip, fill = predicted)) +
  geom_tile() +
  facet_wrap(~veg_type) +
  scale_fill_viridis_c(name = "Predicted\nForage (lbs/ac)") +
  geom_point(data = optima_table,
             aes(x = early_precip, y = late_precip),
             color = "red", size = 2, shape = 4, stroke = 1.5,
             inherit.aes = FALSE) +
  scale_fill_viridis_c(
    limits = c(0, 3000),
    oob = scales::squish,
    name = "lbs/ac"
  ) +
  labs(
    title = "GAM Predicted Forage Production with Ecological Optima",
    subtitle = "Red X = (early, late) combo yielding max predicted forage",
    x = "Early Precip (Jan–Jun, mm)", y = "Late Precip (Jul–Oct, mm)"
  ) +
  theme_minimal()

optima_table %>%
  select(veg_type, early_precip, late_precip, predicted) %>%
  arrange(desc(predicted))


ggplot(boot_thresholds_all, aes(x = early_precip, y = late_precip)) +
  geom_point(alpha = 0.15, size = 0.8, shape = 20, color = "black") +  # lower opacity to avoid masking
  stat_density2d(
    aes(fill = ..level..),
    geom = "polygon",
    bins = 10,
    alpha = 0.6,
    color = NA  # turn off contour lines to reduce clutter
  ) +
  scale_fill_viridis_c(name = "Density", option = "C") +
  facet_wrap(~veg_type) +
  labs(
    title = "Bootstrapped Joint Thresholds: Early vs Late Precipitation",
    subtitle = "Shaded regions show combinations yielding ≥ 90% of max forage production",
    x = "Early Precipitation (Jan–Jun, mm)",
    y = "Late Precipitation (Jul–Oct, mm)"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 10, face = "bold"),
    plot.title = element_text(face = "bold")
  )+
  geom_point(data = optima_table,
               aes(x = early_precip, y = late_precip),
               shape = 4, color = "red", size = 2, stroke = 1,
               inherit.aes = FALSE)


