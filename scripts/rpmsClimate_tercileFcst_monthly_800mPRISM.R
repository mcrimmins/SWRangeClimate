# RPMS climate analysis for KNF
# using 800m PRISM data and EVT
# MAC 12/2/25
##### ORDINAL CATEGORICAL FORECAST 

library(terra)
library(ggplot2)
library(lubridate)
library(RColorBrewer)
library(dplyr)

terraOptions(progress=1)

# Load SpatRaster data -----
#veg_prod <- rast("./data/processed/800m_KNF_RPMS_1984_2024.tif")  # Annual vegetation production
veg_prod <- rast("./data/processed/800m_KNF_RPMS_1984_2024_EVTmask.tif")  # Annual vegetation production with EVT mask
precip <- rast("./data/processed/KNF_PRISM_800m_monthly_pr_1984_2024.tif")  # Monthly precipitation
#veg_type <- rast("./data/processed/800m_KNF_TEU_PNV_CLASS.tif")  # Vegetation type
veg_type <- rast("./data/processed/800m_KNF_LF_EVT_treeMask_top5_2024.tif")  # Vegetation type


# ============================================================
# OPTIONAL: FAST TESTING MODE — crop rasters to small region
# ============================================================

testing <- FALSE     # <---- set to FALSE when running full dataset

if (testing) {
  message("Testing mode ON: Cropping rasters to small test window...")
  
  # Choose a small bounding box inside the Kaibab NF area
  # Modify coordinates if you want a different region
#  test_extent <- ext(-112.5, -112.0,   # xmin, xmax
#                     35.9,  36.2)       # ymin, ymax
  test_extent <- ext(-112.3, -112.2, 36.0, 36.1) ### very small sample
  
  
  # Crop all rasters consistently
  veg_prod <- crop(veg_prod, test_extent)
  precip   <- crop(precip,   test_extent)
  veg_type <- crop(veg_type, test_extent)
  
  message("Test rasters created: ",
          ncell(veg_prod), " pixels instead of ", ncell(rast("./data/processed/800m_KNF_RPMS_1984_2024_EVTmask.tif")),
          " (full).")
}
#########

###############################################################
# OPTIONAL SENSOR-SHIFT CORRECTION + BREAKPOINT DIAGNOSTICS
# Supports: # (1) auto_correction = TRUE → match early/late via mean differences 
# (2) manual_factor = 0.87 → reduce late period by 13%
###############################################################

library(terra)
library(strucchange)

veg_input <- veg_prod

apply_correction <- TRUE # master switch
auto_correction  <- TRUE # TRUE = use automatic shift; FALSE = use manual factor
manual_factor    <- 0.87 # e.g., 0.87 for 13% test; 1.00 = no manual correction

###############################################################
# STEP 1 — Detect breakpoint BEFORE applying correction
###############################################################

years_vec <- as.numeric(names(veg_input))
veg_mean_raw <- global(veg_input, fun="mean", na.rm=TRUE)[,1]

bp_model_raw <- breakpoints(veg_mean_raw ~ 1)
break_index  <- bp_model_raw$breakpoints[1]
break_year   <- years_vec[break_index]

message("Detected breakpoint around year: ", break_year)

###############################################################
# STEP 2 — Apply correction (if enabled)
###############################################################

if (apply_correction) {
  
  message("Applying vegetation sensor-shift correction...")
  
  years <- years_vec
  early_period <- years <= break_year
  late_period  <- years >  break_year
  
  veg_corrected <- veg_input
  
  if (auto_correction) {
    
    message("Using AUTOMATIC correction based on early/late mean differences...")
    
    early_mean <- app(veg_input[[early_period]], mean, na.rm = TRUE)
    late_mean  <- app(veg_input[[late_period]],  mean, na.rm = TRUE)
    
    shift_amount <- early_mean - late_mean
    
    for (i in which(late_period)) {
      veg_corrected[[i]] <- veg_input[[i]] + shift_amount
    }
    
    message("Automatic correction applied.")
  }
  
  if (!auto_correction && manual_factor != 1) {
    
    message("Using MANUAL correction. Factor = ", manual_factor)
    
    for (i in which(late_period)) {
      veg_corrected[[i]] <- veg_input[[i]] * manual_factor
    }
    
    message("Manual correction applied.")
  }
  
  veg_used <- veg_corrected
  
} else {
  
  message("No correction applied.")
  veg_used <- veg_input
}

###############################################################
# STEP 3 — DIAGNOSTICS (using corrected or raw data)
###############################################################

years <- as.numeric(names(veg_used))
veg_mean <- global(veg_used, fun="mean", na.rm=TRUE)[,1]

ts_df <- data.frame(
  year = years,
  mean_prod = veg_mean
)

early_mean_global <- mean(veg_mean[years <= break_year])
late_mean_global  <- mean(veg_mean[years >  break_year])

bp_model <- breakpoints(mean_prod ~ 1, data = ts_df)
break_pos <- bp_model$breakpoints
break_years <- years[break_pos]

cat("Estimated Break Year(s):", break_years, "\n")

# -------------------------------------------------------------
# 4. Plot time series with breakpoints + segmented means
# -------------------------------------------------------------
plot(ts_df$year, ts_df$mean_prod,
     type="b", pch=16, col="black",
     xlab="Year", ylab="Mean Vegetation Production",
     main="Vegetation Production Time Series\nWith Detected Breakpoints")

# Vertical lines at breakpoints
abline(v = break_years, col="red", lwd=3)

# Stepwise mean model
seg_means <- fitted(bp_model)
lines(ts_df$year, seg_means, col="blue", lwd=3)

legend("topleft",
       legend=c("Observed Mean", "Segmented Mean", "Breakpoint"),
       col=c("black", "blue", "red"),
       lwd=c(1,3,3), pch=c(16, NA, NA))

# -------------------------------------------------------------
# 5. Confidence intervals for breakpoints
# -------------------------------------------------------------
print(confint(bp_model))

# reassign veg_prod to the used (corrected or raw) data
veg_prod <- veg_used

###############################################################


#########
# Use vegetation raster because it is guaranteed to match the pixel grid
base <- veg_prod[[1]]

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
#####


## ============================================================
## 0. Packages
## ============================================================
library(ordinal)   # for clm

## ============================================================
## 1. Harmonize precip to veg_prod grid
## ============================================================
precip <- project(precip, veg_prod, method = "bilinear")
precip <- crop(precip, veg_prod)
#precip <- align(precip, veg_prod)

## ============================================================
## 2. Build WY-month index (Oct = 1, …, Sep = 12)
## ============================================================
precip_dates <- as.Date(paste0(names(precip), "01"), "%Y%m%d")
yr  <- as.integer(format(precip_dates, "%Y"))
mon <- as.integer(format(precip_dates, "%m"))

# WY Year
wy_year <- ifelse(mon >= 10, yr + 1L, yr)

# WY Month: Oct=1, …, Sep=12
wy_month <- mon - 9L
wy_month[wy_month <= 0] <- wy_month[wy_month <= 0] + 12L

## ============================================================
## 3. Define cumulative windows (Dec, Mar, Jun, Sep)
## ============================================================
cum_windows <- list(
  Dec = 3,   # Oct–Dec cumulative
  Mar = 6,   # Oct–Mar cumulative
  Jun = 9,   # Oct–Jun cumulative
  Sep = 12   # Oct–Sep cumulative
)

## ============================================================
## 4. Safe cumulative-sum builder
##    Skips years lacking early WY months
## ============================================================
make_cum_stack_safe <- function(precip, wy_year, wy_month, m_max) {
  out_list <- list()
  out_years <- c()
  
  for (y in sort(unique(wy_year))) {
    idx <- wy_year == y & wy_month <= m_max
    if (!any(idx)) {
      message("Skipping WY ", y, ": missing months for m_max=", m_max)
      next
    }
    # sum raster over all months up to m_max
    out_list[[length(out_list) + 1]] <- sum(precip[[idx]], na.rm = TRUE)
    out_years <- c(out_years, y)
  }
  
  r <- rast(out_list)
  names(r) <- out_years
  r
}

## ============================================================
## 5. Build cumulative precip stacks for each window
## ============================================================
Pcum_all <- lapply(cum_windows, function(m) {
  make_cum_stack_safe(precip, wy_year, wy_month, m_max = m)
})

# save Pcum_all to disk
#saveRDS(Pcum_all, file = "./data/processed/KNF_PRISM_800m_cumulative_precip_stacks.rds")
#Pcum_all <- readRDS("~/RProjects/SWRangeClimate/data/processed/KNF_PRISM_800m_cumulative_precip_stacks.rds")

Pcum_Dec <- Pcum_all$Dec
Pcum_Mar <- Pcum_all$Mar
Pcum_Jun <- Pcum_all$Jun
Pcum_Sep <- Pcum_all$Sep

## ============================================================
## 6. Align all cumulative stacks to common usable years
## ============================================================
## ============================================================
## FIXED YEAR ALIGNMENT LOGIC
## ============================================================

# 1. Identify the years present in each cumulative precipitation stack
yrs_dec <- names(Pcum_Dec)
yrs_mar <- names(Pcum_Mar)
yrs_jun <- names(Pcum_Jun)
yrs_sep <- names(Pcum_Sep)
yrs_veg <- names(veg_prod)

# Convert to integers for safe comparison
yrs_dec <- as.integer(yrs_dec)
yrs_mar <- as.integer(yrs_mar)
yrs_jun <- as.integer(yrs_jun)
yrs_sep <- as.integer(yrs_sep)
yrs_veg <- as.integer(yrs_veg)

# 2. Find years available in ALL datasets
common_years <- Reduce(intersect, list(
  yrs_dec, yrs_mar, yrs_jun, yrs_sep, yrs_veg
))

common_years <- sort(common_years)

message("COMMON YEARS: ", paste(common_years, collapse = ", "))

# 3. Subset rasters to common years
Pcum_Dec <- Pcum_Dec[[as.character(common_years)]]
Pcum_Mar <- Pcum_Mar[[as.character(common_years)]]
Pcum_Jun <- Pcum_Jun[[as.character(common_years)]]
Pcum_Sep <- Pcum_Sep[[as.character(common_years)]]

veg <- veg_prod[[as.character(common_years)]]

# 4. Final consistency check
stopifnot(
  identical(names(Pcum_Dec), names(veg)),
  identical(names(Pcum_Mar), names(veg)),
  identical(names(Pcum_Jun), names(veg)),
  identical(names(Pcum_Sep), names(veg))
)

message("All stacks aligned. Years used: ",
        paste(names(veg), collapse = ", "))

veg_yrs <- as.character((as.integer(names(veg))))

## ============================================================
## 4. Production terciles per pixel (1=below,2=near,3=above)
## ============================================================

make_tercile_cat <- function(x, probs = c(1/3, 2/3)) {
  if (all(is.na(x))) return(rep(NA_integer_, length(x)))
  qs <- quantile(x, probs = probs, na.rm = TRUE, type = 8)
  cut(x, breaks = c(-Inf, qs[1], qs[2], Inf), labels = FALSE)
}

veg_terc <- app(veg, make_tercile_cat)
names(veg_terc) <- veg_yrs

## ============================================================
## Per-pixel ordinal modeling via two GLMs + progress bar
## ============================================================

library(terra)
library(future.apply)
library(compiler)
library(progressr)

## 1. Ensure rasters are in memory and aligned
# 
# veg_terc   <- readAll(veg_terc)
# Pcum_Dec   <- readAll(Pcum_Dec)
# Pcum_Mar   <- readAll(Pcum_Mar)
# Pcum_Jun   <- readAll(Pcum_Jun)
# Pcum_Sep   <- readAll(Pcum_Sep)

terraOptions(tempdir = "/dev/shm")   # RAM-based scratch on Linux

# Matrices: rows = pixels, cols = years
Y_mat      <- as.matrix(values(veg_terc))     # tercile categories 1..3
P_Dec_mat  <- as.matrix(values(Pcum_Dec))     # cumulative Oct–Dec (mm)
P_Mar_mat  <- as.matrix(values(Pcum_Mar))     # cumulative Oct–Mar
P_Jun_mat  <- as.matrix(values(Pcum_Jun))     # cumulative Oct–Jun
P_Sep_mat  <- as.matrix(values(Pcum_Sep))     # cumulative Oct–Sep

ncell <- nrow(Y_mat)

## 2. Per-pixel GLM-based ordinal approximation
##    Model 1: I(prod_terc <= 1)  (below vs near/above)
##    Model 2: I(prod_terc <= 2)  (below/near vs above)
##
##    Thresholds (mm) at P = 0.5 are:
##      th1 = -a1 / b1   (boundary below vs near/above)
##      th2 = -a2 / b2   (boundary below/near vs above)

fit_pixel_glm_ord <- function(y, x) {
  ok <- !is.na(y) & !is.na(x)
  if (sum(ok) < 20) return(NULL)  # require some years
  
  y_ok <- y[ok]
  x_ok <- x[ok]
  
  # Model 1: below (1) vs others (0)
  z1 <- as.integer(y_ok <= 1L)
  if (length(unique(z1)) < 2) return(NULL)  # need variation
  
  df1 <- data.frame(z = z1, Pcum = x_ok)
  m1 <- try(glm(z ~ Pcum, data = df1, family = binomial), silent = TRUE)
  if (inherits(m1, "try-error")) return(NULL)
  co1 <- coef(m1)
  a1  <- unname(co1["(Intercept)"])
  b1  <- unname(co1["Pcum"])
  if (is.na(a1) || is.na(b1) || b1 == 0) return(NULL)
  th1 <- -a1 / b1   # mm at which P(z1=1) = 0.5
  
  # Model 2: below or near (1/2) vs above (3)
  z2 <- as.integer(y_ok <= 2L)
  if (length(unique(z2)) < 2) return(NULL)
  
  df2 <- data.frame(z = z2, Pcum = x_ok)
  m2 <- try(glm(z ~ Pcum, data = df2, family = binomial), silent = TRUE)
  if (inherits(m2, "try-error")) return(NULL)
  co2 <- coef(m2)
  a2  <- unname(co2["(Intercept)"])
  b2  <- unname(co2["Pcum"])
  if (is.na(a2) || is.na(b2) || b2 == 0) return(NULL)
  th2 <- -a2 / b2   # mm at which P(z2=1) = 0.5
  
  # simple sensitivity metric (average absolute slope)
  sens <- mean(c(abs(b1), abs(b2)))
  
  list(th1 = th1, th2 = th2, b1 = b1, b2 = b2, sens = sens)
}

fit_pixel_glm_ord_cmp <- compiler::cmpfun(fit_pixel_glm_ord)

## 3. Parallel setup
##    Adjust workers to VM cores (e.g., 4–8)

plan(multisession, workers = 6)

## 4. Progress bar + parallel loop across pixels and lead times

handlers(global = TRUE)

with_progress({
  p <- progressor(steps = ncell)
  
  pixel_results <- future_lapply(
    seq_len(ncell),
    function(i) {
      terraOptions(tempdir = "/dev/shm")
      
      y_i <- Y_mat[i, ]
      if (all(is.na(y_i))) {
        p()
        return(list(Dec = NULL, Mar = NULL, Jun = NULL, Sep = NULL))
      }
      
      res_Dec <- fit_pixel_glm_ord_cmp(y_i, P_Dec_mat[i, ])
      res_Mar <- fit_pixel_glm_ord_cmp(y_i, P_Mar_mat[i, ])
      res_Jun <- fit_pixel_glm_ord_cmp(y_i, P_Jun_mat[i, ])
      res_Sep <- fit_pixel_glm_ord_cmp(y_i, P_Sep_mat[i, ])
      
      p()
      
      list(Dec = res_Dec, Mar = res_Mar, Jun = res_Jun, Sep = res_Sep)
    },
    future.seed = TRUE
  )
})

plan(sequential)

## 5. Unpack results into vectors

th1_Dec <- th2_Dec <- sens_Dec <- rep(NA_real_, ncell)
th1_Mar <- th2_Mar <- sens_Mar <- rep(NA_real_, ncell)
th1_Jun <- th2_Jun <- sens_Jun <- rep(NA_real_, ncell)
th1_Sep <- th2_Sep <- sens_Sep <- rep(NA_real_, ncell)

for (i in seq_len(ncell)) {
  res_i <- pixel_results[[i]]
  # Dec
  if (!is.null(res_i$Dec)) {
    th1_Dec[i]  <- res_i$Dec$th1
    th2_Dec[i]  <- res_i$Dec$th2
    sens_Dec[i] <- res_i$Dec$sens
  }
  # Mar
  if (!is.null(res_i$Mar)) {
    th1_Mar[i]  <- res_i$Mar$th1
    th2_Mar[i]  <- res_i$Mar$th2
    sens_Mar[i] <- res_i$Mar$sens
  }
  # Jun
  if (!is.null(res_i$Jun)) {
    th1_Jun[i]  <- res_i$Jun$th1
    th2_Jun[i]  <- res_i$Jun$th2
    sens_Jun[i] <- res_i$Jun$sens
  }
  # Sep
  if (!is.null(res_i$Sep)) {
    th1_Sep[i]  <- res_i$Sep$th1
    th2_Sep[i]  <- res_i$Sep$th2
    sens_Sep[i] <- res_i$Sep$sens
  }
}

## 6. Quick checks

summary(th1_Dec)
summary(th2_Dec)
summary(th1_Sep)
summary(th2_Sep)
summary(sens_Sep)


## ------------------------------------------------------------
## Sensitivity thresholds
## ------------------------------------------------------------
## ============================================================
## POSTPROCESSING: thresholds, sensitivity, masks, veg summaries,
##                response-type classification
## ============================================================

library(terra)
library(dplyr)

## ------------------------------------------------------------
## 0. Helper to build rasters from vectors
## ------------------------------------------------------------
mk_r <- function(vals, name) {
  r <- base
  values(r) <- vals
  names(r) <- name
  r
}

## ------------------------------------------------------------
## 1. Filtering thresholds & creating clean rasters
## ------------------------------------------------------------

# Sensitivity and threshold filters
sens_cut <- 0.005     # minimum sensitivity for meaningful response
th_min   <- 0         # mm
th_max   <- 2000     # mm (adjust if needed for your climate)

clean_thresholds <- function(th1, th2, sens) {
  good <- !is.na(th1) & !is.na(th2) & !is.na(sens) &
    th1 > th_min & th1 < th_max &
    th2 > th_min & th2 < th_max &
    th2 > th1 &             # enforce proper ordering
    sens > sens_cut         # enforce meaningful response
  
  list(
    th1_clean = ifelse(good, th1, NA_real_),
    th2_clean = ifelse(good, th2, NA_real_),
    mask      = good
  )
}

Dec_clean <- clean_thresholds(th1_Dec, th2_Dec, sens_Dec)
Mar_clean <- clean_thresholds(th1_Mar, th2_Mar, sens_Mar)
Jun_clean <- clean_thresholds(th1_Jun, th2_Jun, sens_Jun)
Sep_clean <- clean_thresholds(th1_Sep, th2_Sep, sens_Sep)

# Cleaned threshold rasters (mm)
th_Dec_bn_r <- mk_r(Dec_clean$th1_clean, "P_Dec_below_near_mm")
th_Dec_na_r <- mk_r(Dec_clean$th2_clean, "P_Dec_near_above_mm")

th_Mar_bn_r <- mk_r(Mar_clean$th1_clean, "P_Mar_below_near_mm")
th_Mar_na_r <- mk_r(Mar_clean$th2_clean, "P_Mar_near_above_mm")

th_Jun_bn_r <- mk_r(Jun_clean$th1_clean, "P_Jun_below_near_mm")
th_Jun_na_r <- mk_r(Jun_clean$th2_clean, "P_Jun_near_above_mm")

th_Sep_bn_r <- mk_r(Sep_clean$th1_clean, "P_Sep_below_near_mm")
th_Sep_na_r <- mk_r(Sep_clean$th2_clean, "P_Sep_near_above_mm")

## Quick sanity checks
summary(values(th_Dec_bn_r))
summary(values(th_Sep_na_r))

## ============================================================
##  Construct slope rasters (b1 and b2) for each season
## ============================================================

make_slope_raster <- function(bvec, mask, name) {
  vals <- ifelse(mask, bvec, NA_real_)
  r <- base
  values(r) <- vals
  names(r) <- name
  r
}

# Helper function
make_slope_raster <- function(bvec, mask, name) {
  r <- base
  vals <- ifelse(mask, bvec, NA_real_)
  values(r) <- vals
  names(r) <- name
  r
}

# Extract raw slope vectors from pixel_results
b1_Dec <- sapply(pixel_results, function(x) if (is.null(x$Dec)) NA else x$Dec$b1)
b2_Dec <- sapply(pixel_results, function(x) if (is.null(x$Dec)) NA else x$Dec$b2)

b1_Mar <- sapply(pixel_results, function(x) if (is.null(x$Mar)) NA else x$Mar$b1)
b2_Mar <- sapply(pixel_results, function(x) if (is.null(x$Mar)) NA else x$Mar$b2)

b1_Jun <- sapply(pixel_results, function(x) if (is.null(x$Jun)) NA else x$Jun$b1)
b2_Jun <- sapply(pixel_results, function(x) if (is.null(x$Jun)) NA else x$Jun$b2)

b1_Sep <- sapply(pixel_results, function(x) if (is.null(x$Sep)) NA else x$Sep$b1)
b2_Sep <- sapply(pixel_results, function(x) if (is.null(x$Sep)) NA else x$Sep$b2)

# Build slope rasters for each season (filtered by their masks)
b1_Dec_r <- make_slope_raster(b1_Dec, Dec_clean$mask, "b1_Dec")
b2_Dec_r <- make_slope_raster(b2_Dec, Dec_clean$mask, "b2_Dec")

b1_Mar_r <- make_slope_raster(b1_Mar, Mar_clean$mask, "b1_Mar")
b2_Mar_r <- make_slope_raster(b2_Mar, Mar_clean$mask, "b2_Mar")

b1_Jun_r <- make_slope_raster(b1_Jun, Jun_clean$mask, "b1_Jun")
b2_Jun_r <- make_slope_raster(b2_Jun, Jun_clean$mask, "b2_Jun")

b1_Sep_r <- make_slope_raster(b1_Sep, Sep_clean$mask, "b1_Sep")
b2_Sep_r <- make_slope_raster(b2_Sep, Sep_clean$mask, "b2_Sep")

# Sanity check
plot(b1_Sep_r, main="Slope b1 (Sep)")
plot(b2_Sep_r, main="Slope b2 (Sep)")

## ------------------------------------------------------------
## 2. Predictability masks
## ------------------------------------------------------------

mask_Dec <- mk_r(Dec_clean$mask, "Dec_predictable")
mask_Mar <- mk_r(Mar_clean$mask, "Mar_predictable")
mask_Jun <- mk_r(Jun_clean$mask, "Jun_predictable")
mask_Sep <- mk_r(Sep_clean$mask, "Sep_predictable")

# Pixels predictable in at least one season
mask_any <- mk_r(
  Dec_clean$mask | Mar_clean$mask | Jun_clean$mask | Sep_clean$mask,
  "predictable_any"
)

# Pixels predictable in all four seasons
mask_all <- mk_r(
  Dec_clean$mask & Mar_clean$mask & Jun_clean$mask & Sep_clean$mask,
  "predictable_all"
)

## ------------------------------------------------------------
## 3. Sensitivity rasters (raw & filtered)
## ------------------------------------------------------------

# Raw sensitivity rasters
sens_Dec_r <- mk_r(sens_Dec, "sens_Dec")
sens_Mar_r <- mk_r(sens_Mar, "sens_Mar")
sens_Jun_r <- mk_r(sens_Jun, "sens_Jun")
sens_Sep_r <- mk_r(sens_Sep, "sens_Sep")

# Filtered sensitivity (only where thresholds are valid)
sens_Dec_clean <- ifelse(Dec_clean$mask, sens_Dec, NA_real_)
sens_Mar_clean <- ifelse(Mar_clean$mask, sens_Mar, NA_real_)
sens_Jun_clean <- ifelse(Jun_clean$mask, sens_Jun, NA_real_)
sens_Sep_clean <- ifelse(Sep_clean$mask, sens_Sep, NA_real_)

sens_Dec_rf <- mk_r(sens_Dec_clean, "sens_Dec_filtered")
sens_Mar_rf <- mk_r(sens_Mar_clean, "sens_Mar_filtered")
sens_Jun_rf <- mk_r(sens_Jun_clean, "sens_Jun_filtered")
sens_Sep_rf <- mk_r(sens_Sep_clean, "sens_Sep_filtered")

## ------------------------------------------------------------
## 4. Vegetation-type summaries
## ------------------------------------------------------------

# Ensure veg_type is on the same grid; if already aligned, skip project()
# veg_type <- project(veg_type, base, method = "near")

v_veg   <- values(veg_type)
v_Dec_t1 <- values(th_Dec_bn_r)
v_Dec_t2 <- values(th_Dec_na_r)
v_Sep_t1 <- values(th_Sep_bn_r)
v_Sep_t2 <- values(th_Sep_na_r)
v_sens_Sep <- values(sens_Sep_rf)

df_summary <- data.frame(
  veg_type = v_veg,
  th1_Dec  = v_Dec_t1,
  th2_Dec  = v_Dec_t2,
  th1_Sep  = v_Sep_t1,
  th2_Sep  = v_Sep_t2,
  sens_Sep = v_sens_Sep
)

df_summary <- df_summary[!is.na(df_summary$EVT_NAME), ]
colnames(df_summary)[1] <- "veg_type"

veg_type_summary <- df_summary %>%
  group_by(veg_type) %>%
  summarise(
    npix             = n(),
    pct_predictable  = mean(!is.na(th1_Dec) | !is.na(th1_Sep)) * 100,
    Dec_th1_median   = median(th1_Dec, na.rm = TRUE),
    Dec_th2_median   = median(th2_Dec, na.rm = TRUE),
    Sep_th1_median   = median(th1_Sep, na.rm = TRUE),
    Sep_th2_median   = median(th2_Sep, na.rm = TRUE),
    Sep_sens_median  = median(sens_Sep, na.rm = TRUE)
  ) %>%
  arrange(desc(pct_predictable))

# View summary table
veg_type_summary

library(ggplot2)

ggplot(df_summary, aes(x = as.factor(veg_type), y = P_Sep_below_near_mm, group=veg_type)) +
  geom_boxplot(outlier.size = 0.5, varwidth = TRUE) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Sep threshold (below→near) by vegetation type")

ggplot(df_summary, aes(x = as.factor(veg_type), y = P_Sep_near_above_mm, group=veg_type)) +
  geom_boxplot(outlier.size = 0.5, varwidth = TRUE) +
  theme(axis.text.x = element_text(angle = 90)) +
  labs(title = "Sep threshold (near→above) by vegetation type")


## ------------------------------------------------------------
## 5. Response-type classification
## ------------------------------------------------------------

# Logical masks as vectors
mD <- Dec_clean$mask
mM <- Mar_clean$mask
mJ <- Jun_clean$mask
mS <- Sep_clean$mask

classify_pixel <- function(d, m, j, s) {
  # All FALSE or NA → weak/insensitive
  if ((is.na(d) | !d) & (is.na(m) | !m) & (is.na(j) | !j) & (is.na(s) | !s)) {
    return(0L)  # weak
  }
  # Monsoon-dominated: only late-season predictability
  if ((is.na(d) | !d) & (is.na(m) | !m) & s) {
    return(1L)  # monsoon-driven
  }
  # Winter-dominated: early-season predictability, weak later
  if ((d | m) & (is.na(s) | !s) & (is.na(j) | !j | !s)) {
    return(2L)  # winter-driven
  }
  # Strong at both winter & summer
  if ((d | m) & s) {
    return(3L)  # dual-control (winter + summer)
  }
  # Everything else
  return(4L)    # other / mixed
}

resp_type <- mapply(classify_pixel, mD, mM, mJ, mS)

resp_type_r <- mk_r(resp_type, "response_type")

# Legend (document for yourself):
# 0 = weak / insensitive
# 1 = monsoon-driven
# 2 = winter-driven
# 3 = dual-control (winter + summer)
# 4 = other / mixed

## ------------------------------------------------------------
## 6. Quick plots (optional, for checking)
## ------------------------------------------------------------

plot(th_Sep_na_r, main = "Oct–Sep near→above threshold (mm)", range=c(0,1000))
  plot(shp, add=TRUE, color=NA)
plot(th_Dec_na_r, main = "Dec–Sep near→above threshold (mm)", range=c(0,500))
  plot(shp, add=TRUE, color=NA)
  
plot(sens_Sep_rf, main = "Sensitivity to Oct–Sep precip (filtered)")
plot(resp_type_r, main = "Response Type")
plot(mask_any, main = "Predictable (any season)")
plot(mask_all, main = "Predictable (all four seasons)")




# predictable in any season vs all seasons
plot(mask_any, 
     main = "Pixels with Predictable Response (Any Season)",
     col = c("gray80","darkgreen"))
plot(shp, add=TRUE, color=NA)


plot(mask_all, 
     main = "Pixels Predictable in All Four Seasons",
     col = c("gray90","navy"))
plot(shp, add=TRUE, color=NA)

df_t <- data.frame(
  veg_type = values(veg_type),
  th_Sep_bn = values(th_Sep_bn_r),
  th_Sep_na = values(th_Sep_na_r),
  th_Dec_bn = values(th_Dec_bn_r),
  th_Dec_na = values(th_Dec_na_r),
  sens_Sep  = values(sens_Sep_rf)
)
df_t <- df_t[complete.cases(df_t), ]

ggplot(df_t, aes(x = P_Sep_below_near_mm)) +
  geom_histogram(bins = 50, fill = "steelblue") +
  labs(title = "Distribution of Sep Threshold (Below→Near)", x = "mm")

ggplot(df_t, aes(x = P_Sep_near_above_mm)) +
  geom_histogram(bins = 50, fill = "darkorange") +
  labs(title = "Distribution of Sep Threshold (Near→Above)", x = "mm")

# 2B. Uncertainty via width of transition zone
# The “transition humidity window” indicates robustness of category transitions.
df_t$transition_width <- df_t$P_Sep_near_above_mm - df_t$P_Sep_below_near_mm

ggplot(df_t, aes(x = transition_width)) +
  geom_histogram(bins = 50, fill = "purple") +
  labs(title = "Width of Precipitation Transition Window (mm)",
       x = "th2 - th1")

# 2C. Scatter of sensitivity vs transition window
# This shows predictive sharpness.
# Interpretation:
# High sensitivity + narrow window = strong predictability
# Low sensitivity + wide window = weak predictability
ggplot(df_t, aes(x = sens_Sep_filtered, y = transition_width)) +
  geom_point(alpha = 0.2, size = 0.7) +
  scale_x_log10() +
  labs(title = "Sensitivity vs Transition Width",
       x = "Sensitivity (slope)", y = "th2 - th1 (mm)")

# spatial uncertainty maps
# 3A. Sensitivity Map (Interpretation: Slope of GLM)
# Shows where production responds steeply vs weakly to precipitation.
plot(sens_Sep_rf,
     main = "Sensitivity to Summer Precipitation",
     col = hcl.colors(50, "Viridis"), range=c(0,0.05))
plot(shp, add=TRUE, color=NA)

plot(sens_Mar_rf,
     main = "Sensitivity to Oct-Mar Precipitation",
     col = hcl.colors(50, "Viridis"), range=c(0,0.05))
plot(shp, add=TRUE, color=NA)

# 3B. Transition Width Map
# This measures uncertainty spatially.
# Interpretation:
#   Large widths = highly uncertain or noisy transitions
# Small widths = sharp, predictable response

transition_width_r <- th_Sep_na_r - th_Sep_bn_r
names(transition_width_r) <- "transition_width"

plot(transition_width_r,
     main="Transition Width (mm): Precipitation Needed to Cross Categories",
     col=hcl.colors(50,"Inferno"), range=c(0,1000))
plot(shp, add=TRUE, color=NA)

df_u <- data.frame(
  Dec = values(th_Dec_na_r - th_Dec_bn_r),
  Mar = values(th_Mar_na_r - th_Mar_bn_r),
  Jun = values(th_Jun_na_r - th_Jun_bn_r),
  Sep = values(th_Sep_na_r - th_Sep_bn_r)
)

df_u <- df_u[complete.cases(df_u), ]

# 5A. Boxplots of transition widths by season
library(reshape2)
dfm <- melt(df_u)

ggplot(dfm, aes(x=variable, y=value)) +
  geom_boxplot() +
  labs(title="Seasonal Uncertainty (Transition Width)",
       x="Season", y="mm of Precipitation")

# Interpretation:
# Wide boxes = high uncertainty
# Short boxes = sharp seasonal thresholds
# Usually Sep has the tightest transitions in monsoon ecosystems

######
# 1. RASTER FORECAST PROBABILITIES FOR A GIVEN PRECIP AMOUNT

# Input precipitation forecast (mm)
P_forecast <- 100   # example: 200 mm forecast of Jun–Sep precip

# Recompute intercepts from slopes + thresholds
a1 <- -b1_Sep_r * th_Sep_bn_r
a2 <- -b2_Sep_r * th_Sep_na_r

names(a1) <- "a1_Sep"
names(a2) <- "a2_Sep"

names(b1_Sep_r) <- "b1_Sep"
names(b2_Sep_r) <- "b2_Sep"

# Cumulative probabilities
P_le_1 <- app(a1 + b1_Sep_r * P_forecast, plogis)
P_le_2 <- app(a2 + b2_Sep_r * P_forecast, plogis)

# Ensure ordinal ordering
P_le_2 <- app(c(P_le_1, P_le_2), fun = function(x) max(x[1], x[2]))

# Category probabilities
Pr_below <- P_le_1
Pr_near  <- P_le_2 - P_le_1
Pr_above <- 1 - P_le_2

# Clean up
Pr_below <- clamp(Pr_below, 0, 1)
Pr_near  <- clamp(Pr_near,  0, 1)
Pr_above <- clamp(Pr_above, 0, 1)

names(Pr_below) <- "Pr_below"
names(Pr_near)  <- "Pr_near"
names(Pr_above) <- "Pr_above"

# Plot probabilities
prob_stack <- c(Pr_below, Pr_near, Pr_above, Pr_below + Pr_near + Pr_above)
names(prob_stack) <- c("Below Normal", "Near Normal", "Above Normal", "Sum Check")

plot(prob_stack)


# Entropy Map (Information Uncertainty)
# ============================================================
# Entropy quantifies uncertainty:
#   H = − Σ p * log(p)
# H is near 0 → very certain
# H is high → very uncertain (probabilities close to equal)

entropy <- -(Pr_below*log(Pr_below) + 
Pr_near*log(Pr_near) + 
  Pr_above*log(Pr_above))

names(entropy) <- "forecast_entropy"

plot(entropy,
     main = "Forecast Uncertainty (Entropy)",
     col = hcl.colors(50, "Plasma"))

# 3. FORECAST UNCERTAINTY METRIC 2:
#   Category Spread (max – min probability)
# ============================================================
# This quantifies “how decisive the forecast is”.
# spread = max(Pr_below, Pr_near, Pr_above) − min(...)
# Spread = 1 → extremely confident
# Spread = 0 → very uncertain

Pr_stack <- c(Pr_below, Pr_near, Pr_above)

maxPr <- app(Pr_stack, fun=max)
minPr <- app(Pr_stack, fun=min)

spread <- maxPr - minPr
names(spread) <- "forecast_spread"

plot(spread,
     main = "Forecast Certainty (Probability Spread)",
     col = hcl.colors(50, "Viridis"))


# 4. FORECAST UNCERTAINTY METRIC 3:
#   Peak Probability (Forecast Confidence)
# ============================================================
#   
#   This is simply the highest probability for any category:

conf <- maxPr
names(conf) <- "forecast_confidence"

plot(conf,
     main = "Forecast Confidence (Peak Probability)",
     col = hcl.colors(50, "Sunset"))

# 5. FORECAST UNCERTAINTY METRIC 4:
#   Ambiguity Map (Lowest category probability)
# ============================================================
#   
#   This highlights pixels where the “least likely” category is still fairly likely.
# Interpretation:
#   
#   High ambiguity = even the least likely category has nontrivial probability → high uncertainty
# 
# Low ambiguity = one category dominates → high confidence


ambig <- minPr
names(ambig) <- "forecast_ambiguity"

plot(ambig,
     main = "Forecast Ambiguity (Min Probability)",
     col = hcl.colors(50, "Cividis"))

# 6. COMPOSITE UNCERTAINTY SUMMARY
# ==============================================

uncert_stack <- c(entropy, spread, conf, ambig)
names(uncert_stack)

plot(uncert_stack)

# 7. OPTIONAL: MAP MOST LIKELY CATEGORY
# ==============================================
# 1 = below
# 2 = near
# 3 = above
winner <- which.max(Pr_stack)
names(winner) <- "most_likely_category"

plot(winner, 
     main = "Most Likely Production Category",
     col=c("firebrick","gold","forestgreen"))

# 8. OPTIONAL: UNCERTAINTY MASKS
# ==============================================

high_uncert_mask <- entropy > 0.9   # adjust threshold 
plot(high_uncert_mask, main="High Uncertainty Areas")

high_conf_mask <- conf > 0.7
plot(high_conf_mask, main="High Confidence Areas")
#######

#### forecast uncertainty metrics function

forecast_probs <- function(P_forecast, th1_r, th2_r, b1_r, b2_r, season="Unknown") {
  # If P_forecast is just a number, convert to a raster with values
  if (is.numeric(P_forecast)) {
    template <- th1_r
    P_forecast <- setValues(template, rep(P_forecast, ncell(template)))
  }
  
  # === Clean rasters to remove inherited varnames/levels ===
  th1_r <- setValues(th1_r, values(th1_r))
  th2_r <- setValues(th2_r, values(th2_r))
  b1_r  <- setValues(b1_r,  values(b1_r))
  b2_r  <- setValues(b2_r,  values(b2_r))
  
  names(th1_r) <- paste0("th1_", season)
  names(th2_r) <- paste0("th2_", season)
  names(b1_r)  <- paste0("b1_", season)
  names(b2_r)  <- paste0("b2_", season)
  
  # === Intercepts ===
  a1 <- -b1_r * th1_r
  a2 <- -b2_r * th2_r
  names(a1) <- paste0("a1_", season)
  names(a2) <- paste0("a2_", season)
  
  # === Linear predictors ===
  eta1 <- a1 + b1_r * P_forecast
  eta2 <- a2 + b2_r * P_forecast
  
  # === Cumulative probabilities ===
  P_le_1 <- app(eta1, plogis)       # P(Y ≤ 1)
  P_le_2 <- app(eta2, plogis)       # P(Y ≤ 2)
  
  # Enforce ordinal structure: ensure P_le_1 ≤ P_le_2
  P_le_2 <- app(c(P_le_1, P_le_2), fun=function(x) max(x[1], x[2]))
  
  # === Category probabilities ===
  Pr_below <- P_le_1
  Pr_near  <- P_le_2 - P_le_1
  Pr_above <- 1 - P_le_2
  
  Pr_below <- clamp(Pr_below, 0, 1)
  Pr_near  <- clamp(Pr_near,  0, 1)
  Pr_above <- clamp(Pr_above, 0, 1)
  
  names(Pr_below) <- paste0("Pr_below_", season)
  names(Pr_near)  <- paste0("Pr_near_",  season)
  names(Pr_above) <- paste0("Pr_above_", season)
  
  # === Uncertainty metrics ===
  # Entropy
  entropy <- -(Pr_below*log(Pr_below) +
                 Pr_near *log(Pr_near)  +
                 Pr_above*log(Pr_above))
  
  names(entropy) <- paste0("entropy_", season)
  
  # Spread
  Pr_stack <- c(Pr_below, Pr_near, Pr_above)
  maxPr <- app(Pr_stack, max)
  minPr <- app(Pr_stack, min)
  spread <- maxPr - minPr
  names(spread) <- paste0("spread_", season)
  
  # Confidence = maxPr
  conf <- maxPr
  names(conf) <- paste0("confidence_", season)
  
  # Ambiguity = minPr
  ambig <- minPr
  names(ambig) <- paste0("ambiguity_", season)
  
  # Most likely category
  winner <- which.max(Pr_stack)
  names(winner) <- paste0("winner_", season)
  
  # === Return everything in a list ===
  return(list(
    Pr_below = Pr_below,
    Pr_near  = Pr_near,
    Pr_above = Pr_above,
    entropy  = entropy,
    spread   = spread,
    confidence = conf,
    ambiguity  = ambig,
    winner     = winner,
    sum_check  = Pr_below + Pr_near + Pr_above
  ))
}

# example usage 
out<- forecast_probs(
  P_forecast = 500,
  th1_r = th_Sep_bn_r,
  th2_r = th_Sep_na_r,
  b1_r  = b1_Sep_r,
  b2_r  = b2_Sep_r,
  season = "Sep"
)

out<- forecast_probs(
  P_forecast = 50,
  th1_r = th_Dec_bn_r,
  th2_r = th_Dec_na_r,
  b1_r  = b1_Dec_r,
  b2_r  = b2_Dec_r,
  season = "Dec"
)

plot(c(out$Pr_below, out$Pr_near, out$Pr_above), range=c(0,1))
plot(out$entropy,
     main="Forecast Entropy (Uncertainty of Tercile Prediction)\nHigher = Less Certain",
     col=hcl.colors(50, "Plasma"))

plot(out$confidence,
     main="Forecast Confidence (Peak Category Probability)\nHigher = More Certain",
     col=hcl.colors(50, "Viridis"),range=c(0,1))

plot(out$winner,
     main="Most Likely Vegetation Production Category\n1=Below, 2=Near, 3=Above",
     col=c("firebrick","gold","forestgreen"))

plot(
  out$spread,
  main = "Forecast Spread (max − min Probability)\nHigher = More Decisive Forecast",
  col = hcl.colors(50, "Viridis")
)

####


##### SINGLE YEAR/SEAS HINDCAST VERIFICATION CODE

hindcast_year <- function(year, season, 
                          th1_r, th2_r, b1_r, b2_r, 
                          Pcum_Dec, Pcum_Mar, Pcum_Jun, Pcum_Sep, 
                          veg_terc) {
  
  # Pick cumulative precip season
  Pobs <- switch(
    season,
    "Dec" = Pcum_Dec[[as.character(year)]],
    "Mar" = Pcum_Mar[[as.character(year)]],
    "Jun" = Pcum_Jun[[as.character(year)]],
    "Sep" = Pcum_Sep[[as.character(year)]],
    stop("Season must be one of: 'Dec','Mar','Jun','Sep'")
  )
  
  if (is.null(Pobs)) stop("No precipitation available for that year/season.")
  
  # Compute forecast probabilities using observed precip
  fc <- forecast_probs(
    P_forecast = Pobs,
    th1_r = th1_r,
    th2_r = th2_r,
    b1_r = b1_r,
    b2_r = b2_r,
    season = season
  )
  
  # Extract observed production tercile
  obs <- veg_terc[[as.character(year)]]
  
  # Winner (most likely forecast category)
  winner <- fc$winner
  
  # Correct/incorrect mask (numeric)
  correct <- ifel(winner == obs, 1, 0)
  names(correct) <- "forecast_correct"
  
  return(list(
    forecast = fc,
    observed = obs,
    winner   = winner,
    correct  = correct
  ))
}

# Example hindcast for 2015 Sep
hc <- hindcast_year(
  year = 2020,
  season = "Jun",
  th1_r = th_Jun_bn_r,
  th2_r = th_Jun_na_r,
  b1_r  = b1_Jun_r,
  b2_r  = b2_Jun_r,
  Pcum_Dec = Pcum_Dec,
  Pcum_Mar = Pcum_Mar,
  Pcum_Jun = Pcum_Jun,
  Pcum_Sep = Pcum_Sep,
  veg_terc = veg_terc
)

# Plot correct/incorrect map
plot(hc$correct,
     main="Hindcast Accuracy (2020, Sep)\n1 = Correct, 0 = Incorrect",
     col=c("red","green"))

plot(c(
  hc$forecast$Pr_below, 
  hc$forecast$Pr_near, 
  hc$forecast$Pr_above
), range=c(0,1))

plot(hc$winner, 
     main="Most Likely Category (Forecast)")

plot(hc$observed, 
     main="Observed Category (Vegetation Tercile)")

plot(hc$correct, 
     main="Forecast Accuracy (Correct / Incorrect)",
     col=c("red", "green"))
######

##### ALL YEARS HINDCAST VERIFICATION CODE -----

hindcast_all_years <- function(
    season,
    th1_r, th2_r, b1_r, b2_r,
    Pcum_Dec, Pcum_Mar, Pcum_Jun, Pcum_Sep,
    veg_terc
) {
  
  years <- intersect(
    as.integer(names(veg_terc)),
    as.integer(names(Pcum_Sep))
  )
  
  results <- list()
  acc_stack <- list()  # will become a SpatRaster later
  
  for (yr in years) {
    message("Hindcasting year ", yr, " (", season, ")...")
    
    hc <- hindcast_year(
      year     = yr,
      season   = season,
      th1_r    = th1_r,
      th2_r    = th2_r,
      b1_r     = b1_r,
      b2_r     = b2_r,
      Pcum_Dec = Pcum_Dec,
      Pcum_Mar = Pcum_Mar,
      Pcum_Jun = Pcum_Jun,
      Pcum_Sep = Pcum_Sep,
      veg_terc = veg_terc
    )
    
    # store
    results[[as.character(yr)]] <- hc
    acc_stack[[as.character(yr)]] <- hc$correct
  }
  
  # Convert to SpatRaster
  acc_stack <- rast(acc_stack)
  
  return(list(
    results = results,
    accuracy_stack = acc_stack,
    years = years
  ))
}

# Example: hindcast all years for Sep
hind_sep <- hindcast_all_years(
  season = "Sep",
  th1_r = th_Sep_bn_r,
  th2_r = th_Sep_na_r,
  b1_r  = b1_Sep_r,
  b2_r  = b2_Sep_r,
  Pcum_Dec = Pcum_Dec,
  Pcum_Mar = Pcum_Mar,
  Pcum_Jun = Pcum_Jun,
  Pcum_Sep = Pcum_Sep,
  veg_terc = veg_terc
)

acc_mean <- app(hind_sep$accuracy_stack, mean, na.rm=TRUE)
names(acc_mean) <- "hindcast_accuracy"

plot(acc_mean,
     main="Hindcast Accuracy (All Years, Sep)\nMean Correct Forecast Probability",
     col=hcl.colors(50, "YlGn"))

acc_count <- app(hind_sep$accuracy_stack, sum, na.rm=TRUE)
names(acc_count) <- "correct_years"

plot(acc_count,
     main="Number of Correct Forecast Years (Sep)",
     col=hcl.colors(50, "Blues"))

hc2020 <- hind_sep$results[["2020"]]

plot(c(hc2020$winner, hc2020$observed),
     main = "2020 Hindcast: Forecast vs Observed")
plot(hc2020$correct, col=c("red","green"),
     main="2020 Forecast Correctness")

acc_class <- classify(
  acc_mean,
  rcl = matrix(c(
    0.00, 0.33, 1,
    0.33, 0.66, 2,
    0.66, 1.00, 3
  ), ncol=3, byrow=TRUE)
)

plot(
  acc_class,
  col = c("red","gold","darkgreen"),
  main = "Hindcast Skill Classes\n1=Low, 2=Medium, 3=High"
)

# all years plot
library(tidyterra)
library(ggplot2)

acc <- hind_sep$accuracy_stack

df <- as.data.frame(acc, xy = TRUE, na.rm = FALSE)
df_long <- tidyr::pivot_longer(df, -c(x,y), names_to="year", values_to="correct")

df_long$correct_f <- factor(df_long$correct, levels=c(0,1),
                            labels=c("Incorrect","Correct"))

ggplot(df_long) +
  geom_raster(aes(x, y, fill = correct_f)) +
  geom_sf(data = shp, fill = NA, color = "black", size = 0.3) +   # <--- polygon outline
  facet_wrap(~year, ncol = 9) +
  scale_fill_manual(values=c("Incorrect"="red", "Correct"="green")) +
  theme_minimal() +
  labs(title="Hindcast Accuracy for All Years",
       fill="Forecast") +
  coord_sf()


##### TREND and BREAKPOINT ANALYSIS -----

# plot veg terciles
library(tidyterra)
library(ggplot2)
library(dplyr)
library(tidyr)

veg_df <- as.data.frame(veg_terc, xy = TRUE, na.rm = FALSE)
veg_long <- pivot_longer(
  veg_df,
  cols = -c(x, y),
  names_to = "year",
  values_to = "tercile"
)

veg_long$tercile_f <- factor(
  veg_long$tercile,
  levels = c(1, 2, 3),
  labels = c("Below Normal", "Near Normal", "Above Normal")
)

ggplot(veg_long) +
  geom_raster(aes(x, y, fill = tercile_f)) +
  geom_sf(data = shp, fill = NA, color = "black", size = 0.3) +   # <--- polygon outline
  facet_wrap(~year, ncol = 9) +
  scale_fill_manual(
    values = c(
      "Below Normal" = "firebrick",
      "Near Normal"  = "gold",
      "Above Normal" = "forestgreen"
    )
  ) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    title = "Vegetation Production Terciles (All Years)",
    fill = "Tercile Category"
  ) +
  coord_sf()

# pixel wise trends
library(terra)
library(trend)

mk_trend <- app(veg_prod, fun = function(x) {
  if(all(is.na(x))) return(NA)
  mk <- mk.test(x)
  mk$estimates[1]  # Sen slope
})

pct_by_year <- apply(values(veg_terc), 2, function(v) {
  c(
    below = mean(v == 1, na.rm=TRUE),
    near  = mean(v == 2, na.rm=TRUE),
    above = mean(v == 3, na.rm=TRUE)
  )
})
matplot(t(pct_by_year), type="l", lwd=2)

# Extract correct years from veg_terc
years <- as.numeric(names(veg_terc))

# Matrix of tercile fractions (year × 3)
pct_mat <- t(pct_by_year)

# Colors
cols <- c(
  "firebrick",   # below
  "gold",        # near
  "forestgreen"  # above
)

# Plot
matplot(
  years, pct_mat,
  type = "l", lwd = 2, col = cols,
  xlab = "Year",
  ylab = "Fraction of Pixels",
  main = "Vegetation Production Tercile Fractions Over Time",
  ylim = c(0,1)
)

legend(
  "topleft",
  legend = c("Below Normal", "Near Normal", "Above Normal"),
  col = cols,
  lwd = 2,
  bty = "n"
)



#####

###############################################################
# BREAKPOINT DIAGNOSTICS FOR VEGETATION PRODUCTION TIME SERIES
# -------------------------------------------------------------
# This script:
#   1. Extracts annual mean vegetation production from raster stack
#   2. Runs structural breakpoint detection (strucchange)
#   3. Plots BIC/RSS diagnostics
#   4. Identifies optimal breakpoint year(s)
#   5. Plots observed vs segmented (stepwise) means with breakpoints
###############################################################

library(terra)
library(strucchange)

# -------------------------------------------------------------
# 1. Extract spatial mean vegetation production per year
# -------------------------------------------------------------
veg <- veg_prod

years <- as.numeric(names(veg))    # numeric year labels
veg_mean <- global(veg, fun="mean", na.rm=TRUE)[,1]

ts_df <- data.frame(
  year = years,
  mean_prod = veg_mean
)

# -------------------------------------------------------------
# 2. Fit structural break model (intercept-only)
# -------------------------------------------------------------
# This tests whether the mean level of vegetation production 
# changes abruptly at any point in time.
bp_model <- breakpoints(mean_prod ~ 1, data = ts_df)

print(bp_model)

# -------------------------------------------------------------
# 3. Plot BIC and RSS diagnostics
# -------------------------------------------------------------
# This shows how model fit improves as breakpoints are added.
# Look for a sharp decrease at k=1 → strong evidence of a jump.
plot(bp_model,
     main = "Breakpoint Diagnostics (RSS and BIC)\nIdentify Number of Structural Breaks")

# -------------------------------------------------------------
# 4. Extract estimated breakpoints (index positions)
# -------------------------------------------------------------
best_break_index <- bp_model$breakpoints
break_years <- years[best_break_index]

cat("Estimated Break Year(s):", break_years, "\n")

# -------------------------------------------------------------
# 5. Plot time series with fitted segmented means + break lines
# -------------------------------------------------------------
plot(ts_df$year, ts_df$mean_prod,
     type="b", pch=16, col="black",
     xlab="Year", ylab="Mean Vegetation Production",
     main="Vegetation Production Time Series with Detected Breakpoints")

# Add breakpoint vertical lines
abline(v = break_years, col="red", lwd=3)

# Add segmented mean estimate
segmented_means <- fitted(bp_model)
lines(ts_df$year, segmented_means, col="blue", lwd=3)

legend("topleft",
       legend=c("Observed Mean", "Segmented Mean", "Breakpoint"),
       col=c("black", "blue", "red"),
       lwd=c(1,3,3), pch=c(16, NA, NA))

# -------------------------------------------------------------
# 6. Optional: Confidence intervals for breakpoint location
# -------------------------------------------------------------
print(confint(bp_model))
#####

###### Precip threshold maps ------

library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)

# --- Put all threshold rasters in a named list ---
threshold_rasters <- list(
  Dec_BN = th_Dec_bn_r,
  Dec_NA = th_Dec_na_r,
  Mar_BN = th_Mar_bn_r,
  Mar_NA = th_Mar_na_r,
  Jun_BN = th_Jun_bn_r,
  Jun_NA = th_Jun_na_r,
  Sep_BN = th_Sep_bn_r,
  Sep_NA = th_Sep_na_r
)

# Convert rasters → data frames and combine
threshold_df <- lapply(names(threshold_rasters), function(nm) {
  df <- as.data.frame(threshold_rasters[[nm]], xy = TRUE, na.rm = TRUE)
  names(df)[3] <- "threshold"
  
  # Split name into Season + Type
  parts <- strsplit(nm, "_")[[1]]
  df$Season <- parts[1]
  df$Type   <- ifelse(parts[2] == "BN", "Below→Near", "Near→Above")
  
  df
}) %>% bind_rows()


library(rasterVis)

bwplot(veg_prod, main="KNF RPMS Vegetation Production Time Series", do.out=FALSE,
       ylim = c(0, 1200) )

library(terra)
library(dplyr)
library(tidyr)
library(ggplot2)

# Convert raster to long table: year, value
df <- as.data.frame(veg_prod, xy = FALSE) %>%
  pivot_longer(cols = everything(),
               names_to = "year",
               values_to = "production") %>%
  mutate(year = as.numeric(year))

# Plot with no whiskers (geom_boxplot with coef = 0 or boxplotstats override)
ggplot(df, aes(x = factor(year), y = production)) +
  geom_boxplot(outlier.shape = NA, coef = 0) +    
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  labs(title = "Vegetation Production Time Series",
       x = "Year",
       y = "Production Value") +
  coord_cartesian(ylim = c(250, 1000))   # <--- This zooms without dropping data



