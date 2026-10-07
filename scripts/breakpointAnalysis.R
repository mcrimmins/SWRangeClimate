
library(terra)
library(strucchange)


# set AOI
shp <- sf::st_read(dsn = "./data/shapes/AdministrativeForest.gdb")
shp<- sf::st_transform(shp, crs = sf::st_crs("+proj=longlat +datum=WGS84"))

# subset to selected NF
#shp<-subset(shp, shp$FORESTNAME=="Kaibab National Forest")
#shp<-subset(shp, shp$FORESTNAME=="Tonto National Forest")
shp<-subset(shp, shp$FORESTNAME=="Coronado National Forest")

# Load SpatRaster data -----
#veg_prod <- rast("./data/processed/800m_KNF_RPMS_1984_2024.tif")  # Annual vegetation production
#veg_prod <- rast("./data/processed/800m_KNF_RPMS_1984_2024_EVTmask.tif")  # Annual vegetation production with EVT mask
#veg_prod<- rast("./data/processed/KNF_RPMS_1984_2024_EVT_treeMask.tif")

veg_prod<-rast("./data/processed/CoronadoNF_RPMS_1984_2024_no_mask.tif")


###############################################################
# A. OPTIONAL SENSOR-SHIFT CORRECTION
###############################################################

veg_input <- veg_prod  # keep a clean copy of the original
rm(veg_prod)

# ============================================================
# Detect breakpoint (required!)
# ============================================================

years_vec <- as.numeric(names(veg_input))
veg_mean_ts <- global(veg_input, fun = "mean", na.rm = TRUE)[,1]

bp_model <- breakpoints(veg_mean_ts ~ 1)
bp_index <- bp_model$breakpoints[1]
break_year <- years_vec[bp_index]

message("Detected breakpoint around year: ", break_year)

# ============================================================
# FLEXIBLE SENSOR-SHIFT CORRECTION BLOCK
# ============================================================

apply_correction <- TRUE          # master switch
auto_correction  <- TRUE          # TRUE = automatic shift; FALSE = manual factor
manual_factor    <- 0.87          # e.g., reduce late period by 13%

if (apply_correction) {
  
  message("Applying vegetation sensor-shift correction...")
  
  years <- as.numeric(names(veg_input))
  early_period <- years <= break_year
  late_period  <- years >  break_year
  
  veg_corrected <- veg_input
  
  # -------------------------------
  # Mode 1 — Automatic correction
  # -------------------------------
  if (auto_correction) {
    
    message("Using AUTOMATIC correction based on early/late mean differences...")
    
    early_mean <- app(veg_input[[early_period]], mean, na.rm = TRUE)
    late_mean  <- app(veg_input[[late_period]],  mean, na.rm = TRUE)
    
    shift_amount <- early_mean - late_mean
    
    for (i in which(late_period)) {
      veg_corrected[[i]] <- veg_input[[i]] + shift_amount
    }
    
    
    early_mean_global <- global(early_mean, fun = "mean", na.rm = TRUE)[1,1]
    late_mean_global  <- global(late_mean,  fun = "mean", na.rm = TRUE)[1,1]
    
    early_mean_global
    late_mean_global
    
    
    message("Automatic correction applied.")
  }
  
  # -------------------------------
  # Mode 2 — Manual multiplicative test
  # -------------------------------
  if (!auto_correction && manual_factor != 1.00) {
    
    message("Using MANUAL multiplicative correction. Factor = ", manual_factor)
    
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
# B. BREAKPOINT DIAGNOSTICS (using corrected or raw data)
###############################################################

# Extract annual mean values
years <- as.numeric(names(veg_used))
#veg_mean <- global(veg_used, fun="mean", na.rm=TRUE)[,1]

ts_df <- data.frame(
  year = years,
  mean_prod = veg_mean_ts
)

early_mean_global <- mean(veg_mean_ts[years <= break_year])
late_mean_global  <- mean(veg_mean_ts[years > break_year])

early_mean_global
late_mean_global


# -------------------------------------------------------------
# 1. Fit breakpoints model
# -------------------------------------------------------------
bp_model <- breakpoints(mean_prod ~ 1, data = ts_df)
print(bp_model)

# -------------------------------------------------------------
# 2. Plot BIC and RSS diagnostics
# -------------------------------------------------------------
plot(bp_model,
     main = "Breakpoint Diagnostics (RSS & BIC)\nRaw or Corrected Vegetation Time Series")

# -------------------------------------------------------------
# 3. Extract breakpoint years
# -------------------------------------------------------------
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

veg_input<-veg_prod
names(veg_input)

plot(veg_input[[38]], range=c(0,5000),main=paste0(names(veg_input[[38]])), col = rev(hcl.colors(50, "Inferno")))
plot(shp, add=TRUE, color=NA)

plot(veg_input[[26:41]], range=c(0,5000))

