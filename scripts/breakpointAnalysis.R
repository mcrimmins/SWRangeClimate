
library(terra)
library(strucchange)


# set AOI
shp <- sf::st_read(dsn = "./data/shapes/AdministrativeForest.gdb")
shp<- sf::st_transform(shp, crs = sf::st_crs("+proj=longlat +datum=WGS84"))

# subset to selected NF
shp<-subset(shp, shp$FORESTNAME=="Kaibab National Forest")
#shp<-subset(shp, shp$FORESTNAME=="Tonto National Forest")
#shp<-subset(shp, shp$FORESTNAME=="Coronado National Forest")

# Load SpatRaster data -----
#veg_prod <- rast("./data/processed/800m_KNF_RPMS_1984_2024.tif")  # Annual vegetation production
#veg_prod <- rast("./data/processed/800m_KNF_RPMS_1984_2024_EVTmask.tif")  # Annual vegetation production with EVT mask
veg_prod<- rast("./data/processed/KNF_RPMS_1984_2024_EVT_treeMask.tif")

#veg_prod<-rast("./data/processed/TontoNF_RPMS_1984_2024_no_mask.tif")


###############################################################
# A. OPTIONAL SENSOR-SHIFT CORRECTION
###############################################################

veg_input <- veg_prod  # keep a clean copy of the original
rm(veg_prod)

# ============================================================
# FLEXIBLE SENSOR-SHIFT CORRECTION BLOCK
# Supports:
#   (1) auto_correction = TRUE  → match early/late via mean differences
#   (2) manual_factor   = 0.87  → reduce late period by 13%
# ============================================================

apply_correction <- FALSE          # master switch
auto_correction  <- FALSE          # TRUE = use automatic shift; FALSE = use manual factor
manual_factor    <- 0.87         # e.g., 0.87 for 13% test; 1.00 = no manual correction

if (apply_correction) {
  
  message("Applying vegetation sensor-shift correction...")
  
  years <- as.numeric(names(veg_input))
  early_period <- years <= break_year
  late_period  <- years >  break_year
  
  veg_corrected <- veg_input
  
  if (auto_correction) {
    
    message("Using AUTOMATIC correction based on early/late mean differences...")
    
    early_mean <- app(veg_input[[early_period]], mean, na.rm = TRUE)
    late_mean  <- app(veg_input[[late_period]],  mean, na.rm = TRUE)
    
    shift_amount <- early_mean - late_mean
    
    # Apply pixel-wise shift to each late-period layer
    for (i in which(late_period)) {
      veg_corrected[[i]] <- veg_input[[i]] + shift_amount
    }
    
    message("Automatic correction applied.")
  }
  
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
veg_mean <- global(veg_used, fun="mean", na.rm=TRUE)[,1]

ts_df <- data.frame(
  year = years,
  mean_prod = veg_mean
)

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

