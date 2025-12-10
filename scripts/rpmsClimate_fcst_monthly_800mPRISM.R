# RPMS climate analysis for KNF
# using 800m PRISM data and EVT
# MAC 11/5/25

library(terra)
library(ggplot2)
library(lubridate)
library(RColorBrewer)
library(dplyr)

terraOptions(progress=1)

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


##### 0. Harmonize rasters and time indices

# project precip to veg_prod grid/CRS
precip <- project(precip, veg_prod, method = "bilinear")

# make sure extents, resolution, alignment all match
precip <- crop(precip, veg_prod)
#precip <- align(precip, veg_prod)  # only if needed

# precip layer names are "YYYYMM"
precip_dates <- as.Date(paste0(names(precip), "01"), format = "%Y%m%d")
yr  <- year(precip_dates)
mon <- month(precip_dates)

# water-year tied to October peak:
# Oct–Dec belong to WY = calendar year + 1
# Jan–Sep belong to WY = calendar year
wy_year  <- ifelse(mon >= 10, yr + 1, yr)
wy_month <- mon

# veg production years available (drop 2012)
veg_years <- as.integer(names(veg_prod))
veg_years_use <- setdiff(veg_years, 2012L)

# keep only precip months whose WY is in veg_years_use
keep <- wy_year %in% veg_years_use
precip   <- precip[[keep]]
wy_year  <- wy_year[keep]
wy_month <- wy_month[keep]

#####


##### 1. Define water-year windows and build precip summaries

make_window_sum <- function(precip, wy_year, wy_month, months, years_out) {
  sel <- wy_month %in% months & wy_year %in% years_out
  g   <- wy_year[sel]  # group index = water year
  
  r <- tapp(precip[[sel]], index = g, fun = sum, na.rm = TRUE)
  names(r) <- sort(unique(g))
  r
}

years_out <- veg_years_use

precip_cool   <- make_window_sum(precip, wy_year, wy_month,
                                 months = c(10, 11, 12, 1, 2, 3),
                                 years_out = years_out)

precip_spring <- make_window_sum(precip, wy_year, wy_month,
                                 months = 4:6,
                                 years_out = years_out)

precip_monsoon <- make_window_sum(precip, wy_year, wy_month,
                                  months = 7:9,
                                  years_out = years_out)

precip_wy     <- make_window_sum(precip, wy_year, wy_month,
                                 months = 1:12,  # all WY months (since wy_year already assigned)
                                 years_out = years_out)

order_idx <- match(as.character(veg_years_use), names(precip_wy))

precip_cool    <- precip_cool[[order_idx]]
precip_spring  <- precip_spring[[order_idx]]
precip_monsoon <- precip_monsoon[[order_idx]]
precip_wy      <- precip_wy[[order_idx]]

names(precip_cool)    <- veg_years_use
names(precip_spring)  <- veg_years_use
names(precip_monsoon) <- veg_years_use
names(precip_wy)      <- veg_years_use

#####

##### 2. Standardize local precip as percentiles (per pixel)

ecdf_rank <- function(x) {
  if (all(is.na(x))) return(x)
  r <- rank(x, ties.method = "average", na.last = "keep")
  r / (sum(!is.na(x)) + 1)   # in (0,1)
}

precip_cool_pct    <- app(precip_cool,    ecdf_rank)
precip_spring_pct  <- app(precip_spring,  ecdf_rank)
precip_monsoon_pct <- app(precip_monsoon, ecdf_rank)
precip_wy_pct      <- app(precip_wy,      ecdf_rank)

#####

##### 3. Define production categories (binary or terciles)

veg <- veg_prod[[names(veg_prod) %in% veg_years_use]]
veg_years_use <- as.integer(names(veg))  # refresh, just in case

# check alignment: should match
stopifnot(all(veg_years_use == veg_years_use))

# 3a. Binary above / below threshold (e.g., median)
make_binary_cat <- function(x, q = 0.5) {
  if (all(is.na(x))) return(x)
  thr <- quantile(x, probs = q, na.rm = TRUE, type = 8)
  as.integer(x > thr)  # 1 = above, 0 = at or below
}

veg_bin <- app(veg, make_binary_cat)
names(veg_bin) <- veg_years_use

# 3b. Terciles (1 = below, 2 = near, 3 = above)
make_tercile_cat <- function(x, probs = c(1/3, 2/3)) {
  if (all(is.na(x))) return(x)
  qs <- quantile(x, probs = probs, na.rm = TRUE, type = 8)
  cut(x,
      breaks = c(-Inf, qs[1], qs[2], Inf),
      labels = FALSE)  # 1,2,3
}

veg_terc <- app(veg, make_tercile_cat)
names(veg_terc) <- veg_years_use

# 4. Add vegetation type
veg_type <- project(veg_type, veg, method = "near")  # categorical, so nearest

# optional: restrict to dominant classes, set NA in others, etc.
######

###### 5. Build a modeling dataset (pixel × year)
# align all rasters to the same set of years
years_str <- as.character(veg_years_use)

veg_bin         <- veg_bin[[years_str]]
veg_terc        <- veg_terc[[years_str]]
precip_cool_pct <- precip_cool_pct[[years_str]]
precip_monsoon_pct <- precip_monsoon_pct[[years_str]]
precip_wy_pct   <- precip_wy_pct[[years_str]]

# rename layers with a prefix for clarity
names(veg_bin)            <- paste0("yb_", years_str)
names(veg_terc)           <- paste0("yt_", years_str)
names(precip_cool_pct)    <- paste0("pcool_", years_str)
names(precip_monsoon_pct) <- paste0("pmon_", years_str)
names(precip_wy_pct)      <- paste0("pwy_", years_str)

build_df_for_year <- function(y) {
  ys <- as.character(y)
  
  r_stack <- c(
    veg_bin[[paste0("yb_", ys)]],
    veg_terc[[paste0("yt_", ys)]],
    precip_cool_pct[[paste0("pcool_", ys)]],
    precip_monsoon_pct[[paste0("pmon_", ys)]],
    precip_wy_pct[[paste0("pwy_", ys)]],
    veg_type
  )
  
  names(r_stack) <- c("prod_bin", "prod_terc",
                      "pcool_pct", "pmon_pct", "pwy_pct",
                      "veg_type")
  
  df <- as.data.frame(r_stack, xy = FALSE, cells = TRUE, na.rm = TRUE)
  df$year <- y
  df
}

df_list <- lapply(veg_years_use, build_df_for_year)
dat_all <- bind_rows(df_list)

# turn into factors where appropriate
dat_all$prod_bin  <- factor(dat_all$prod_bin, levels = c(0, 1), labels = c("below_or_equal", "above"))
dat_all$prod_terc <- factor(dat_all$prod_terc, levels = c(1, 2, 3),
                            labels = c("below", "near", "above"))

dat_all$veg_type  <- as.factor(dat_all$veg_type)
######

###### 6. Fit probability models
# 6a. Binary logistic (above vs below-median production)
m_bin <- glm(prod_bin ~ pcool_pct + pmon_pct + pwy_pct + veg_type,
             data = dat_all,
             family = binomial)

summary(m_bin)

# interpretive plots

library(ggplot2)

# Generate a grid of predictor values
newdat <- data.frame(
  pcool_pct = seq(0, 1, length = 100),
  pmon_pct  = mean(dat_all$pmon_pct,  na.rm = TRUE),
  pwy_pct   = mean(dat_all$pwy_pct,   na.rm = TRUE),
  veg_type  = names(sort(table(dat_all$veg_type), decreasing = TRUE))[1]
)

newdat$prob <- predict(m_bin, newdata = newdat, type = "response")

ggplot(newdat, aes(x = pcool_pct, y = prob)) +
  geom_line(size = 1.2) +
  labs(
    x = "Cool-season precipitation percentile",
    y = "Probability of above-median annual production",
    title = "Effect of Cool-Season Precipitation"
  ) +
  ylim(0, 1)

# for monsoon season
newdat2 <- newdat
newdat2$pmon_pct <- seq(0, 1, length = 100)
newdat2$pcool_pct <- mean(dat_all$pcool_pct, na.rm = TRUE)

newdat2$prob <- predict(m_bin, newdata = newdat2, type = "response")

ggplot(newdat2, aes(x = pmon_pct, y = prob)) +
  geom_line(size = 1.2, color = "firebrick") +
  labs(
    x = "Monsoon precipitation percentile",
    y = "Probability of above-median annual production",
    title = "Effect of Monsoon Precipitation"
  ) +
  ylim(0, 1)

## combine effect
# grid of predictors
grid <- expand.grid(
  pcool_pct = seq(0,1,length=50),
  pmon_pct  = seq(0,1,length=50)
)

grid$pwy_pct <- mean(dat_all$pwy_pct, na.rm = TRUE)
grid$veg_type <- names(sort(table(dat_all$veg_type), decreasing = TRUE))[1]

grid$prob <- predict(m_bin, newdata = grid, type = "response")

ggplot(grid, aes(pcool_pct, pmon_pct, fill = prob)) +
  geom_tile() +
  scale_fill_viridis_c(option = "C", limits = c(0,1)) +
  labs(
    x = "Cool-season precipitation percentile",
    y = "Monsoon precipitation percentile",
    fill = "P(above)",
    title = "Combined Effect of Cool and Monsoon Precip Percentiles"
  )

# coefficient sizes
library(broom)

coef_df <- broom::tidy(m_bin) %>% 
  filter(term %in% c("pcool_pct", "pmon_pct", "pwy_pct"))

coef_df$OR <- exp(coef_df$estimate)

ggplot(coef_df, aes(x = term, y = OR)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = exp(estimate - 1.96*std.error),
                    ymax = exp(estimate + 1.96*std.error)),
                width = 0.1, size = 1) +
  labs(
    x = "Predictor",
    y = "Odds Ratio (±95% CI)",
    title = "Seasonal Precipitation Effects on Biomass Probability"
  ) +
  theme_minimal()




## does not work
# vegetation type
# models_by_type <- dat_all %>%
#   group_split(veg_type) %>%
#   purrr::map(~ glm(prod_bin ~ pcool_pct + pmon_pct + pwy_pct,
#                    data = .x, family = binomial))
# names(models_by_type) <- levels(dat_all$veg_type)

# 6b. Terciles (below / near / above) via multinomial or ordinal models
# install.packages("ordinal")
library(ordinal)

m_terc <- clm(prod_terc ~ pcool_pct + pmon_pct + pwy_pct + veg_type,
              data = dat_all,
              link = "logit")
summary(m_terc)
######

###### 7. Map predicted probabilities (per pixel, for each year)
predict_year_raster <- function(y, model, type = c("binary", "tercile")) {
  type <- match.arg(type)
  ys <- as.character(y)
  
  # stack predictors for that year
  pred_stack <- c(
    precip_cool_pct[[paste0("pcool_", ys)]],
    precip_monsoon_pct[[paste0("pmon_", ys)]],
    precip_wy_pct[[paste0("pwy_", ys)]],
    veg_type
  )
  
  names(pred_stack) <- c("pcool_pct", "pmon_pct", "pwy_pct", "veg_type")
  
  df_pred <- as.data.frame(pred_stack, cells = TRUE, na.rm = FALSE)
  
  if (type == "binary") {
    p <- predict(model, newdata = df_pred, type = "response")  # P(above)
    out <- rast(veg[[1]])
    values(out) <- p
    names(out) <- paste0("p_above_", ys)
  } else {
    # ordinal: need matrix of probabilities for each class
    pmat <- predict(model, newdata = df_pred, type = "prob")
    out <- rast(veg[[1]], nlyrs = 3)
    values(out) <- as.matrix(pmat)
    names(out) <- paste0(c("p_below_", "p_near_", "p_above_"), ys)
  }
  
  out
}

p_map_2019 <- predict_year_raster(2019, m_bin, type = "binary")
plot(p_map_2019, main = "P(above-median production) for 2019")
#####

##### 8. Pixel-level “confidence” metrics
# ===============================================================
#  Packages
# ===============================================================
library(terra)
library(dplyr)
library(ggplot2)
library(broom)

# ===============================================================
#  User Inputs (already built earlier in workflow)
# ===============================================================
# Required objects already in your environment:
#   veg         — SpatRaster (vegetation production grid, template)
#   dat_all     — data frame with columns:
#                 cell, year, prod_bin, pcool_pct, pmon_pct,
#                 pwy_pct, veg_type
#
#   prod_bin    — factor with "above" and "below"
#
# Assumes dat_all was constructed correctly using the per-year builder.

# ===============================================================
#  Add row/column indices to dat_all
# ===============================================================
base <- veg[[1]]

coords <- terra::rowColFromCell(base, dat_all$cell)
dat_all$row <- coords[,1]
dat_all$col <- coords[,2]

# Ensure observed binary is numeric 0/1
dat_all$y_obs <- as.numeric(dat_all$prod_bin == "above")

# ===============================================================
#  LOYO Cross Validation (produces cv_df)
# ===============================================================
years_use <- sort(unique(dat_all$year))

loyo_results <- lapply(years_use, function(ytest) {
  
  train <- dat_all %>% filter(year != ytest)
  test  <- dat_all %>% filter(year == ytest)
  
  # logistic regression model
  m <- glm(prod_bin ~ pcool_pct + pmon_pct + pwy_pct + veg_type,
           data = train,
           family = binomial)
  
  # Predict probability for withheld year
  test$pred <- predict(m, newdata = test, type = "response")
  test$year_test <- ytest
  
  test
})

cv_df <- bind_rows(loyo_results)

# ===============================================================
#  Skill Raster Builder Function
# ===============================================================
make_skill_raster <- function(cv_df, veg, stat = c("cor", "hitrate", "both")) {
  stat <- match.arg(stat)
  
  base <- veg[[1]]
  
  # Build lookup table for all pixels in the raster
  base_cells <- 1:ncell(base)
  base_rc <- terra::rowColFromCell(base, base_cells)
  base_map <- data.frame(
    cell = base_cells,
    row  = base_rc[,1],
    col  = base_rc[,2]
  )
  
  # Ensure row/col exist
  if (!("row" %in% names(cv_df))) {
    rc <- terra::rowColFromCell(base, cv_df$cell)
    cv_df$row <- rc[,1]
    cv_df$col <- rc[,2]
  }
  
  # Observed binary outcome
  cv_df$y_obs <- as.numeric(cv_df$prod_bin == "above")
  
  # Compute skill per (row,col)
  skill_df <- cv_df %>%
    group_by(row, col) %>%
    summarize(
      cor_py = if (stat %in% c("cor","both") &&
                   length(unique(y_obs)) > 1)
        cor(pred, y_obs, use = "complete.obs") else NA_real_,
      
      hitrate = if (stat %in% c("hitrate","both"))
        mean((pred > 0.5) == (y_obs == 1), na.rm = TRUE) else NA_real_,
      .groups = "drop"
    )
  
  # Merge skill values to full pixel grid
  merged <- left_join(base_map, skill_df, by = c("row","col"))
  
  out_list <- list()
  
  # Correlation raster
  if (stat %in% c("cor","both")) {
    r_cor <- base
    terra::values(r_cor) <- merged$cor_py
    names(r_cor) <- "skill_cor"
    out_list$cor <- r_cor
  }
  
  # Hit rate raster
  if (stat %in% c("hitrate","both")) {
    r_hit <- base
    terra::values(r_hit) <- merged$hitrate
    names(r_hit) <- "skill_hitrate"
    out_list$hitrate <- r_hit
  }
  
  if (stat == "cor")     return(out_list$cor)
  if (stat == "hitrate") return(out_list$hitrate)
  return(out_list)
}

# ===============================================================
#  Build Skill Rasters
# ===============================================================
skill_cor_r <- make_skill_raster(cv_df, veg, stat = "cor")
skill_hit_r <- make_skill_raster(cv_df, veg, stat = "hitrate")

# ===============================================================
#  Plot Skill Maps
# ===============================================================
plot(skill_cor_r, main = "Pixel-Level Skill (Correlation)")
plot(skill_hit_r, main = "Pixel-Level Skill (Hit Rate)")

# ===============================================================
#  Save (optional)
# ===============================================================
# terra::writeRaster(skill_cor_r, "skill_correlation.tif", overwrite=TRUE)
# terra::writeRaster(skill_hit_r, "skill_hitrate.tif", overwrite=TRUE)

# 1. Join skill values to vegetation type (per-pixel overlay)
# Extract skill and vegetation type values
skill_cor_df <- as.data.frame(skill_cor_r, cells = TRUE, xy = FALSE)
#colnames(skill_cor_df) <- c("skill_cor", "cell")
veg_df <- as.data.frame(veg_type, cells = TRUE, xy = FALSE)
skill_df <- left_join(skill_cor_df, veg_df, by = "cell")

skill_hit_df <- as.data.frame(skill_hit_r, cells = TRUE, xy = FALSE)
#colnames(skill_cor_df) <- c("skill_cor", "cell")
#veg_df <- as.data.frame(veg_type, cells = TRUE, xy = FALSE)
skill_df <- left_join(skill_hit_df, veg_df, by = "cell")

skill_by_type <- skill_df %>%
  group_by(EVT_NAME) %>%
  summarize(
    n_pixels = sum(!is.na(skill_hitrate)),
    mean_hit = mean(skill_hitrate, na.rm = TRUE),
    sd_hit   = sd(skill_hitrate, na.rm = TRUE)
  ) %>%
  arrange(desc(mean_hit))

