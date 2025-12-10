# RPMS climate analysis for KNF
# using 800m PRISM data and EVT
# MAC 11/5/25

library(terra)
library(ggplot2)
library(lubridate)
library(RColorBrewer)

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

# global relationships
#precip_annual <- tapp(precip, index = year(dates(precip)), fun = sum)
#cor_vals <- global(c(precip_annual, veg_prod), fun = cor, na.rm = TRUE)

##### analysis code
library(terra)
library(dplyr)
library(tidyr)
library(purrr)

# 0a) Align CRS and extents (reproject precip to WGS84 to match veg_prod if needed)
# if (!identical(crs(precip), crs(veg_prod))) {
#   precip <- project(precip, crs(veg_prod))
# }

# 0b) Crop/mask to common extent
ext_common <- intersect(ext(veg_prod), ext(precip))
veg_prod  <- crop(veg_prod,  ext_common)
precip    <- crop(precip,    ext_common)
veg_type  <- crop(veg_type,  ext_common)

# apply veg_type as mask to precip and to veg_prod
precip <- mask(precip, veg_type)
veg_prod <- mask(veg_prod, veg_type)


# 0c) Attach time to precip from names like "YYYYMM"
get_dates_from_names <- function(nm) {
  as.Date(paste0(substr(nm,1,4), "-", substr(nm,5,6), "-01"))
}
time(precip) <- get_dates_from_names(names(precip))

# 0d) Helper indices for seasons; water year = Oct–Sep
p_dates <- time(precip)
p_year  <- as.integer(format(p_dates, "%Y"))
p_month <- as.integer(format(p_dates, "%m"))

is_monsoon <- p_month %in% 6:9            # Jun–Sep (adjust if you prefer Jul–Sep)
is_cool    <- p_month %in% c(11,12,1,2,3,4)
wy_year    <- ifelse(p_month >= 10, p_year + 1L, p_year)   # WY year label

#### simple correlational analysis
# 1a) Annual total precip (calendar year) to match veg_prod years
precip_annual <- tapp(precip, index = p_year, fun = sum)
names(precip_annual) <- sort(unique(p_year))

# Ensure years align with veg_prod layers
common_years <- intersect(names(veg_prod), names(precip_annual))
precip_annual <- precip_annual[[common_years]]
veg_prod_aln  <- veg_prod[[common_years]]

# 1b) Seasonal partitions (cool vs. monsoon)
precip_cool <- tapp(precip[[is_cool]],    index = p_year[is_cool],    fun = sum)
precip_mons <- tapp(precip[[is_monsoon]], index = p_year[is_monsoon], fun = sum)

# strip out X in names
names(veg_prod)     <- sub("^X", "", names(veg_prod))
names(precip_cool)  <- sub("^X", "", names(precip_cool))
names(precip_mons)  <- sub("^X", "", names(precip_mons))
names(precip_annual) <- sub("^X", "", names(precip_annual))

precip_cool <- precip_cool[[common_years]]
precip_mons <- precip_mons[[common_years]]

# 1c) Per-pixel correlation with annual totals
corr_ann <- app(c(precip_annual, veg_prod_aln), fun = function(v) {
  n <- length(v) / 2
  x <- v[1:n]; y <- v[(n+1):(2*n)]
  if (all(is.na(x)) || all(is.na(y))) return(NA)
  suppressWarnings(cor(x, y, use = "pairwise.complete.obs"))
})
names(corr_ann) <- "r_annual"
    # plot up correlation
    # Diverging palette (e.g., blue-white-red)
    pal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))

    # Determine color breaks centered at 0
    rng <- range(values(corr_ann), na.rm = TRUE)
    max_abs <- max(abs(rng))
    breaks <- seq(-max_abs, max_abs, length.out = 10)
    
    # Plot correlation map
    plot(corr_ann,
         col = pal(length(breaks) - 1),
         breaks = breaks,
         main = "Correlation: Annual Precipitation vs Vegetation Production",
         axes = FALSE)
    plot(shp, add=TRUE, col=NA, border="black")


# 1d) Per-pixel multiple regression: prod ~ cool + monsoon (R² map)
r2_map <- app(c(precip_cool, precip_mons, veg_prod_aln), fun = function(v){
  n <- length(v)/3
  pc <- v[1:n]; pm <- v[(n+1):(2*n)]; y <- v[(2*n+1):(3*n)]
  if (sum(is.finite(pc) & is.finite(pm) & is.finite(y)) < 8) return(NA)
  df <- data.frame(y, pc, pm)
  fit <- try(lm(y ~ pc + pm, df), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA)
  summary(fit)$r.squared
})
names(r2_map) <- "R2_cool_monsoon"

  # plot r2 map
  plot(r2_map,
       col = rev(terrain.colors(10)),
       main = "R²: Veg Production ~ Cool + Monsoon Precipitation",
       axes = FALSE)
  plot(shp, add=TRUE, col=NA, border="black")

# beta coefficients of predictors
  beta_map <- app(c(precip_cool, precip_mons, veg_prod_aln), fun = function(v){
    n <- length(v)/3
    pc <- v[1:n]; pm <- v[(n+1):(2*n)]; y <- v[(2*n+1):(3*n)]
    ok <- is.finite(pc) & is.finite(pm) & is.finite(y)
    if (sum(ok) < 8) return(c(NA, NA))
    df <- data.frame(y = scale(y[ok]), pc = scale(pc[ok]), pm = scale(pm[ok]))
    fit <- lm(y ~ pc + pm, df)
    c(coef(fit)[2], coef(fit)[3])
  })
  
  names(beta_map) <- c("beta_cool", "beta_monsoon")
  
  # 1 = Cool season dominant, 2 = Monsoon dominant
  dominant <- app(beta_map, fun = function(x){
    if (any(is.na(x))) return(NA)
    if (abs(x[1]) > abs(x[2])) return(1) else return(2)
  })
  
  levels(dominant) <- data.frame(ID = 1:2, season = c("Cool", "Monsoon"))
  
  plot(dominant,
       col = c("steelblue3", "darkorange2"),
       main = "Dominant Precipitation Season (by standardized β)",
       legend = TRUE)
  plot(shp, add=TRUE, col=NA, border="black")
  
  dominant_masked <- mask(dominant, r2_map < 0.2, maskvalues = TRUE)
  plot(dominant_masked,
       col = c("steelblue3", "darkorange2"),
       main = "Dominant Season (β), masked where R² < 0.2")
  plot(shp, add=TRUE, col=NA, border="black")
  
# 2a) Build monthly-to-annual design matrix by year with custom month weights
# Here: three windows as examples — early monsoon (Jun–Jul), late monsoon (Aug–Sep), fall (Oct–Nov)
win_ey <- p_month %in% 6:7
win_lm <- p_month %in% 8:9
win_fa <- p_month %in% 10:11

sum_by <- function(mask, label) {
  out <- tapp(precip[[mask]], index = p_year[mask], fun = sum)
  names(out) <- sort(unique(p_year[mask]))
  out[[common_years]]
}
P_early <- sum_by(win_ey, "early")
P_late  <- sum_by(win_lm, "late")
P_fall  <- sum_by(win_fa, "fall")

# 2b) Pixel-wise regression: prod ~ early + late + fall (Adj R² and beta for 'late')
reg_maps <- app(c(P_early, P_late, P_fall, veg_prod_aln), fun = function(v){
  n <- length(v)/4
  e <- v[1:n]; l <- v[(n+1):(2*n)]; f <- v[(2*n+1):(3*n)]; y <- v[(3*n+1):(4*n)]
  ok <- is.finite(e) & is.finite(l) & is.finite(f) & is.finite(y)
  if (sum(ok) < 8) return(rep(NA, 4))
  fit <- try(lm(y ~ e + l + f), silent = TRUE)
  if (inherits(fit, "try-error")) return(rep(NA, 4))
  coefs <- coef(fit)
  c(summary(fit)$adj.r.squared, coefs["e"], coefs["l"], coefs["f"])
})
names(reg_maps) <- c("AdjR2", "beta_early", "beta_late", "beta_fall")
plot(reg_maps)

library(RColorBrewer)
# symmetric limits
lo <- -3; hi <- 3
nb <- 10
brks <- seq(lo, hi, length.out = nb + 1)
# diverging palette centered at 0
pal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))  # or "BrBG", "PiYG", "RdYlBu"

plot(reg_maps[[2:4]],
     col = pal(nb),
     breaks = brks,
     #main = "β Fall Precip",
     axes = FALSE)
plot(shp, add = TRUE, col = NA, border = "black")

  # plot beta on common colorramp
  # Combine just the beta layers
  betas <- reg_maps[[c("beta_early", "beta_late", "beta_fall")]]
  
  brks <- seq(-3, 3, length.out = 11)
  pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))
  
  # Plot each layer with identical color scale
  par(mfrow = c(1, 3))
  plot(betas$beta_early,  breaks = brks, col = pal(length(brks)-1),
       main = "β Early Monsoon",  axes = FALSE)
        plot(shp, add = TRUE, col = NA, border = "black")
  plot(betas$beta_late,   breaks = brks, col = pal(length(brks)-1),
       main = "β Late Monsoon",   axes = FALSE)
        plot(shp, add = TRUE, col = NA, border = "black")
  plot(betas$beta_fall,   breaks = brks, col = pal(length(brks)-1),
       main = "β Fall Precip",    axes = FALSE)
        plot(shp, add = TRUE, col = NA, border = "black")
  par(mfrow = c(1,1))# end of script
  
##### stratify by vegetation type
# analysis function  
  library(terra)
  library(dplyr)
  library(tidyr)
  library(purrr)
  
  # --- Align rasters ------------------------------------------------------------
  #veg_type_f <- as.factor(veg_type)
  veg_type_f<-veg_type
  
  precip_cool  <- project(precip_cool,  crs(veg_type_f))
  precip_mons  <- project(precip_mons,  crs(veg_type_f))
  veg_prod_aln <- project(veg_prod_aln, crs(veg_type_f))
  
  precip_cool  <- mask(crop(precip_cool,  veg_type_f), veg_type_f)
  precip_mons  <- mask(crop(precip_mons,  veg_type_f), veg_type_f)
  veg_prod_aln <- mask(crop(veg_prod_aln, veg_type_f), veg_type_f)
  
  # --- 1. Zonal (by vegetation type) mean per year ------------------------------
  z_cool <- zonal(precip_cool,  veg_type_f, fun = "mean", na.rm = TRUE)
  z_mons <- zonal(precip_mons,  veg_type_f, fun = "mean", na.rm = TRUE)
  z_prod <- zonal(veg_prod_aln, veg_type_f, fun = "mean", na.rm = TRUE)
  
  # Each is a matrix: first column = veg_type ID, others = yearly means
  # Convert to tidy data frames
  to_long <- function(z, varname){
    as_tibble(z) |>
      rename(type = 1) |>
      pivot_longer(-type, names_to = "year", values_to = varname) |>
      mutate(year = as.integer(gsub("\\D", "", year)))
  }
  
  df_cool <- to_long(z_cool, "Pcool")
  df_mons <- to_long(z_mons, "Pmons")
  df_prod <- to_long(z_prod, "Prod")
  
  # --- 2. Merge to one tidy table -----------------------------------------------
  DF <- reduce(list(df_cool, df_mons, df_prod), full_join, by = c("type","year")) |>
    drop_na()
  
  # --- Attach vegetation labels (robust) ---
  # Extract the lookup table properly
  # lut <- levels(veg_type_f)[[1]]
  # colnames(lut) <- c("ID", "veg_name")
  # 
  # # Ensure matching types for join
  # lut <- lut %>% mutate(type = as.character(ID))
  # 
  # # Join labels onto DF
  # DF <- DF %>%
  #   left_join(lut %>% select(type, veg_name), by = "type")
  
  # Sanity check
  unique(DF$type)
  colnames(DF)[1]<-c("veg_name")
  
  # --- 3. Model fits by vegetation type -----------------------------------------
  mods <- DF |>
    group_by(veg_name) |>
    nest() |>
    mutate(model = map(data, ~lm(Prod ~ Pcool + Pmons, data = .x)),
           tidy  = map(model, broom::tidy),
           glance= map(model, broom::glance))
  
  # --- 4. Summarize coefficients & R² -------------------------------------------
  coef_tbl <- mods |>
    select(veg_name, tidy) |>
    unnest(tidy) |>
    filter(term %in% c("Pcool","Pmons"))
  
  r2_tbl <- mods |>
    select(veg_name, glance) |>
    unnest(glance) |>
    select(veg_name, r.squared, adj.r.squared)
  
  # --- 5. Quick visualization ---------------------------------------------------
  library(ggplot2)
  
  ggplot(coef_tbl, aes(x = veg_name, y = estimate,
                       fill = term)) +
    geom_bar(stat="identity", position="dodge") +
    coord_flip() +
    labs(y = "β coefficient", x = "Vegetation type",
         title = "Seasonal precipitation sensitivity by vegetation type") +
    theme_minimal()
  
  ggplot(r2_tbl, aes(x = reorder(veg_name, adj.r.squared), y = adj.r.squared)) +
    geom_col(fill="steelblue") +
    coord_flip() +
    labs(y = "Adjusted R²", x = "Vegetation type",
         title = "Model fit (Prod ~ Pcool + Pmons)") +
    theme_minimal()
  
  #############
  
  # 4a) Compute API per month with exponential decay, then summarize to year
  api_monthly <- app(precip, fun = function(x) {
    # x is vector over months for a pixel
    k <- 0.9  # memory parameter (try 0.8–0.95)
    y <- rep(NA_real_, length(x))
    for (i in seq_along(x)) {
      if (i == 1) y[i] <- x[i]
      else y[i] <- k*y[i-1] + x[i]
    }
    y
  })
  
  # 4b) Summarize API over critical window (e.g., Aug–Oct) and correlate with production
  api_AugOct <- tapp(api_monthly[[p_month %in% 8:10]], index = p_year[p_month %in% 8:10], fun = mean)
  names(api_AugOct) <- sub("^X", "", names(api_AugOct))
  api_AugOct <- api_AugOct[[common_years]]
  
  api_r2 <- app(c(api_AugOct, veg_prod_aln), fun = function(v){
    n <- length(v)/2
    a <- v[1:n]; y <- v[(n+1):(2*n)]
    if (sum(is.finite(a) & is.finite(y)) < 8) return(NA)
    summary(lm(y ~ a))$r.squared
  })
  names(api_r2) <- "R2_API_AugOct"
  
  # with residuals
  
  api_out <- app(c(api_AugOct, veg_prod_aln), fun = function(v){
    n <- length(v)/2
    api <- v[1:n]
    y   <- v[(n+1):(2*n)]
    ok <- is.finite(api) & is.finite(y)
    if (sum(ok) < 8) return(c(NA, NA))
    fit <- try(lm(y ~ api), silent = TRUE)
    if (inherits(fit, "try-error")) return(c(NA, NA))
    pred <- predict(fit)
    resid <- y - pred
    r2 <- summary(fit)$r.squared
    rmse <- sqrt(mean(resid^2))
    c(r2, rmse)
  })
  names(api_out) <- c("R2_API_AugOct", "RMSE_API_AugOct")
  
  plot(api_out$R2_API_AugOct, col = rev(terrain.colors(20)),
       main = "R2: API Aug–Oct → Production",
       axes = FALSE)
  plot(api_out$RMSE_API_AugOct, col = rev(terrain.colors(20)),
       main = "Model RMSE: API Aug–Oct → Production",
       axes = FALSE)
  plot(shp, add=TRUE, col=NA, border="black")
  

  ############################################################
  # 5) DIMENSIONALITY REDUCTION OF MONTHLY PRECIP SIGNALS
  ############################################################
  
  library(terra)
  library(dplyr)
  library(pls)
  
  set.seed(42)
  
  ## ---------------------------------------------------------
  ## A. Sample pixels within vegetation domain
  ## ---------------------------------------------------------
  
  samp2 <- spatSample(
    veg_type,
    size   = 5000,
    method = "regular",
    as.points = TRUE,
    na.rm  = TRUE
  )
  
  ## ---------------------------------------------------------
  ## B. Extract MONTHLY precip at sample points
  ## ---------------------------------------------------------
  
  P_mon_raw <- terra::extract(precip, samp2, bind = TRUE)
  P_mon     <- as.data.frame(P_mon_raw)
  
  # Keep only monthly precip columns: names like "198401" or "X198401"
  is_month_col <- grepl("^X?[0-9]{6}$", names(P_mon))
  P_mon <- P_mon[, is_month_col, drop = FALSE]
  
  # Standardize names to "YYYYMM"
  names(P_mon) <- sub("^X", "", names(P_mon))
  
  # Build year/month vectors from these names
  pr_names <- names(P_mon)                       # e.g., "198401"
  p_year   <- as.integer(substr(pr_names, 1, 4))
  p_month  <- as.integer(substr(pr_names, 5, 6))
  
  stopifnot(ncol(P_mon) == length(p_year))
  
  precip_years <- sort(unique(p_year))
  
  ## ---------------------------------------------------------
  ## C. Extract ANNUAL vegetation production at same points
  ## ---------------------------------------------------------
  
  Y_ann_raw <- terra::extract(veg_prod_aln, samp2, bind = TRUE)
  Y_ann     <- as.data.frame(Y_ann_raw)
  
  # Keep only annual production columns: "1984", "X1984", etc.
  is_year_col <- grepl("^X?[0-9]{4}$", names(Y_ann))
  Y_ann <- Y_ann[, is_year_col, drop = FALSE]
  
  # Standardize names to "YYYY"
  names(Y_ann) <- sub("^X", "", names(Y_ann))
  
  prod_years <- as.integer(names(veg_prod_aln))
  
  ## ---------------------------------------------------------
  ## D. Define COMMON YEARS (drops 2012 automatically)
  ## ---------------------------------------------------------
  
  common_years <- sort(intersect(precip_years, prod_years))
  # e.g. 1984–2011, 2013–2024
  
  ## ---------------------------------------------------------
  ## E. Align PRECIP to common years
  ## ---------------------------------------------------------
  
  good_idx <- which(p_year %in% common_years)
  
  P_mon_aligned <- P_mon[, good_idx, drop = FALSE]
  p_year        <- p_year[good_idx]
  p_month       <- p_month[good_idx]
  
  stopifnot(ncol(P_mon_aligned) == length(p_year))
  
  ## ---------------------------------------------------------
  ## F. Align PRODUCTION to common years
  ## ---------------------------------------------------------
  
  # Reorder / subset annual production columns to common_years
  Y_ann_aligned <- Y_ann[, match(as.character(common_years), names(Y_ann)), drop = FALSE]
  stopifnot(ncol(Y_ann_aligned) == length(common_years))
  
  ## ---------------------------------------------------------
  ## G. Monthly feature builder (base R; always years × 12)
  ## ---------------------------------------------------------
  
  monthly_by_year <- function(row_vec) {
    ny <- length(common_years)
    mat <- matrix(NA_real_, nrow = ny, ncol = 12,
                  dimnames = list(common_years, as.character(1:12)))
    
    for (i in seq_along(common_years)) {
      yr  <- common_years[i]
      idx <- which(p_year == yr)
      # Expect exactly 12 months per year
      if (length(idx) == 12) {
        mat[i, ] <- row_vec[idx]
      } else {
        # If something is off, leave NA for that year
        mat[i, ] <- NA_real_
      }
    }
    
    mat
  }
  
  ## ---------------------------------------------------------
  ## H. Build list of [years × 12] matrices (one per pixel)
  ## ---------------------------------------------------------
  
  X_list <- lapply(seq_len(nrow(P_mon_aligned)), function(i) {
    monthly_by_year(as.numeric(P_mon_aligned[i, ]))
  })
  
  # Sanity: all matrices should have same dimensions
  dim_tab <- table(sapply(X_list, function(m) paste(dim(m), collapse = "x")))
  print(dim_tab)   # expect "40x12" with count ~5000
  
  ## Production matrix [samples × years]
  Y_mat <- as.matrix(Y_ann_aligned)  # columns correspond to common_years
  
  ## ---------------------------------------------------------
  ## I. Aggregate across pixels for exploratory PLS
  ## ---------------------------------------------------------
  
  # Mean monthly pattern across pixels: [n_years × 12]
  X_year_mean <- Reduce("+", X_list) / length(X_list)
  
  # Mean annual production across pixels: [n_years]
  y_year_mean <- colMeans(Y_mat, na.rm = TRUE)
  
  # Check shapes
  print(dim(X_year_mean))          # expect c(length(common_years), 12)
  print(length(y_year_mean))       # expect length(common_years)
  
  stopifnot(nrow(X_year_mean) == length(y_year_mean))
  
  ## ---------------------------------------------------------
  ## J. PLS regression (mean pattern, 12 monthly predictors)
  ## ---------------------------------------------------------
  
  df_pls <- data.frame(
    Prod = y_year_mean,
    X_year_mean
  )
  
  # Name the 12 month columns M1..M12
  colnames(df_pls)[-1] <- paste0("M", 1:12)
  
  # Ensure we don't ask for more components than allowed
  max_ncomp <- min(3, ncol(df_pls) - 1, nrow(df_pls) - 1)
  
  pls_fit <- plsr(
    Prod ~ .,
    data       = df_pls,
    ncomp      = max_ncomp,
    validation = "LOO"
  )
  
  summary(pls_fit)
  loadings(pls_fit)
  scores(pls_fit)
  
  