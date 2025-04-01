
library(terra)   # For working with SpatRaster
library(dplyr)   # Data manipulation
library(tidyr)   # Data reshaping
library(lme4)    # Linear mixed-effects models
library(mgcv)    # Generalized Additive Models
library(ranger)  # Random forest modeling
library(ggplot2) # Visualization
library(purrr)   # Functional programming

# Load SpatRaster data
veg_prod <- rast("./data/processed/4K_KNF_RPMS_1984_2024.tif")  # Annual vegetation production
precip <- rast("./data/processed/KNF_GridMet_monthly_pr_1984_2024.tif")  # Monthly precipitation
veg_type <- rast("./data/processed/4K_KNF_200BPS_GROUPVEG.tif")  # Vegetation type
activeCat(veg_type)<-5

# Create a vector of years
years <- 1984:2024

# Aggregate monthly precipitation to annual precipitation
annual_precip <- tapp(precip, rep(1:41, each=12), fun=sum)
names(annual_precip) <- years

# Compute antecedent precipitation patterns
antecedent_precip <- function(precip, months) {
  tapply(precip, rep(1:(492 - months), each=months, length.out=492), sum, na.rm=TRUE)
}

precip_3mo <- tapp(precip, rep(1:(41*12 - 3), each=3, length.out=492), fun=sum)
precip_6mo <- tapp(precip, rep(1:(41*12 - 6), each=6, length.out=492), fun=sum)
precip_12mo <- annual_precip  # 12-month antecedent precipitation

# Convert SpatRasters to Data Frames
veg_prod_df <- as.data.frame(veg_prod, xy=TRUE) %>% pivot_longer(-c(x, y), names_to="year", values_to="veg_prod")
precip_df <- as.data.frame(precip_12mo, xy=TRUE) %>% pivot_longer(-c(x, y), names_to="year", values_to="precip_12mo")
precip_6mo_df <- as.data.frame(precip_6mo, xy=TRUE) %>% pivot_longer(-c(x, y), names_to="year", values_to="precip_6mo")
precip_3mo_df <- as.data.frame(precip_3mo, xy=TRUE) %>% pivot_longer(-c(x, y), names_to="year", values_to="precip_3mo")
veg_type_df <- as.data.frame(veg_type, xy=TRUE) %>% rename(veg_type = GROUPVEG)

# Merge all data into one dataframe
df <- full_join(veg_prod_df, precip_df, by=c("x", "y", "year")) %>%
  full_join(precip_6mo_df, by=c("x", "y", "year")) %>%
  full_join(precip_3mo_df, by=c("x", "y", "year")) %>%
  full_join(veg_type_df, by=c("x", "y"))

df <- df %>% mutate(year = as.numeric(year))

df <- df %>% drop_na(veg_prod, precip_12mo, precip_6mo, precip_3mo, veg_type)

# Exploratory Data Analysis
summary(df)
ggplot(df, aes(x=precip_12mo, y=veg_prod, color=veg_type)) +
  geom_point(alpha=0.5) + geom_smooth(method="lm") +
  labs(title="Vegetation Production vs. Precipitation", x="12-Month Precip", y="Vegetation Production")

# Linear Mixed Model
lmm <- lmer(veg_prod ~ precip_12mo + precip_6mo + precip_3mo + (1 | veg_type) + (1 | x/y), data=df)
summary(lmm)

# Generalized Additive Model
gam_model <- gam(veg_prod ~ s(precip_12mo) + s(precip_6mo) + s(precip_3mo) + veg_type, data=df)
summary(gam_model)

# Random Forest Model
rf_model <- ranger(veg_prod ~ precip_12mo + precip_6mo + precip_3mo + veg_type, data=df, importance="impurity")
print(rf_model)

# Variable Importance
importance(rf_model)

# Model Validation
df$pred_lmm <- predict(lmm, df)
df$pred_gam <- predict(gam_model, df)
df$pred_rf <- predict(rf_model, df)

# Compare Predictions
rmse <- function(actual, predicted) sqrt(mean((actual - predicted)^2, na.rm=TRUE))
cat("LMM RMSE:", rmse(df$veg_prod, df$pred_lmm), "\n")
cat("GAM RMSE:", rmse(df$veg_prod, df$pred_gam), "\n")
cat("RF RMSE:", rmse(df$veg_prod, df$pred_rf), "\n")

# Visualization of Predictions
ggplot(df, aes(x=veg_prod, y=pred_lmm)) + geom_point() + geom_smooth(method="lm") + labs(title="LMM Predictions")
ggplot(df, aes(x=veg_prod, y=pred_gam)) + geom_point() + geom_smooth(method="lm") + labs(title="GAM Predictions")
ggplot(df, aes(x=veg_prod, y=pred_rf)) + geom_point() + geom_smooth(method="lm") + labs(title="RF Predictions")
