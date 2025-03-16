# randomForest analysis of RPMS-climate relationships
# MAC 01/30/25


library(terra)
library(dplyr)
library(tidyr)
library(caret)
library(ggplot2)
library(Metrics)

# Load raster layers
prod_stack <- rast("./data/4K_KNF_RPMS_1984_2024.tif")
precip_stack <- rast("./data/KNF_GridMet_monthly_pr_1984_2024.tif")

# Convert Precipitation Data to Data Frame
precip_df <- as.data.frame(precip_stack, xy = TRUE, na.rm = TRUE) %>%
  pivot_longer(cols = -c(x, y), names_to = "date", values_to = "precip") %>%
  mutate(date = gsub("^X", "", date)) %>%
  separate(date, into = c("year", "month"), sep = "\\.", convert = TRUE) %>%
  mutate(month = as.integer(month), year = as.integer(year))

# Convert Production Data to Data Frame
prod_df <- as.data.frame(prod_stack, xy = TRUE, na.rm = TRUE) %>%
  pivot_longer(cols = -c(x, y), names_to = "year", values_to = "production") %>%
  mutate(year = as.integer(year))

# Ensure Each (x, y, year) Has 12 Months
precip_df <- precip_df %>%
  complete(x, y, year, month = 1:12, fill = list(precip = 0))

# Merge Precipitation and Production Data
merged_df <- precip_df %>%
  left_join(prod_df, by = c("x", "y", "year")) %>%
  group_by(x, y, year) %>%
  arrange(month, .by_group = TRUE) %>%  # Ensure correct order
  mutate(
    precip_total = sum(precip, na.rm = TRUE),
    jfm_precip = sum(precip[month %in% c(1:3)], na.rm = TRUE),
    amj_precip = sum(precip[month %in% c(4:6)], na.rm = TRUE),
    jas_precip = sum(precip[month %in% c(7:9)], na.rm = TRUE),
    ond_precip = sum(precip[month %in% c(10:12)], na.rm = TRUE),
    cum_precip = cumsum(precip),  
    lag_precip = lag(precip, default = 0) 
  ) %>%
  ungroup() %>%
  na.omit()

# Ensure Month is a Factor
merged_df <- merged_df %>%
  mutate(month = factor(month, levels = 1:12))

# Spatially Informed Train-Test Split
set.seed(123)
sampled_locs <- merged_df %>% distinct(x, y) %>% sample_frac(0.1)
train_data <- merged_df %>% semi_join(sampled_locs, by = c("x", "y"))
test_data <- anti_join(merged_df, train_data, by = c("x", "y", "year", "month"))

# Train Model
model_formula <- production ~ precip_total + jfm_precip + amj_precip + jas_precip + ond_precip + cum_precip + lag_precip + month

rf_model <- train(
  model_formula,
  data = train_data,
  method = "ranger",
  trControl = trainControl(method = "cv", number = 5),
  tuneLength = 3,
  num.trees = 10,
  importance = "impurity"  # Enables variable importance calculation
)

# Evaluate Model
predictions <- predict(rf_model, newdata = test_data)
results <- data.frame(Observed = test_data$production, Predicted = predictions)

# Compute Evaluation Metrics
correlation <- cor(results$Observed, results$Predicted)
rmse_value <- rmse(results$Observed, results$Predicted)
mae_value <- mae(results$Observed, results$Predicted)

cat("Correlation:", correlation, "\n")
cat("RMSE:", rmse_value, "\n")
cat("MAE:", mae_value, "\n")

#####

library(reshape2)

# Compute correlation matrix
cor_matrix <- cor(train_data %>% select(-x, -y, -year, -month))

# Convert to long format for plotting
cor_melted <- melt(cor_matrix)

# Plot correlation heatmap
ggplot(cor_melted, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Correlation Heatmap of Predictors")


library(ggplot2)
library(viridis)  # For a continuous color scale

# Calculate prediction error and add it to the test dataset
test_data <- test_data %>%
  mutate(Predicted = predictions,
         Error = production - Predicted)

# Optionally, filter by a specific year/month if needed:
error_map_df <- test_data %>% filter(year == 2020, month == "1")

# Create the spatial error map:
ggplot(error_map_df, aes(x = x, y = y, fill = Error)) +
  geom_raster() +  # Assumes the grid is regular; use geom_tile() otherwise.
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  labs(
    title = "Spatial Map of Prediction Error",
    x = "Longitude",
    y = "Latitude",
    fill = "Error"
  ) +
  theme_minimal()

#####
# map by pixel
library(dplyr)

# Calculate the error and aggregate by pixel (x, y)
avg_error_df <- test_data %>%
  mutate(Error = production - predictions) %>%
  group_by(x, y) %>%
  summarize(Avg_Error = mean(Error, na.rm = TRUE)) %>%
  ungroup()

library(ggplot2)
library(viridis)  # Optional, for a pleasing color scale

ggplot(avg_error_df, aes(x = x, y = y, fill = Avg_Error)) +
  geom_raster() +  # Use geom_tile() if your grid isn't perfectly regular
  scale_fill_viridis_c(option = "plasma", direction = -1) +
  labs(
    title = "Average Prediction Error by Pixel",
    x = "Longitude",
    y = "Latitude",
    fill = "Avg Error"
  ) +
  theme_minimal()


