# experimental ML code on monthly precipitation and vegetation production data
# MAC 3/18/25

library(terra)   # For spatial data
library(dplyr)   # For data manipulation
library(tidyr)   # For reshaping data
library(ranger)  # Random forest modeling
library(ggplot2) # Visualization
# library(viridis) # For color scales
# library(xgboost) # For XGBoost modeling
# library(mgcv)    # For Generalized Additive Models (GAMs)
# library(GWmodel) # For Geographically Weighted Regression (GWR)
# library(stringr) # For string manipulation
library(caret)
library(SPEI)

terraOptions(progress=1)

# Load SpatRaster data
veg_prod <- rast("./data/processed/4K_KNF_RPMS_1984_2024.tif")  # Annual vegetation production
precip <- rast("./data/processed/KNF_GridMet_monthly_pr_1984_2024.tif")  # Monthly precipitation
veg_type <- rast("./data/processed/4K_KNF_200BPS_GROUPVEG.tif")  # Vegetation type
activeCat(veg_type)<-5

# set new extent to all spatRast
newExt<-ext(veg_prod)
newExt<-ext(newExt[1],newExt[2],newExt[3],newExt[4]-1)
  veg_prod<-crop(veg_prod,newExt)
  precip<-crop(precip,newExt)
  veg_type<-crop(veg_type,newExt)

# mask out veg_type by veg_prod
veg_type<-mask(veg_type, veg_prod[[40]])
  
# Extract years and months
years <- 1984:2024
months <- rep(1:12, length.out=492)

# Compute rolling 3-month seasonal precipitation for each pixel over time
rolling_sum <- function(x) {
  zoo::rollapply(x, width=3, FUN=sum, align='right', fill=NA, na.rm=TRUE)
}
precip_3mo <- app(precip, rolling_sum)# Compute rolling 3-month seasonal precipitation using tapp() in terra
names(precip_3mo)<-names(precip)
# calculate spi
funSPI <- function(x, scale=y, na.rm=TRUE,...) as.numeric((SPEI::spi(x, scale=scale, 
                                                                     na.rm=na.rm, verbose = FALSE, ...))$fitted)
spi3 <- app(precip, fun=function(x) funSPI(x,3))
names(spi3)<-names(precip)
#  limit values to -3 to 3
spi3<-clamp(spi3,-3,3)

##### convert RPMS to different anomalies
## Std anomalies
# Compute the mean across all years (pixel-wise)
mean_veg <- app(veg_prod, mean, na.rm=TRUE)
# Compute the standard deviation across all years (pixel-wise)
sd_veg <- app(veg_prod, sd, na.rm=TRUE)
# Compute the standardized values (Z-score: (x - mean) / sd)
veg_prod <- (veg_prod - mean_veg) / sd_veg
veg_prod <- clamp(veg_prod, -3, 3)  # Limit values to -3 to 3
##


##### convert to data frames
# Convert SpatRasters to Data Frames
veg_prod_df <- as.data.frame(veg_prod, xy=TRUE) %>% pivot_longer(-c(x, y), names_to="year", values_to="veg_prod")
  veg_prod_df$year <- as.numeric(veg_prod_df$year)
#precip_df <- as.data.frame(precip, xy=TRUE) %>% pivot_longer(-c(x, y), names_to="month_year", values_to="precip")
precip_3mo_df <- as.data.frame(spi3, xy=TRUE) %>% pivot_longer(-c(x, y), names_to="month_year", values_to="precip_3mo") %>%
  mutate(year = as.numeric(substr(month_year, 2, 5)),
         month = as.numeric(substr(month_year, 7, 8)))
veg_type_df <- as.data.frame(veg_type, xy=TRUE) %>% rename(veg_type = GROUPVEG)

# Extract year and month
# precip_df <- precip_df %>% mutate(year = as.numeric(substr(month_year, 2, 5)),
#                                   month = as.numeric(substr(month_year, 7, 8)))
precip_3mo_df <- precip_3mo_df %>% mutate(year = as.numeric(substr(month_year, 2, 5)))

# Spread monthly precipitation into separate columns
# precip_wide <- precip_df %>%
#   pivot_wider(names_from = month, values_from = precip, names_prefix = "precip_")

precip_3mo_wide <- precip_3mo_df %>%
  select(-month_year) %>%
  pivot_wider(names_from = month, values_from = precip_3mo, names_prefix = "precip3mo_")

# Merge all data
df <- full_join(veg_prod_df, precip_3mo_wide, by=c("x", "y", "year")) %>%
  full_join(veg_type_df, by=c("x", "y"))

# Clean Data
df <- df %>% drop_na()
df <- df %>% mutate(veg_type = as.factor(veg_type))

# clean up memory
rm(veg_prod, precip_3mo, veg_type, veg_prod_df, precip_3mo_df, veg_type_df)
gc()

# keep subset of veg_type in df
df<-df %>% filter(veg_type %in% c("Shrubland", "Conifer", "Grassland"))

#### models by veg type
rf_models <- df %>%
  group_split(veg_type) %>%
  setNames(unique(df$veg_type)) %>%
  lapply(function(df_sub) {
    ranger(veg_prod ~ ., data = df_sub %>% select(-x, -y, -year, -veg_type),
           importance="impurity", 
           num.trees=500,  # Reduce trees for faster execution
           #mtry=2,  # Consider fewer variables per split
           #write.forest=FALSE,  # Disable tree storage to save memory
           #min.node.size=10,  # Prevent deep trees for speedup
           num.threads=0,
           keep.inbag = TRUE)  # Use multi-threading for parallel execution
  })

# Compare R² values
sapply(rf_models, function(model) model$r.squared)

# Compare top 5 most important variables for each vegetation type
importance_list <- lapply(rf_models, function(model) {
  as.data.frame(model$variable.importance) %>%
    tibble::rownames_to_column("Variable") %>%
    arrange(desc(model$variable.importance))
})
lapply(importance_list, head, 5)  # Show top 5 most important variables

# Convert importance_list into a single dataframe
importance_df <- bind_rows(
  lapply(names(importance_list), function(veg) {
    importance_list[[veg]] %>%
      mutate(veg_type = veg)  # Add vegetation type as a column
  })
)

# Rename columns for clarity
colnames(importance_df) <- c("Variable", "Importance", "veg_type")

importance_df <- importance_df %>%
  group_by(veg_type) %>%
  mutate(Month = as.numeric(gsub("precip3mo_", "", Variable))) %>%
  mutate(Importance = Importance / max(Importance) * 100)

ggplot(importance_df, aes(y = reorder(Variable, Month), x = Importance, fill = veg_type)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  facet_wrap(~veg_type, ncol=1) +
  scale_fill_viridis_d() +
  labs(title = "Normalized Variable Importance by Vegetation Type",
       y = "Month (Precipitation Lag)",
       x = "Relative Importance (Scaled)") +
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# partial dependence plots
library(pdp)
# Choose a vegetation type
veg_choice <- "Shrubland"
  # Extract model
  rf_model <- rf_models[[veg_choice]]
  # Subset the data for this vegetation type
  df_sub <- df %>% filter(veg_type == veg_choice)
  # Compute PDP for a single predictor (e.g., precip3mo_6)
  pdp_precip <- partial(rf_model, pred.var = "precip3mo_8", train = df_sub, grid.resolution = 100, progress = TRUE)
  # Convert to dataframe and plot
  pdp_precip_df <- as.data.frame(pdp_precip)
  
  ggplot(pdp_precip_df, aes(x = precip3mo_8, y = yhat)) +
    geom_line() +
    labs(title = paste("PDP for precip3mo_8 in", veg_choice),
         x = "Jun-Jul-Aug 3mo SPI", y = "Predicted Vegetation Production") +
    theme_minimal()

##### Gradient Boosting approach
  
  library(xgboost)
  library(dplyr)
  
  # Split data by vegetation type and train separate XGBoost models
  xgb_models <- df %>%
    group_split(veg_type) %>%
    setNames(unique(df$veg_type)) %>%
    lapply(function(df_sub) {
      X <- as.matrix(df_sub %>% select(-x, -y, -year, -veg_type, -veg_prod))  # Features
      y <- df_sub$veg_prod  # Target variable
      
      xgboost(
        data = X,
        label = y,
        objective = "reg:squarederror",  # Regression task
        nrounds = 100,
        eta = 0.05,
        max_depth = 6,
        subsample = 0.8,
        colsample_bytree = 0.8,
        eval_metric = "rmse",
        verbose = 1  # 0=Silence output
      )
    })
  
  # Extract variable importance for each vegetation type
  importance_list <- lapply(names(xgb_models), function(veg) {
    model <- xgb_models[[veg]]
    importance_df <- xgb.importance(model = model, feature_names = colnames(model$data))
    importance_df <- importance_df %>% mutate(veg_type = veg)  # Add vegetation type column
    return(importance_df)
  })
  
  # Combine into one dataframe for plotting
  importance_df <- bind_rows(importance_list)
  
  importance_df <- importance_df %>%
    mutate(Month = as.numeric(gsub("precip3mo_", "", Feature)))  # Extract month number
  
  # Plot with correct chronological order
  ggplot(importance_df, aes(y = reorder(Feature, Month), x = Gain, fill = veg_type)) +
    geom_col(show.legend = FALSE) +
    coord_flip() +
    facet_wrap(~veg_type, ncol=1) +
    scale_fill_viridis_d() +
    labs(title = "XGBoost Variable Importance by Vegetation Type",
         x = "Month (Precipitation Lag)",
         y = "Importance (Gain)") +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # SHAP values
  library(shapviz)
  
  shap_results <- lapply(names(xgb_models), function(veg) {
    model <- xgb_models[[veg]]
    
    # Get feature matrix for this vegetation type
    df_sub <- df %>% filter(veg_type == veg) %>% select(-x, -y, -year, -veg_type, -veg_prod)
    X <- as.matrix(df_sub)
    
    # Compute SHAP values
    shap_values <- shapviz(model, X)
    
    return(list(veg_type = veg, shap_values = shap_values))
  })
  
  shap_plots <- lapply(shap_results, function(result) {
    shapviz_obj <- result$shap_values
    veg_type <- result$veg_type
    
    sv_importance(shapviz_obj) + ggtitle(paste("SHAP Importance for", veg_type))
  })
  
  # show plots
  shap_plots[[1]]
  shap_plots[[2]]  
  shap_plots[[3]]
    
  # Extract and reshape SHAP values
  shap_df <- bind_rows(lapply(shap_results, function(result) {
    shapviz_obj <- result$shap_values
    veg_type <- result$veg_type
    
    # Extract SHAP values from S matrix
    shap_long <- as.data.frame(shapviz_obj$S) %>%
      mutate(Observation = row_number()) %>%
      pivot_longer(cols = -Observation, names_to = "Feature", values_to = "SHAP") %>%
      mutate(veg_type = veg_type)
    
    return(shap_long)
  }))
  
  # shap_df <- shap_df %>%
  #   mutate(Month = as.numeric(gsub("precip3mo_", "", Feature)))  # Extract month number
  # 
  # ggplot(shap_df, aes(x = reorder(Feature, Month), y = SHAP, color = veg_type)) +
  #   geom_jitter(alpha = 0.3, width = 0.2) +
  #   facet_wrap(~veg_type, scales = "free_y") +
  #   labs(
  #     title = "SHAP Values for Precipitation Variables by Vegetation Type",
  #     x = "Precipitation Month",
  #     y = "SHAP Value (Feature Contribution)"
  #   ) +
  #   theme_minimal()
  
  shap_df_sampled <- shap_df %>%
    group_by(veg_type, Feature) %>%
    slice_sample(n = 500) %>%  # Adjust sample size as needed
    ungroup()
  
  ggplot(shap_df_sampled, aes(x = reorder(Feature, as.numeric(gsub("precip3mo_", "", Feature))), y = SHAP)) +
    geom_boxplot(outlier.alpha = 0.1, linewidth = 0.3) +
    facet_wrap(~veg_type, scales = "free_y") +
    labs(
      title = "SHAP Value Distribution per Precipitation Month (Boxplot)",
      x = "Precipitation Month",
      y = "SHAP Value"
    ) +
    theme_minimal()+
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  
  #✔ This shows how each precipitation variable contributes to predictions across different vegetation types.
  #✔ A large positive SHAP value means high precipitation in that month increases vegetation production.
  # Positive SHAP value → Feature increases predicted vegetation production.
  # Negative SHAP value → Feature decreases the prediction.
  # Spread (vertical) → Indicates variability in how that feature affects different observations.
  # Color → Useful if you group by vegetation type or another category.
  

  # Compute mean absolute SHAP values per feature and vegetation type
  shap_summary <- shap_df %>%
    group_by(veg_type, Feature) %>%
    summarise(
      mean_SHAP = mean(SHAP, na.rm = TRUE),   # Average SHAP value (directional)
      abs_mean_SHAP = mean(abs(SHAP), na.rm = TRUE),  # Importance metric
      .groups = "drop"
    ) %>%
    mutate(Month = as.numeric(gsub("precip3mo_", "", Feature)))  # Extract numeric month
  
  # View summary table
  print(shap_summary)
  
  ggplot(shap_summary, aes(x = reorder(Feature, Month), y = mean_SHAP, fill = veg_type)) +
    geom_col(show.legend = TRUE, position = "dodge") +
    facet_wrap(~veg_type, scales = "free_y") +
    labs(
      title = "Directional Mean SHAP Values by Precipitation Month and Vegetation Type",
      x = "Precipitation Month",
      y = "Mean SHAP (Effect on Vegetation Production)"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ###### simple EDA
  # Compute correlation between veg_prod and each monthly precipitation variable, grouped by vegetation type
  correlations <- df %>%
    dplyr::select(veg_type, veg_prod, starts_with("precip3mo_")) %>%
    group_by(veg_type) %>%
    summarise(across(starts_with("precip3mo_"), ~ cor(.x, veg_prod, use = "complete.obs"), .names = "cor_{.col}")) %>%
    pivot_longer(cols = starts_with("cor_"), names_to = "Month", values_to = "Correlation") %>%
    mutate(Month = as.numeric(gsub("cor_precip3mo_", "", Month)))  # Extract numeric month
  
  # View results
  print(correlations)
  
  ggplot(correlations, aes(x = as.factor(Month), y = veg_type, fill = Correlation)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    labs(title = "Correlation Between Monthly Precipitation & Vegetation Production",
         x = "Month",
         y = "Vegetation Type") +
    theme_minimal()
  
  # Run linear regression models by vegetation type
  lm_results <- df %>%
    group_by(veg_type) %>%
    do(model = lm(veg_prod ~ ., data = select(., -x, -y, -year, -veg_type)))
  
  # Extract model summaries
  lm_summaries <- lapply(lm_results$model, summary)
  
  # View R² and coefficients for each vegetation type
  lm_summaries
  
  # plot R2 value
  r2_values <- data.frame(
    veg_type = lm_results$veg_type,
    R2 = sapply(lm_summaries, function(x) x$r.squared)
  )
  
  # Plot R² values
  ggplot(r2_values, aes(x = veg_type, y = R2, fill = veg_type)) +
    geom_col(show.legend = FALSE) +
    coord_flip() +
    labs(title = "Regression R² by Vegetation Type",
         x = "Vegetation Type",
         y = "R² (Model Fit)") +
    theme_minimal()
  
  # screen for interactions
  # library(dplyr)
  # library(broom)  # For tidy model summaries
  # 
  # # Create interaction terms for all precipitation month pairs
  # df_interactions <- df %>%
  #   mutate(across(starts_with("precip3mo_"), scale)) %>%  # Standardize
  #   mutate(
  #     precip3mo_1_7 = precip3mo_1 * precip3mo_7,  # Example: June x July interaction
  #     precip3mo_2_9 = precip3mo_2 * precip3mo_9,  # Example: March x September interaction
  #     precip3mo_3_8 = precip3mo_3 * precip3mo_8   # Example: May x August interaction
  #   )
  # 
  # # Fit model including interaction terms
  # interaction_model <- lm(veg_prod ~ ., data = df_interactions %>% select(-x, -y, -year, -veg_type))
  # summary(interaction_model)
  
  ##### screening variables
 
  library(glmnet)
  
  # Function to fit LASSO regression per vegetation type
  fit_lasso_model <- function(df_sub) {
    df_sub <- as_tibble(df_sub)
    
    # Select only the precipitation predictors
    X_df <- dplyr::select(df_sub, starts_with("precip3mo_"))
    X <- as.matrix(X_df)
    y <- df_sub$veg_prod
    
    lasso_model <- cv.glmnet(X, y, alpha = 1, nfolds = 5)
    best_lambda <- lasso_model$lambda.min
    
    lasso_final <- glmnet(X, y, alpha = 1, lambda = best_lambda)
    
    return(coef(lasso_final))
  }
  
  # Apply to each veg type
  lasso_results <- df %>%
    group_split(veg_type, .keep = FALSE) %>%
    setNames(unique(df$veg_type)) %>%
    lapply(fit_lasso_model)
  
  # Print results
  lasso_results
  
  library(tibble)
  library(dplyr)
  
  # Convert each LASSO coefficient object into a tidy dataframe
  lasso_summary <- bind_rows(
    lapply(names(lasso_results), function(veg) {
      coef_df <- as.matrix(lasso_results[[veg]]) %>%
        as.data.frame() %>%
        rownames_to_column("Feature")
      
      # Rename the coefficient column (always the second column)
      names(coef_df)[2] <- "Coefficient"
      
      coef_df %>%
        filter(Coefficient != 0) %>%
        mutate(veg_type = veg)
    }),
    .id = NULL
  )
  
  lasso_summary <- lasso_summary %>%
    mutate(Month = as.numeric(gsub("precip3mo_", "", Feature))) %>%
    filter(Feature != "(Intercept)")
  
  ggplot(lasso_summary, aes(x = reorder(Feature, Month), y = Coefficient, fill = veg_type)) +
    geom_col(position = "dodge") +
    facet_wrap(~veg_type, scales = "free_y") +
    labs(
      title = "LASSO-Selected Precipitation Predictors by Vegetation Type",
      x = "Precipitation Month",
      y = "LASSO Coefficient"
    ) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  #### Stepwise regression
  fit_limited_stepwise_model <- function(df_sub, max_vars = 3) {
    df_sub <- as_tibble(df_sub) %>%
      dplyr::select(veg_prod, starts_with("precip3mo_"))  # Keep only precipitation predictors
    
    # Fit full model
    full_model <- lm(veg_prod ~ ., data = df_sub)
    
    # Run stepwise selection
    step_model <- stepAIC(full_model, direction = "both", trace = FALSE)
    
    # Extract selected variables
    selected_vars <- broom::tidy(step_model) %>%
      filter(term != "(Intercept)") %>%
      arrange(desc(abs(estimate))) %>%  # Sort by absolute effect size
      head(max_vars)  # Keep only top 2-3 predictors
    
    # Refit final model using only the best predictors
    final_formula <- as.formula(
      paste("veg_prod ~", paste(selected_vars$term, collapse = " + "))
    )
    final_model <- lm(final_formula, data = df_sub)
    
    # Return final model and details
    return(list(
      model = final_model,
      summary = broom::tidy(final_model),
      R2 = summary(final_model)$r.squared,
      AIC = AIC(final_model)
    ))
  }
  
  best_limited_models <- df %>%
    group_split(veg_type, .keep = FALSE) %>%
    setNames(unique(df$veg_type)) %>%
    lapply(fit_limited_stepwise_model, max_vars = 3)  # Choose top 3 predictors
  
  r2_aic_limited <- tibble(
    veg_type = names(best_limited_models),
    R2 = sapply(best_limited_models, function(x) x$R2),
    AIC = sapply(best_limited_models, function(x) x$AIC)
  )
  
  # View results
  print(r2_aic_limited)
  
  # Plot R² & AIC for model comparison
  ggplot(r2_aic_limited, aes(x = veg_type)) +
    geom_col(aes(y = R2, fill = "R²"), position = "dodge") +
    geom_col(aes(y = -AIC / 100, fill = "AIC"), position = "dodge") +  # Scale AIC for visualization
    labs(title = "Best 2-3 Predictor Model Performance by Vegetation Type",
         y = "Higher R² (Better Fit) / Lower AIC (Simpler Model)",
         x = "Vegetation Type") +
    theme_minimal() +
    scale_fill_manual(values = c("R²" = "blue", "AIC" = "red"))
  
  best_limited_summary <- bind_rows(
    lapply(names(best_limited_models), function(veg) {
      best_limited_models[[veg]]$summary %>%
        mutate(veg_type = veg)  # Explicitly add vegetation type
    }),
    .id = NULL
  )
  
  best_limited_summary <- best_limited_summary %>%
    mutate(Month = as.numeric(gsub("precip3mo_", "", term))) %>%  # Extract month number
    arrange(veg_type, Month)  # Ensure proper order within each vegetation type
  
  ggplot(best_limited_summary, aes(x = reorder(term, Month), y = estimate, fill = veg_type)) +
    geom_col(position = "dodge") +
    facet_wrap(~veg_type) +
    labs(title = "Best 3 Predictors by Vegetation Type - Stepwise Regression",
         x = "Precipitation Month",
         y = "Regression Coefficient") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  