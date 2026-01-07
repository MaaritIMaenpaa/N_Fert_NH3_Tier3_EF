# Packages ----
library(lme4)

# Data and model ----
# Replace with data and model to use for validation
data <- NH3predict_mALL
model_formula <- formula(mALL_lmer)

# Leave-One-Out Cross-Validation ----
n <- nrow(data)
predictions <- numeric(n)

for (i in 1:n) {
  # training data selection
  train_data <- data[-i, ]
  test_data <- data[i, , drop = FALSE]
  
  # training model
  loo_model <- lmer(model_formula, data = train_data)
  
  # Predict for left-out observation
  predictions[i] <- predict(loo_model, 
                            newdata = test_data, 
                            re.form = NA)
}

# Compare predictions to actual values
actual <- sqrt(data$NH3loss + 4)
results <- data.frame(actual = actual, predicted = predictions)
saveRDS(results, "output/LOO_model_A1Base")

results <- readRDS("output/LOO_model_A1Base")

# Performance metrics

# MSE
mse <-mean((results$actual - results$predicted)^2)
mse 

# RMSE
rmse <- sqrt(mse)
rmse/mean(NH3predict$NH3loss)

# Correlation
cor(results$actual, results$predicted)
ggplot(results, aes(x=actual, y=predicted)) +
  geom_point()

#
#
#

# END ----

