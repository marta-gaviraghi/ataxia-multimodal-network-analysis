############################################
# Regression Analysis Pipeline

# Stepwise linear regression models computed independently
# for each functional network and clinical outcome
############################################

library(readxl)
library(MASS)
library(dplyr)

############################################
# Functions
############################################

# Variance Inflation Factor
calcola_vif <- function(model) {
  X <- model.matrix(model)[, -1]
  vif_values <- diag(solve(cor(X)))
  return(vif_values)
}

# Normalized RMSE
rmse <- function(actual, predicted) {
  rmse_val <- sqrt(mean((actual - predicted)^2, na.rm = TRUE))
  range_y <- max(actual, na.rm = TRUE) - min(actual, na.rm = TRUE)
  rmse_val / range_y
}

############################################
# Load data
############################################

data <- read_excel("data/file_regression.xlsx")
data <- data.frame(lapply(data, function(x) if(is.list(x)) unlist(x) else x))

############################################
# Column selection
############################################

test_cols  <- which(colnames(data) == "SARA"):(ncol(data) - 1)
param_cols <- 3:12

############################################
# Results container
############################################

risultati <- data.frame(
  Test = character(),
  Network = character(),
  R2_adj = numeric(),
  RMSE = numeric(),
  Variabili = character(),
  P_values = character(),
  VIF = character(),
  stringsAsFactors = FALSE
)

############################################
# Regression loop
############################################

for (test_idx in test_cols) {
  
  test_name <- colnames(data)[test_idx]
  
  for (network in unique(data$Network)) {
    
    subset_data <- data[data$Network == network, ]
    
    X <- subset_data[, param_cols]
    y <- subset_data[, test_idx]
    
    reg_data <- data.frame(y = y, X)
    
    tryCatch({
      
      # Full model
      formula <- as.formula(
        paste("y ~", paste(names(X), collapse = " + "))
      )
      
      full_model <- lm(formula, data = reg_data)
      
      # Stepwise selection
      step_model <- stepAIC(full_model,
                            direction = "both",
                            trace = FALSE)
      
      summary_model <- summary(step_model)
      
      predictions <- predict(step_model, newdata = reg_data)
      
      rmse_value <- rmse(y, predictions)
      
      coef_info <- summary_model$coefficients
      
      significant_vars <- rownames(coef_info)[
        coef_info[,4] < 0.05 &
          rownames(coef_info) != "(Intercept)"
      ]
      
      vif_values <- calcola_vif(step_model)
      
      if (length(significant_vars) > 0) {
        
        risultati <- rbind(
          risultati,
          data.frame(
            Test = test_name,
            Network = network,
            R2_adj = summary_model$adj.r.squared,
            RMSE = rmse_value,
            Variabili = paste(significant_vars, collapse = ", "),
            P_values = paste(
              round(coef_info[significant_vars,4],4),
              collapse = ", "
            ),
            VIF = paste(round(vif_values,2), collapse = ", ")
          )
        )
      }
      
    }, error = function(e) {
      cat("ERROR:", conditionMessage(e),
          "Test:", test_name,
          "Network:", network, "\n")
    })
  }
}

############################################
# Save results
############################################

write.csv(risultati,
          "results/regression_results.csv",
          row.names = FALSE)

############################################
# Best models selection
############################################

soglia_rmse <- median(risultati$RMSE, na.rm = TRUE)

migliori_risultati <- risultati[
  risultati$R2_adj > 0.6 &
    risultati$RMSE < soglia_rmse,
]

write.csv(migliori_risultati,
          "results/best_models.csv",
          row.names = FALSE)