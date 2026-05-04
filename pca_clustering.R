############################################
# PCA-based Clustering per Network
############################################

library(readxl)
library(cluster)
library(dplyr)

############################################
# Load data
############################################

data <- read_excel("data/database_parametri.xlsx")
data <- data.frame(lapply(data, function(x)
  if(is.list(x)) unlist(x) else x))

networks <- unique(data$Network)

############################################
# Function: optimal k via Gap Statistic
############################################

calculate_optimal_k <- function(gap_stat, k_min = 2, k_max = 5) {
  
  gap_values <- gap_stat$Tab[k_min:k_max, "gap"]
  se_values  <- gap_stat$Tab[k_min:k_max, "SE.sim"]
  
  n_k <- length(gap_values)
  
  for (i in 1:(n_k - 1)) {
    if (gap_values[i] >= gap_values[i + 1] - se_values[i + 1]) {
      return(i + k_min - 1)
    }
  }
  
  return(which.max(gap_values) + k_min - 1)
}

############################################
# Results container
############################################

cluster_results <- data.frame()
pca_variance <- data.frame()

############################################
# Loop per network
############################################

for (network in networks) {
  
  cat("\n==== PCA clustering:", network, "====\n")
  
  subset_data <- data[data$Network == network, ]
  
  numeric_cols <- sapply(subset_data, is.numeric)
  param_cols <- setdiff(names(subset_data)[numeric_cols],
                        c("ID", "Network"))
  
  param_data <- subset_data[, param_cols, drop = FALSE]
  
  # remove empty columns
  param_data <- param_data[, colSums(!is.na(param_data)) > 0,
                           drop = FALSE]
  
  if (nrow(param_data) < 3) next
  
  ############################################
  # PCA
  ############################################
  
  param_scaled <- scale(param_data)
  
  pca_model <- prcomp(param_scaled,
                      center = TRUE,
                      scale. = TRUE)
  
  variance_explained <- summary(pca_model)$importance[2, 1:2]
  
  pca_variance <- rbind(
    pca_variance,
    data.frame(
      Network = network,
      PC1_variance = variance_explained[1],
      PC2_variance = variance_explained[2]
    )
  )
  
  ############################################
  # PCA scores
  ############################################
  
  scores <- as.data.frame(pca_model$x[,1:2])
  
  ############################################
  # Optimal number of clusters
  ############################################
  
  set.seed(123)
  
  gap_stat <- clusGap(
    scores,
    FUN = kmeans,
    K.max = 5,
    B = 1000,
    nstart = 200
  )
  
  optimal_k <- calculate_optimal_k(gap_stat)
  
  ############################################
  # Final k-means clustering
  ############################################
  
  set.seed(123)
  
  final_kmeans <- kmeans(
    scores,
    centers = optimal_k,
    nstart = 100
  )
  
  ############################################
  # Save cluster assignment
  ############################################
  
  cluster_df <- data.frame(
    ID = subset_data$ID,
    Network = network,
    Cluster = final_kmeans$cluster
  )
  
  cluster_results <- rbind(cluster_results, cluster_df)
}

############################################
# Save outputs
############################################

write.csv(cluster_results,
          "results/pca_clusters.csv",
          row.names = FALSE)

write.csv(pca_variance,
          "results/pca_variance_explained.csv",
          row.names = FALSE)

cat("\n==== ANALYSIS COMPLETED ====\n")