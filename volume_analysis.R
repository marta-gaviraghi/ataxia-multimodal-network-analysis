############################################################
# Ataxia Network Volume Analysis
#
# Statistical pipeline for group comparison of regional
# brain volumes across networks (somatomotor and ventral attention networks).
#
# Groups:
# HC = Healthy Controls
# JS = Joubert Syndrome
# SP = Spastic Paraplegia
#
# Pipeline:
# 1. Shapiro–Wilk normality test (per group)
# 2. ANOVA (normal data) or Kruskal–Wallis (non-normal data)
# 3. Post-hoc testing:
#       - Tukey HSD (ANOVA)
#       - Dunn test with Holm correction (Kruskal–Wallis)
# 4. Boxplots generated only for significant regions
#
# NOTE:
# This script represents the statistical workflow used in
# the manuscript. Example datasets are provided for
# reproducibility purposes.
############################################################


############################
# Load required libraries
############################

library(ggplot2)
library(dplyr)
library(dunn.test)
library(readxl)
library(tidyr)


############################
# Define input data paths
############################
# Example data included in repository

network_files <- c(
  "somatomotor" = "example_data/example_volumes.xlsx",
  "ventral_att" = "example_data/example_volumes.xlsx"
)


############################
# Define groups and colors
############################

groups <- c("HC", "JS", "SP")

colors <- c(
  "#2E86C1",  # HC
  "#E67E22",  # JS
  "#8E44AD"   # SP
)


############################
# Output directory
############################

output_base_dir <- "results"
dir.create(output_base_dir, showWarnings = FALSE)


############################
# Utility function:
# p-value formatting
############################

format_pvalue <- function(p) {
  if (p < 0.001) return("p < 0.001")
  if (p < 0.01) return(sprintf("p = %.3f", p))
  return(sprintf("p = %.2f", p))
}

############################
# Dataframe collecting
# all network results
############################

all_results_df <- data.frame(
  network = character(),
  region = character(),
  test_used = character(),
  overall_p = numeric(),
  significant = logical(),
  significant_comparisons = character(),
  stringsAsFactors = FALSE
)

############################################################
# MAIN ANALYSIS LOOP
############################################################

for (network_name in names(network_files)) {
  
  cat("\n===========================================================\n")
  cat("ANALYSIS NETWORK:", network_name, "\n")
  cat("===========================================================\n")
  
  network_output_dir <- file.path(output_base_dir, network_name)
  dir.create(network_output_dir, showWarnings = FALSE)
  
  file_path <- network_files[network_name]
  data <- read_excel(file_path)
  
  data$GROUP <- factor(data$GROUP, levels = groups)
  
  region_columns <- colnames(data)[
    !colnames(data) %in% c("subject_names", "GROUP")
  ]
  
  network_results_df <- data.frame(
    network = character(),
    region = character(),
    test_used = character(),
    overall_p = numeric(),
    significant = logical(),
    significant_comparisons = character(),
    stringsAsFactors = FALSE
  )
  
  total_regions <- length(region_columns)
  significant_regions_count <- 0
  
  
  ##########################################################
  # REGION-WISE ANALYSIS
  ##########################################################
  
  for (region in region_columns) {
    
    cat("\n--------------------------\n")
    cat("Region:", region, "\n")
    cat("--------------------------\n")
    
    ####################################
    # Normality testing
    ####################################
    
    normality_results <- sapply(groups, function(grp) {
      
      group_data <- data %>%
        filter(GROUP == grp) %>%
        pull(region)
      
      normality_test <- shapiro.test(group_data)
      
      cat(
        paste0(
          "Shapiro-Wilk p-value (", grp, "): ",
          normality_test$p.value, "\n"
        )
      )
      
      return(normality_test$p.value >= 0.05)
    })
    
    all_normal <- all(normality_results)
    
    significant_comparisons <- data.frame(
      comparisons = character(),
      p_values = numeric(),
      stringsAsFactors = FALSE
    )
    
    overall_p <- NA
    test_used <- ""
    
    
    ####################################
    # Statistical testing
    ####################################
    
    if (!all_normal) {
      
      cat("Non-normal distribution detected → Kruskal-Wallis test\n")
      
      kw_result <- kruskal.test(
        as.formula(paste0("`", region, "` ~ GROUP")),
        data = data
      )
      
      overall_p <- kw_result$p.value
      test_used <- "Kruskal-Wallis"
      
      if (overall_p < 0.05) {
        
        dunn_result <- dunn.test(
          data[[region]],
          data$GROUP,
          method = "holm"
        )
        
        dunn_df <- data.frame(
          comparisons = dunn_result$comparisons,
          p_values = dunn_result$P.adjusted
        )
        
        significant_comparisons <- dunn_df[
          dunn_df$p_values < 0.05, ]
      }
      
    } else {
      
      cat("Normal distribution → ANOVA\n")
      
      anova_result <- aov(
        as.formula(paste0("`", region, "` ~ GROUP")),
        data = data
      )
      
      overall_p <- summary(anova_result)[[1]][["Pr(>F)"]][1]
      test_used <- "ANOVA"
      
      if (overall_p < 0.05) {
        
        tukey_result <- TukeyHSD(anova_result, "GROUP")
        
        tukey_df <- data.frame(
          comparisons = rownames(tukey_result$GROUP),
          p_values = tukey_result$GROUP[, "p adj"]
        )
        
        significant_comparisons <- tukey_df[
          tukey_df$p_values < 0.05, ]
      }
    }
    
    
    ####################################
    # Save results
    ####################################
    
    sig_comparisons_text <- ""
    
    if (nrow(significant_comparisons) > 0) {
      sig_comparisons_text <- paste(
        apply(significant_comparisons, 1, function(row) {
          paste0(
            row["comparisons"],
            " (",
            format_pvalue(as.numeric(row["p_values"])),
            ")"
          )
        }),
        collapse = "; "
      )
    }
    
    result_row <- data.frame(
      network = network_name,
      region = region,
      test_used = test_used,
      overall_p = overall_p,
      significant = overall_p < 0.05,
      significant_comparisons = sig_comparisons_text
    )
    
    network_results_df <- rbind(network_results_df, result_row)
    
    
    ####################################
    # Plot significant regions only
    ####################################
    
    if (overall_p < 0.05) {
      
      significant_regions_count <- significant_regions_count + 1
      
      p <- ggplot(
        data,
        aes(x = GROUP, y = .data[[region]], fill = GROUP)
      ) +
        geom_boxplot(alpha = 0.8, outlier.shape = NA) +
        geom_jitter(width = 0.2, alpha = 0.7) +
        scale_fill_manual(values = colors) +
        theme_bw() +
        labs(
          title = paste("Volume in", region),
          subtitle = paste("Network:", network_name),
          x = "Groups",
          y = "Volume"
        )
      
      filename <- file.path(
        network_output_dir,
        paste0(
          "Boxplot_Volume_",
          gsub("[^a-zA-Z0-9]", "_", region),
          ".png"
        )
      )
      
      ggsave(filename, plot = p,
             width = 8, height = 6, dpi = 300)
    }
  }
  
  
  ##########################################################
  # Save network results
  ##########################################################
  
  all_results_df <- rbind(all_results_df, network_results_df)
  
  write.csv(
    network_results_df,
    file.path(
      network_output_dir,
      paste0(network_name, "_volume_analysis_results.csv")
    ),
    row.names = FALSE
  )
  
  significant_results <- network_results_df %>%
    filter(significant == TRUE)
  
  write.csv(
    significant_results,
    file.path(
      network_output_dir,
      paste0(network_name, "_significant_volume_results.csv")
    ),
    row.names = FALSE
  )
}


############################################################
# GLOBAL RESULTS
############################################################

write.csv(
  all_results_df,
  file.path(
    output_base_dir,
    "all_networks_volume_analysis_results.csv"
  ),
  row.names = FALSE
)

all_significant_results <- all_results_df %>%
  filter(significant == TRUE)

write.csv(
  all_significant_results,
  file.path(
    output_base_dir,
    "all_networks_significant_volume_results.csv"
  ),
  row.names = FALSE
)

cat("\n=== ANALYSIS COMPLETED ===\n")