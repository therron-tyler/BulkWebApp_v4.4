# Gene Visualization: Winter Lab Database of bulk RNA-seq Data 1/22/2025
# Tyler Therron, MS
# Deborah Winter, PhD
# Version 4.3.1.9

library(shiny)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(readr)
library(shiny.fluent)
library(htmlwidgets) 
library(webshot) 
library(future)
library(DT)
library(future.apply)
library(ggdendro)
library(stats)
library(RColorBrewer)
library(logger)
library(ragg)

# Function to perform k-means clustering
k_means_clustering <- function(expression_data, num_clusters) {
  # Scale the data before clustering
  
  if (nrow(expression_data) == 0 || ncol(expression_data) == 0) {
    stop("Error: The input expression data is empty. Please provide valid gene expression data.")
  }
  
  # log_info("Number of clusters for expression data: {num_clusters}")
  # log_info("expression_data in k_means_clustering FXN: {paste(capture.output(print(expression_data)), collapse = '\n')}")
  
  # Convert all columns to numeric
  expression_data[] <- lapply(expression_data, function(x) as.numeric(trimws(as.character(x))))
  
  normalized_data_t <- apply(expression_data, 1, min_max_normalization, min_value = -1, max_value = 1)
  
  if (nrow(normalized_data_t) == 0 || ncol(normalized_data_t) == 0) {
    stop("Error: The normalization process returned an empty dataset.")
  }
  
  # Convert the matrix back to a data frame with the same column names
  normalized_data <- t(as.data.frame(normalized_data_t))
  
  # log_info("expression_data in normalized_data dataframe: {paste(capture.output(print(normalized_data)), collapse = '\n')}")
  
  num_unique_data_points <- nrow(unique(normalized_data))
  
  # log_info("Number of unique data points: {num_unique_data_points}")
  
  # Identify valid rows (i.e., rows without NA, NaN, or Inf)
  valid_rows <- apply(normalized_data, 1, function(row) !any(is.na(row) | is.nan(row) | is.infinite(row)))
  
  # Subset the matrix to keep only valid rows
  normalized_data <- normalized_data[valid_rows, , drop = FALSE]
  
  # Log removed genes if any
  removed_genes <- setdiff(rownames(normalized_data), rownames(normalized_data))
  if (length(removed_genes) > 0) {
    # log_info("Removed genes due to NA/NaN/Inf values: {paste(removed_genes, collapse = ', ')}")
  }
  
  # Ensure at least one gene remains before calling kmeans_heatmap
  if (nrow(normalized_data) == 0) {
    # log_error("All genes were removed due to NA/NaN/Inf values. Cannot generate heatmap.")
    return(list(error = TRUE, message = "All genes contained invalid values and were removed. Heatmap cannot be generated."))
  }
  
  if (num_unique_data_points < num_clusters) {
    stop(paste("Error: There are more cluster centers (", num_clusters, 
               ") than distinct data points (", num_unique_data_points, 
               "). Please reduce the number of clusters or broaden your filtering criteria."))
  }
  
  # Perform k-means clustering
  set.seed(40) # For reproducibility
  
  kmeans_result <- kmeans(normalized_data, centers = num_clusters, nstart = 25)
  
  Clusters <- data.frame(Cluster = kmeans_result$cluster)
  
  norm_data_and_clusters <- cbind(normalized_data, Clusters)
  
  # log_info("expression_data in norm_data_and_clusters dataframe: {paste(capture.output(print(norm_data_and_clusters)), collapse = '\n')}")
  return(norm_data_and_clusters)
}
