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

source("./Sourced_Functions/HM_PreDefined_Functions.R")

prepare_heatmap_data <- function(data, grps, Group_or_Sample, exclude_cols = NULL) {
  # Remove specified columns if needed
  if (!is.null(exclude_cols)) {
    data <- data %>% select(-one_of(exclude_cols))
  }
  
  # log_info("data in prepare_heatmap_data: {paste(capture.output(print(data)), collapse = '\n')}")
  # log_info("grps in prepare_heatmap_data: {paste(capture.output(print(grps)), collapse = '\n')}")
  
  # Determine col_ID based on Group_or_Sample
  if (Group_or_Sample == "Group") {
    col_ID <- unique(grps$Group)
  } else {
    col_ID <- grps$Sample
  }
  
  # Order columns
  column_order <- match(col_ID, colnames(data))
  data_ordered <- data[, column_order, drop = FALSE]
  
  # Convert columns to numeric
  data_ordered[] <- lapply(data_ordered, function(x) as.numeric(trimws(as.character(x))))
  
  return(data_ordered)
}

generate_heatmap_annotations <- function(grps, Group_or_Sample) {
  if (Group_or_Sample == "Group") {
    column_annotation_input <- unique(grps$Group)
  } else {
    column_annotation_input <- grps$Group
  }
  
  annotations <- generate_annotations(column_annotation_input, viridis::turbo)
  EXP_groups <- annotations$exp_groups
  col_annotations <- annotations$col_annotations_exp
  
  return(list(EXP_groups = EXP_groups, col_annotations = col_annotations))
}

# II
calculate_heatmap_parameters <- function(num_genes, num_columns) {
  # Dynamically adjust spacing strategy
  dynamic_space <- if (num_columns <= 30) {
    generate_space_text(num_columns, base_space = 12)
  } else {
    remove_space_text(num_columns, base_space = 8)
  }
  
  # Set adaptive plot height
  plot_height <- 800 + (min(40 * num_genes, 3000))  # Cap max height to prevent excessive scrolling
  
  # Smart font scaling: prevents tiny unreadable text
  font_size <- case_when(
    num_genes > 100 ~ 6,
    num_genes > 50 ~ 8,
    num_genes > 25 ~ 10,
    num_genes > 10 ~ 12,
    num_genes > 5  ~ 14,
    TRUE ~ 16
  )
  
  return(list(dynamic_space = dynamic_space, plot_height = plot_height, font_size = font_size))
}

# Helper function to calculate dynamic column annotation size
calculate_annotation_size <- function(column_names) {
  num_columns <- length(column_names)  # Get number of columns
  
  # Adaptive font scaling strategy
  font_size <- case_when(
    num_columns > 50 ~ 8,       # Many columns, small font
    num_columns > 45 ~ 10,      # Medium number of columns
    num_columns > 40 ~ 11,       # Many columns, small font
    num_columns > 35 ~ 12,
    num_columns > 30 ~ 13,       # Many columns, small font
    num_columns > 25 ~ 14,      # Medium number of columns
    num_columns > 20 ~ 15,       # Many columns, small font
    num_columns > 15 ~ 16,
    num_columns > 10 ~ 17,
    TRUE ~ 18                   # Default readable size
  )
  
  return(font_size)
}

# Helper function to calculate dynamic font size
calculate_legend_font_size <- function(column_names) {
  max_name_length <- max(nchar(column_names))  # Find the longest column name
  
  print("character max length")
  print(max_name_length)
  # Adaptive font scaling strategy
  font_size <- case_when(
    max_name_length > 50 ~ 7,       # Many columns, small font
    max_name_length > 45 ~ 8,      # Medium number of columns
    max_name_length > 40 ~ 9,       # Many columns, small font
    max_name_length > 35 ~ 10,
    max_name_length > 30 ~ 11,       # Many columns, small font
    max_name_length > 25 ~ 12,      # Medium number of columns
    max_name_length > 20 ~ 13,       # Many columns, small font
    max_name_length > 15 ~ 14,
    max_name_length > 10 ~ 15,
    TRUE ~ 16                   # Default readable size
  )
  
  return(font_size)
}

# helper function that calculates the width that the heatmap legend should be based on the font size
calculate_legend_width <- function(column_names) {
  max_name_length <- max(nchar(column_names))  # Find the longest column name
  
  # Adaptive font scaling strategy
  font_size <- case_when(
    max_name_length > 50 ~ 230,       # Many columns, small font
    max_name_length > 45 ~ 220,      # Medium number of columns
    max_name_length > 40 ~ 210,       # Many columns, small font
    max_name_length > 35 ~ 200,
    max_name_length > 30 ~ 180,       # Many columns, small font
    max_name_length > 25 ~ 170,      # Medium number of columns
    max_name_length > 20 ~ 170,       # Many columns, small font
    max_name_length > 15 ~ 160,
    max_name_length > 10 ~ 150,
    TRUE ~ 140                   # Default readable size
  )
  
  return(font_size)
}

create_heatmap_plot <- function(normalized_data, col_annotations, dynamic_space, plot_height) {
  
  # Apply dynamic font size
  annotation_font_size <- calculate_legend_font_size(colnames(normalized_data))
  
  heatmap_plot <- plot_ly(
    x = colnames(normalized_data),
    y = rownames(normalized_data),
    z = as.matrix(normalized_data),
    width = "2000",
    height = plot_height,
    type = "heatmap",
    colors = colorRamp(c("blue", "white", "red")),
    showscale = TRUE,
    colorbar = list(
      title = list(
        text = "Min-Max Normalized Gene Expression",
        font = list(size = 14, family = "Trebuchet MS", color = "grey10"),
        side = "right"
      ),
      x = 1.065,  # Moves the color bar to the right
      tickfont = list(size = 12),
      len = 0.25,
      y = 0.65,
      lenmode = "fraction",
      titleside = "right",
      title_standoff = 45
    ),
    hovertemplate = paste(
      "Sample: %{x}<br>",
      "Gene: %{y}<br>",
      "Normalized Expression: %{z:.2f}<extra></extra>"
    )
  ) %>%
    add_annotations(
      x = colnames(normalized_data),
      y = 1.001,
      text = dynamic_space,
      showarrow = FALSE,
      font = list(size = calculate_annotation_size, color = as.character(col_annotations)),
      bgcolor = col_annotations,
      xref = "x",
      yref = "paper",
      xanchor = "center",
      yanchor = "bottom",
      align = "center"
    )
  
  return(heatmap_plot)
}

#VI
calculate_y_positions <- function(num_items, y_start = 0.98, fixed_gap = 0.025, y_split = 0.705) {
  if (num_items == 1) {
    print("y positions")
    print(y_start)
    return(y_start)
  }
  
  print("y positions (before adjustment):")
  
  # Compute initial y-positions using fixed gap
  y_positions <- y_start - (seq_len(num_items) - 1) * fixed_gap
  print(y_positions)
  return(y_positions)
}

#III
calculate_legend_x_coords <- function(column_names, y_positions, x_base = 1.2, x_shift = 0.165, y_split = 0.705) {
  # Ensure column_names is not empty
  if (length(column_names) == 0) {
    stop("Error: column_names is empty. Cannot calculate legend positions.")
  }
  
  # Ensure y_positions and column_names have the same length
  if (length(column_names) != length(y_positions)) {
    stop("Error: Mismatch between column_names and y_positions.")
  }
  
  print("Original y_positions:")
  print(y_positions)
  
  # Initialize all x_positions at x_base
  x_positions <- rep(x_base, length(y_positions))
  
  # Identify indices that belong to the second column
  second_col_indices <- which(y_positions <= y_split)
  
  if (length(second_col_indices) > 0) {
    print("Shifting second column labels to the right...")
    
    # Assign second column a higher x position
    x_positions[second_col_indices] <- x_base + x_shift
    
    # Adjust y_positions for the second column to maintain even spacing
    y_positions[second_col_indices] <- max(y_positions) - 
      (seq_along(second_col_indices) - 1) * 0.025
  }
  
  print("Final x_positions:")
  print(x_positions)
  
  print("Updated y_positions:")
  print(y_positions)
  
  return(list(x_positions = x_positions, y_positions = y_positions))
}

# II
# add_legends_to_heatmap <- function(plot, EXP_groups, col_annotations) {
#   unique_groups <- unique(EXP_groups)
#   unique_colors <- unique(col_annotations)
#   num_items <- length(unique_groups)
# 
#   # Compute y positions dynamically
#   print("num_items")
#   log_info(" num_items: {paste(capture.output(print(unique(num_items))), collapse = '\n')}")
#   log_info(" unique_groups: {paste(capture.output(print(unique(unique_groups))), collapse = '\n')}")
#   y_positions <- calculate_y_positions(num_items)
# 
#   # Compute x positions dynamically
#   x_and_y_positions <- calculate_legend_x_coords(unique_groups, y_positions)
#   
#   x <- x_and_y_positions$x_positions
#   y <- x_and_y_positions$y_positions
# 
#   # Adaptive font size & legend width
#   annotation_font_size <- calculate_legend_font_size(unique_groups)
#   legend_width <- calculate_legend_width(unique_groups)
# 
#   # Generate list of annotations
#   exp_annotations <- lapply(seq_along(unique_groups), function(i) {
#     list(
#       x = x[i],  # ✅ Correct X position
#       y = y[i],  # ✅ Correct Y position
#       xref = "paper",
#       yref = "paper",
#       text = unique_groups[i],
#       showarrow = FALSE,
#       font = list(color = "#000000", size = annotation_font_size),
#       bgcolor = unique_colors[i],
#       borderpad = 4,
#       bordercolor = "#FFFFFF",
#       borderwidth = 1,
#       align = "center",
#       width = legend_width,
#       height = 20
#     )
#   })
# 
#   # Apply legend annotations to the plot
#   plot <- plot %>%
#     layout(
#       xaxis = list(tickangle = -45, tickfont = list(size = 13)),
#       margin = list(l = 100, r = 600, b = 50, t = 80),
#       xaxis = list(scaleanchor = "x", scaleratio = 1),
#       yaxis = list(scaleanchor = "y", scaleratio = 1),
#       legend = list(
#         itemsizing = "constant",
#         traceorder = "normal",
#         font = list(family = "sans-serif", size = 12, color = "#000"),
#         bgcolor = "#E2E2E2",
#         bordercolor = "#FFFFFF",
#         borderwidth = 2,
#         orientation = "v"
#       ),
#       annotations = exp_annotations
#     )
# 
#   return(plot)
# }

# III
add_legends_to_heatmap <- function(plot, EXP_groups, col_annotations) {
  unique_groups <- unique(EXP_groups)
  unique_colors <- unique(col_annotations)
  num_items <- length(unique_groups)
  
  # Compute y positions dynamically
  y_positions <- calculate_y_positions(num_items)
  
  # Compute x positions dynamically
  x_and_y_positions <- calculate_legend_x_coords(unique_groups, y_positions)
  
  x <- x_and_y_positions$x_positions
  y <- x_and_y_positions$y_positions
  
  # Adaptive font size & legend width
  annotation_font_size <- calculate_legend_font_size(unique_groups)
  legend_width <- calculate_legend_width(unique_groups)
  
  # 🔹 **Insert a Title for the Experimental Groups**
  legend_title <- list(
    x = 1.185,  # Centered above the experimental groups
    y = max(y) + 0.02,  # Slightly above the topmost group
    text = "<b>Experimental Groups</b>",  # Bold text
    showarrow = FALSE,
    font = list(size = 14, color = "#000000", family = "Trebuchet MS"),
    xref = "paper",
    yref = "paper",
    align = "left"
  )
  
  # Generate list of annotations
  exp_annotations <- lapply(seq_along(unique_groups), function(i) {
    list(
      x = x[i],  # ✅ Correct X position
      y = y[i],  # ✅ Correct Y position
      xref = "paper",
      yref = "paper",
      text = unique_groups[i],
      showarrow = FALSE,
      font = list(color = "#000000", size = annotation_font_size),
      bgcolor = unique_colors[i],
      borderpad = 4,
      bordercolor = "#FFFFFF",
      borderwidth = 1,
      align = "center",
      width = legend_width,
      height = 20
    )
  })
  
  # Add the title annotation to the list
  exp_annotations <- c(list(legend_title), exp_annotations)
  
  # Apply legend annotations to the plot
  plot <- plot %>%
    layout(
      xaxis = list(tickangle = -45, tickfont = list(size = 13)),
      margin = list(l = 100, r = 600, b = 50, t = 80),
      xaxis = list(scaleanchor = "x", scaleratio = 1),
      yaxis = list(scaleanchor = "y", scaleratio = 1),
      legend = list(
        itemsizing = "constant",
        traceorder = "normal",
        font = list(family = "sans-serif", size = 12, color = "#000"),
        bgcolor = "#E2E2E2",
        bordercolor = "#FFFFFF",
        borderwidth = 2,
        orientation = "v"
      ),
      annotations = exp_annotations  # ✅ Includes the title now!
    )
  
  return(plot)
}

kmeans_heatmap <- function(HM_expression_data_with_cluster_info, grps, Group_or_Sample) {
  # Prepare data and remove clustering columns
  HM_expression_data_ordered <- prepare_heatmap_data(
    HM_expression_data_with_cluster_info, grps, Group_or_Sample, exclude_cols = c("Cluster", "ClusterColor")
  )
  
  # Extract cluster information
  cluster_info <- HM_expression_data_with_cluster_info %>% select(Cluster, ClusterColor)
  
  # Generate annotations
  annotations_list <- generate_heatmap_annotations(grps, Group_or_Sample)
  EXP_groups <- annotations_list$EXP_groups
  col_annotations <- annotations_list$col_annotations
  
  num_genes <- nrow(HM_expression_data_ordered)
  num_columns <- ncol(HM_expression_data_ordered)
  
  # Calculate dynamic parameters
  params <- calculate_heatmap_parameters(num_genes, num_columns)
  dynamic_space <- params$dynamic_space
  plot_height <- params$plot_height
  font_size <- params$font_size
  
  # Create the base heatmap
  p <- create_heatmap_plot(HM_expression_data_ordered, col_annotations, dynamic_space, plot_height)
  
  # Add cluster annotations
  p <- p %>%
    add_trace(
      y = seq_along(cluster_info$Cluster),
      x = rep(ncol(HM_expression_data_ordered) + 1, length(cluster_info$Cluster)),
      z = matrix(cluster_info$ClusterColor, nrow = length(cluster_info$Cluster), ncol = 1),
      type = "heatmap",
      colorscale = "Greys",
      showscale = FALSE
    ) %>%
    layout(
      annotations = lapply(seq_along(cluster_info$Cluster), function(i) {
        list(
          x = 1.05,
          y = i - 1,
          text = paste0("Cluster ", cluster_info$Cluster[i]),
          xref = "paper",
          yref = "y",
          showarrow = FALSE,
          font = list(size = font_size),
          bgcolor = cluster_info$ClusterColor[i]
        )
      })
    )
  
  # Add legends
  p <- add_legends_to_heatmap(p, EXP_groups, col_annotations)
  
  return(p)
}

HierClust_HM <- function(clean_variable_genes, grps, Group_or_Sample) {
  # Prepare data
  clean_variable_genes_ordered <- prepare_heatmap_data(clean_variable_genes, grps, Group_or_Sample)
  
  # Perform hierarchical clustering
  hclust_rows <- hierarchical_clustering(clean_variable_genes_ordered)
  ordered_GE_data <- clean_variable_genes_ordered[hclust_rows$order, ]
  
  # Min-max normalization
  normalized_data_t <- apply(ordered_GE_data, 1, min_max_normalization, min_value = -1, max_value = 1)
  normalized_data <- t(as.data.frame(normalized_data_t))
  
  # Generate annotations
  annotations_list <- generate_heatmap_annotations(grps, Group_or_Sample)
  EXP_groups <- annotations_list$EXP_groups
  col_annotations <- annotations_list$col_annotations
  
  num_genes <- nrow(normalized_data)
  num_columns <- ncol(normalized_data)
  
  # Calculate dynamic parameters
  params <- calculate_heatmap_parameters(num_genes, num_columns)
  dynamic_space <- params$dynamic_space
  plot_height <- params$plot_height
  font_size <- params$font_size
  
  # Create the base heatmap
  heatmap_plot <- create_heatmap_plot(normalized_data, col_annotations, dynamic_space, plot_height)
  
  # Add legends
  heatmap_plot <- add_legends_to_heatmap(heatmap_plot, EXP_groups, col_annotations)
  
  return(heatmap_plot)
}

no_clustering_heatmap <- function(clean_variable_genes, grps, Group_or_Sample) {
  # Prepare data
  clean_variable_genes_ordered <- prepare_heatmap_data(clean_variable_genes, grps, Group_or_Sample)
  
  # Min-max normalization
  normalized_data_t <- apply(clean_variable_genes_ordered, 1, min_max_normalization, min_value = -1, max_value = 1)
  normalized_data <- t(as.data.frame(normalized_data_t))
  
  # Generate annotations
  annotations_list <- generate_heatmap_annotations(grps, Group_or_Sample)
  EXP_groups <- annotations_list$EXP_groups
  col_annotations <- annotations_list$col_annotations
  
  num_genes <- nrow(normalized_data)
  num_columns <- ncol(normalized_data)
  
  # Calculate dynamic parameters
  params <- calculate_heatmap_parameters(num_genes, num_columns)
  dynamic_space <- params$dynamic_space
  plot_height <- params$plot_height
  font_size <- params$font_size
  
  # Create the base heatmap
  heatmap_plot <- create_heatmap_plot(normalized_data, col_annotations, dynamic_space, plot_height)
  
  # Add legends
  heatmap_plot <- add_legends_to_heatmap(heatmap_plot, EXP_groups, col_annotations)
  
  return(heatmap_plot)
}