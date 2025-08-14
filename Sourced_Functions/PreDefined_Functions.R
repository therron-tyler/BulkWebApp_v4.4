#  ======================================. Libraries ==================
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
library(shiny)
library(shiny.fluent)
library(htmlwidgets) 
library(webshot) 
library(future)
library(DT)
library(future.apply)
library(data.table)

print(getOption("future.globals.maxSize"))  # Should return 8589934592 (8 GiB)

# UI list checker for data so the CPM and FPKM option can be updated accordingly
datasets_info %<-% function(fpkm_path, cpm_path) {
  available_datasets <- list.files(fpkm_path, full.names = FALSE, recursive = FALSE)
  all_datasets <- list.files(cpm_path, full.names = FALSE, recursive = FALSE)
  
  unavailable_datasets <- setdiff(all_datasets, available_datasets)
  
  list(available = available_datasets, unavailable = unavailable_datasets)
}



format_and_merge_UI_chosen_genes <- function(gene_expression_multiple_genes, groupfile_in) {
  df_in <- gene_expression_multiple_genes
  
  df_in_Count_column <- future_lapply(df_in, function(x) {
    x <- as.data.frame(t(data.frame(x)))
    colnames(x) <- rep("Count", ncol(x))
    return(x)
  })
  
  full_group_in_list <- future_lapply(df_in_Count_column, function(x) {
    full_group_in <- data.frame(Group = NA, Sample = row.names(x))
    row.names(full_group_in) <- row.names(x)
    return(full_group_in)
  })
  
  ## ---------------------------------------- Group File Processing 
  group_in <- groupfile_in
  
  group_in <- as.data.frame(t(group_in)) %>% 
    dplyr::rename(group_color = V2)
  
  groupfile_order <- t(unique(t(group_in))) #10/24/23
  
  # make column color vector an object for later
  color_vector_column <- group_in$group_color # new 10/9
  
  # remove color vector column
  group_in <- group_in %>% dplyr::select(-group_color) # new 10/9 - can move this up to line 86
  
  colnames(group_in) <- c("Group")
  group_in$Sample <- row.names(group_in)
  
  # Keep only rows that exist in group_in for both df_in_Count_column and full_group_in_list
  matching_samples <- group_in$Sample
  
  df_in_Count_column <- future_lapply(df_in_Count_column, function(x) {
    x <- x[matching_samples, , drop = FALSE]  # Filter to only keep rows with matching samples
    return(x)
  })
  
  # ---------------------------------------- Group File Processing 
  
  full_group_in_list_joined <- future_lapply(full_group_in_list, function(x) {
    x[row.names(group_in), ] <- group_in
  })
  
  group_and_data_combine <- future_mapply(function(counts, groupfile_info) {
    cbind(counts, groupfile_info)
  }, df_in_Count_column, full_group_in_list_joined, SIMPLIFY = FALSE)
  
  # - wrap in lapply
  group_and_data_list <- future_lapply(group_and_data_combine, function(x) {
    group_and_data <- x[, !duplicated(names(x))] %>% 
      filter(!is.na(Group)) %>% 
      arrange(match(Group, groupfile_order))
    group_and_data$color_column <- color_vector_column
    return(group_and_data)
  })
  
  return(group_and_data_list)
}

# 20231204 - this function was refactored to only contain 1 for loop instead of 2. 
add_group_traces %<-% function(p, data, name, ordered_groups, color_mapping) {
  name_pieces = strsplit(name, " - ")
  print("ADD_GROUP_TRACES name")
  pieces = name_pieces[[1]]
  citation = pieces[length(pieces)]
  
  for (grp in ordered_groups) {
    grp_data <- subset(data, Group == grp)
    # grp_hover_text <- paste("Replicate: ", grp_data$Sample, "\nCitation: ", citation)
    grp_hover_text = paste("Citation:", citation, "<br>","Replicate:", grp_data$Sample)
    
    # citation_hover
    # 
    # Add box plot trace
    p <- add_trace(p, data = grp_data, x = ~Group, y = ~Count, type = "box",
                   color = I(color_mapping[grp]), name = grp,
                   text = grp_hover_text, hoverinfo = 'text+y')
    
    # Add scatter plot points
    grp_colors <- rep(color_mapping[grp], times = nrow(grp_data))
    
    p <- add_trace(p, data = grp_data, x = ~Group, y = ~Count,
                   type = 'scatter', mode = 'markers',
                   marker = list(color = I(grp_colors), size = 10, opacity = 0.6),
                   text = grp_hover_text, hoverinfo = 'text+y',
                   showlegend = FALSE)
  }
  return(p)
}



# 12/8/23 - refactored plotting code - also changing the titles in this one
plot_the_data <- function(data, name, cpm_fpkm_label, gene_name) {
  print("PLOT_THE_DATA gene_name")
  print(gene_name)
  print("PLOT_THE_DATA name")
  print(name)
  
  data$Group <- as.factor(data$Group) # convert to factor for plotting purposes
  ordered_groups <- unique(data$Group)
  upper_limit <- max(data$Count) + 20
  
  color_mapping <- setNames(object = unique(data$color_column), nm = unique(data$Group))
  hover_text <- paste("Replicate: ", data$Sample) # assuming 'SampleName' is the column with individual sample names
  
  gene_name <- paste(toupper(substring(gene_name, 1, 1)), 
                     substring(gene_name, 2), 
                     sep = "")
  
  # Initialize an empty plotly object
  p <- plot_ly()
  
  p <- add_group_traces(p, data, name, ordered_groups, color_mapping) %>%
    layout(
      title = list(text = paste0("<span style='color:blue;'>", gene_name, "</span>"), 
                   font = list(size = 20),
                   x = 0,  # Align title to the left
                   xref = 'paper' ),
      # annotations = list(
      #   list(
      #     text = name,
      #     font = list(size = 15, family = "Arial Italic"),
      #     x = 0.5,
      #     y = -0.5,
      #     xref = 'paper',
      #     yref = 'paper',
      #     showarrow = FALSE,
      #     xanchor = 'center',
      #     yanchor = 'bottom'
      #   )
      # ),
      xaxis = list(title = "", categoryorder = "array", categoryarray = ordered_groups),
      yaxis = list(title = cpm_fpkm_label,
                   zeroline = TRUE,
                   zerolinewidth = 2,
                   zerolinecolor = "#000000",
                   range = c(0, upper_limit)),
      margin = list(t = 100, b = 300, l = 50) # Adjust bottom margin to give space for subtitle
    )
  
  return(p)
}

# df_pri - calls this function
process_gene_data <- function(data_in, gene_of_interest, selNum_function) {
  data_in %>%
    filter(tolower(data_in[, selNum_function()]) == tolower(as.character(gene_of_interest))) %>%
    distinct(.keep_all = TRUE) %>%
    dplyr::select(3:ncol(data_in))
}



# Define a function to get the title based on user's choice
get_data_title %<-% function(selected_data, uploaded_data_title) {
  if (selected_data == 'Upload your Own Dataset.') {
    return(uploaded_data_title)  # If user uploaded data, return the input title
  } else {
    return(selected_data)  # Otherwise, return the selected data
  }
}

split_name_at_first_dash <- function(name) {
  # Check if there's a dash in the name
  if (grepl(" - ", name, fixed = TRUE)) {
    # Split the string at the first dash, and limit to split only once
    parts <- strsplit(name, " - ", fixed = TRUE)[[1]]
    # If the split resulted in more than one part, concatenate them with newline characters
    if (length(parts) > 1) {
      new_name <- paste(parts, collapse = "\n")
    } else {
      # If there's only one part, no dash was found, return the original name
      new_name <- name
    }
  } else {
    # If there is no dash, just return the original name
    new_name <- name
  }
  return(new_name)
}

# sort datasets by numeric rather than lexigraphical order - week of halloween
sort_numeric %<-% function(files) {
  files[order(as.numeric(gsub("^(\\d+).*", "\\1", files)))]
}


# Dataset Info read in code - early november
dataset_info %<-% readr::read_csv("./Sourced_CSVs/Dataset_Information_BulkGEapp.csv")

# Define the base URL for GEO IDs -20231117
base_url <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="

# Publication URL -20231117
PM_id_base_url <- "https://pubmed.ncbi.nlm.nih.gov/"

geo_link_creator <- function(id) {
  if(id != "There is no GEO accession number.") {
    # Create the hyperlink HTML tag
    return(HTML(paste0("GEO: ",'<a href="', base_url, id, '" target="_blank">', id, '</a>')))
  } else {
    return(id) # Leave "N/A" as is
  }
}

# Modify the GEO IDs to include hyperlinks-20231117
dataset_info$GEO <- future_sapply(dataset_info$GEO, function(id) {
  if(id != "N/A") {
    # Create the hyperlink HTML tag
    paste0('<a href="', base_url, id, '" target="_blank">', id, '</a>')
  } else {
    id # Leave "N/A" as is
  }
}, USE.NAMES = FALSE)

# Modify the PMC IDs to include hyperlinks-20231117
dataset_info$Publication <- future_sapply(dataset_info$Publication, function(id) {
  if(!(id %in% c("In preparation", "Unanalyzed","bioRxiv"))) {
    # Create the hyperlink HTML tag
    paste0("PM ID: ",'<a href="', PM_id_base_url, id, '" target="_blank">', id, '</a>')
  } else {
    id # Leave "N/A" as is
  }
}, USE.NAMES = FALSE)

# Legend table read in
legend_info %<-% readr::read_csv("./Sourced_CSVs/Legend_Table_Bulk_WebApp.csv")

## ----------------------------- initial function in the heatmap viz function

# original
merge_gene_data <- function(gene_data_list) {
  # Step 1: Rename the `Count` column in each dataframe to its gene name
  gene_data_list <- lapply(names(gene_data_list), function(gene) {
    df <- gene_data_list[[gene]]
    colnames(df)[colnames(df) == "Count"] <- gene  # Rename the Count column to the gene name
    return(df)
  })
  
  # Step 2: Merge all dataframes in the list into a single dataframe
  # We'll use `Reduce()` with `merge()` to recursively combine them
  merged_df <- Reduce(function(x, y) {
    merge(x, y, by = c("Group", "Sample", "color_column"), all = TRUE)
  }, gene_data_list)
  
  return(merged_df)
}

## ----------------------------- initial function in the heatmap viz function

## ----------------------------- min-max norm fxn in the heatmap viz function

min_max_normalization <- function(x, min_value = -1, max_value = 1) {
  min_x <- min(x)
  max_x <- max(x)
  normalized_x <- (x - min_x) / (max_x - min_x) * (max_value - min_value) + min_value
  return(normalized_x)
}

## ----------------------------- min-max norm fxn in the heatmap viz function

# parallelized Function to assign names to the sub-list items
assign_names_to_sublist <- function(data_list, names_list) {
  future_lapply(data_list, function(sub_list) {
    names(sub_list) <- names_list
    return(sub_list)
  })
}

## ----------------------------- load datasets function

# 20241216 - Gene choices to speed up laoding times
gene_choices_for_plots %<-% readr::read_csv("./Sourced_CSVs/GEV_app_GeneChoices.csv")

transpose_groupfile <- function(groupfile) {
  group_t <- data.frame(t(groupfile), row.names = NULL)
  colnames(group_t) <- group_t[1,]
  group_t2 <- group_t[c(2,3), ]
  row.names(group_t2) <- c("V1", "V2")
  return(group_t2)
}

load_datasets <- function(datasets, normalization_mode, selected_gene_id, selNum) {
  # selNum = 1 means filter by Symbol (Ensembl ID)
  # selNum = 2 means filter by Gene_Symbol
  
  # Determine which column to filter on:
  filter_column <- if (selNum == 1) "Symbol" else "Gene_Symbol"
  
  data_list <- list()
  group_list <- list()
  
  for (dataset in datasets) {
    data_path <- if (normalization_mode == "cpm") {
      paste0('./Primary/', dataset, '/Data.csv')
    } else {
      paste0('./FPKM_expression_files/', dataset, '/Data.csv')
    }
    
    # Read just the header first
    header <- names(fread(data_path, nrows = 0))
    
    # Set column classes:
    # First two columns are character (Symbol, Gene_Symbol), rest are numeric.
    col_classes <- c("character", "character", rep("numeric", length(header) - 2))
    
    full_data <- fread(data_path, colClasses = col_classes)
    
    # Filter based on selected column:
    single_gene_data <- full_data[get(filter_column) %in% selected_gene_id]
    single_gene_data <- as.data.frame(single_gene_data)
    
    data_list[[dataset]] <- single_gene_data
    
    group_path <- if (normalization_mode == "cpm") {
      paste0('./Primary/', dataset, '/GroupFile.csv')
    } else {
      paste0('./FPKM_expression_files/', dataset, '/GroupFile.csv')
    }
    
    group_data <- read.csv(group_path)
    
    data_list[[dataset]] <- as.data.frame(data_list[[dataset]])
    group_list[[dataset]] <- as.data.frame(group_data)
  }
  
  return(list(data = data_list, groups = group_list))
}
## ----------------------------- load datasets function

## ------ check for a match between user uploaded GE table and group file Sample column 

validate_uploaded_data <- function(expr_data, group_file) {
  # Extract sample column from the group file
  group_samples <- group_file$Sample
  
  # Extract column names from the expression data (excluding metadata columns)
  expr_samples <- colnames(expr_data)[!colnames(expr_data) %in% c("Gene_Symbol", "Symbol")]
  
  # Check if all group samples exist in expression data columns
  missing_samples <- setdiff(group_samples, expr_samples)
  
  # Return results
  if (length(missing_samples) > 0) {
    return(list(
      valid = FALSE,
      message = paste("The following samples in the group file are missing from the expression file:", 
                      paste(missing_samples, collapse = ", "))
    ))
  }
  
  return(list(valid = TRUE, message = "Group file and expression file match successfully."))
}

## process only one dataset at a time functions ------

# Process added datasets
process_added_datasets <- function(datasets) {
  # Load and process added datasets
  results_list <- preprocessing_selected_genes_datasets()[datasets]
  groups <- processed_datasets$groups[datasets]
  
  group_and_data <- future_mapply(function(dataset_genes_list, groupfile_list) {
    format_and_merge_UI_chosen_genes(dataset_genes_list, groupfile_list)
  }, results_list, groups, SIMPLIFY = FALSE)
  
  # Update the formatted_data reactive
  formatted_data$data <- c(formatted_data$data, group_and_data)
}

# Process removed datasets
process_removed_datasets <- function(datasets) {
  # Remove datasets from formatted_data reactive
  formatted_data$data <- formatted_data$data[!(names(formatted_data$data) %in% datasets)]
}

