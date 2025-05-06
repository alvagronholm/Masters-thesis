# Code for composition analysis of global and immune cell types

# Author: Alva Grönholm
# Creation date: 15.8.2024
# Modification date: 3.5.2025

# Purpose
  # This code generates heatmaps and barplots to visualize the composition of global and immune cell types across patients. 
  # Heatmaps are annotated with progression-free survival status. 

#_______________________________________________________________________________
#Load libraries
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(patchwork)
library(tidyr)
library(tidyverse)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)

#Download data
df <- read.csv("L:/Projects/4_CellCycle/Aleksandras_pipeline_running/full_data_with_clinical_info_13122024.csv")

# Changing the names in the phenotype column
df <- df %>%
  mutate(phenotype = recode(phenotype,
                            "CD163_neg_Macrophages" = "CD163- Macrophages",
                            "CD163_pos_Macrophages" = "CD163+ Macrophages",
                            "CD4_T_cells" = "CD4+ T cells",
                            'CD4_FOXP3_T_cells' = 'FOXP3+ CD4+ T cells',
                            'CD8_T_cells' = 'CD8+ T cells',
                            'B_cells' = 'B cells',
                            'CD163_pos_DC' = 'CD163+ DC'))

# Own df for PFS status
PFS_data <- df %>%
  select(imageid, PFS_33p) %>%  
  distinct(imageid, .keep_all = TRUE)

#_______________________________________________________________________________

# Color scales

globalcelltype_colors <- c(
  "Tumor" = "#B1CCE3", 
  "Stromal" = "#F7F9A3",  
  "Immune" = "#DA3683"
)

all_celltype_colors <- c(
  "Tumor" = "#B1CCE3",       
  "CD4+ T cells" = "#DA3683",
  "FOXP3+ CD4+ T cells" = "#FFC2D7",
  "CD8+ T cells" = "#EAA09A",
  "B cells" = "#045275",              
  "CD163+ DC" = "#fedc80",             
  "DC" = "#BADE8C",                    
  "CD163+ Macrophages" = "#C6B4D5",    
  "CD163- Macrophages" = "#C060BD",    
  "Stromal" = "#F7F9A3",              
  "Endothelial" = "#EB881F"
)

# Annotation colors
PFS_color <- c("Short" = "#A4C8E1", "Middle class" = "#0076A8", "Long" = "#003D5B")

#_______________________________________________________________________________
# Cell type percentages

# Compute percentage of each phenotype
df_percentage <- df %>%
  count(phenotype) %>%
  mutate(percentage = (n / sum(n)) * 100)

# Compute percentage of each phenotype
df_percentage_global <- df %>%
  count(GlobalCellType) %>%
  mutate(percentage = (n / sum(n)) * 100)

# Count percetage of NK cells removed
phenotype_df <- read.csv("Z:/afarkkila_microscopy/Data/CellCycle/analysis/00-phenotyping/labels/phenotype_full_data.csv")
# Filter phenotype_df to keep only matching imageIDs
filtered_phenotype_df <- phenotype_df %>%
  filter(imageid %in% df$imageid)

# Compute percentage of each phenotype
df_percentage_all_pheno <- filtered_phenotype_df %>%
  count(phenotype) %>%
  mutate(percentage = (n / sum(n)) * 100)

lymphocytes <- filtered_phenotype_df[filtered_phenotype_df$phenotype %in% 
                                                      c("NK cells", "CD8+ T cells", "CD4+ T cells", 
                                                        "B cells", "CD4+ FOXP3+ T cells"), ]

df_percentage_all_lympho <- lymphocytes %>%
  count(phenotype) %>%
  mutate(percentage = (n / sum(n)) * 100)

#_______________________________________________________________________________
# Median PFS
PFS <- df %>%
  select(imageid, PFS) %>%  
  distinct(imageid, .keep_all = TRUE)
median(PFS$PFS)

#_______________________________________________________________________________
# Heatmaps for all cell types and global cell types 

# Heatmap colors
vmin=0
vmax=45 

# Heatmap for all cell types with PFS annotation
heatmap_all_cell_types <- function(df, PFS_data, PFS_color) {
  # Aggregate data by imageid and phenotype and calculate percentages
  agg_data <- df %>%
    group_by(imageid, phenotype) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(imageid) %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ungroup()
  
  # Spread data to wide format for heatmap input
  wide_data <- agg_data %>%
    select(imageid, phenotype, percentage) %>%
    pivot_wider(names_from = phenotype, values_from = percentage, values_fill = list(percentage = 0))
  
  # Merge with PFS and PFS_progression data
  merged_data <- left_join(wide_data, PFS_data, by = "imageid")
  color_mapping <- colorRamp2(c(vmin, vmax), c("white", "#D73027"))
  
  # Prepare data for heatmap
  heatmap_matrix <- as.matrix(merged_data %>% select(-imageid, -PFS_33p))
  rownames(heatmap_matrix) <- merged_data$imageid  # Set row names to imageid
  
  row_annotations <- merged_data %>% select(PFS_33p) %>% as.data.frame()
  
  
  # Create row annotation with PFS and PFS_progression
  row_annotation <- rowAnnotation(
    PFS = row_annotations$PFS_33p,
    col = list(PFS = PFS_color),
    annotation_name_rot = 55,
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
    annotation_legend_param = list(
      PFS = list(title = "PFS"))
  )
  
  # Create heatmap
  heatmap <- Heatmap(
    heatmap_matrix,
    name = "% of cell type per patient",
    col = color_mapping,
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),
    clustering_method_rows = "ward.D",
    clustering_distance_rows = "euclidean",
    clustering_method_columns = "ward.D",
    clustering_distance_columns = "euclidean",
    row_dend_width = unit(1.5, "cm"),
    column_dend_height = unit(1.5, "cm"),,
    show_row_names = TRUE,  # Show row names
    row_title_gp = gpar(fontsize = 13, fontface = "bold"),
    row_names_gp = gpar(fontsize = 10, fontface = "bold"),
    left_annotation = row_annotation,
    show_row_dend = TRUE, show_column_dend = TRUE,
    column_title = "Composition of cell types across patients",
    row_title = "Patients",
    column_names_rot = 55,
    #rect_gp = gpar(col = "grey90"), # possibility to add boarders around each cell
    border = "white"
  )
  
  # Draw the heatmap to finalize clustering computation
  hm <- draw(heatmap)
  return(hm)
}
heatmap_all_cell_types(df, PFS_data, PFS_color)

#########

# Heatmap for globalcelltypes
heatmap_global_cell_types <- function(df, PFS_data, PFS_color) {
  # Aggregate data by imageid and phenotype and calculate percentages
  agg_data <- df %>%
    group_by(imageid, GlobalCellType) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(imageid) %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ungroup()
  
  # Spread data to wide format for heatmap input
  wide_data <- agg_data %>%
    select(imageid, GlobalCellType, percentage) %>%
    pivot_wider(names_from = GlobalCellType, values_from = percentage, values_fill = list(percentage = 0))
  
  # Merge with PFS and PFS_progression data
  merged_data <- left_join(wide_data, PFS_data, by = "imageid")
  color_mapping <- colorRamp2(c(vmin, vmax), c("white", "#D73027"))
  
  # Prepare data for heatmap
  heatmap_matrix <- as.matrix(merged_data %>% select(-imageid, -PFS_33p))
  rownames(heatmap_matrix) <- merged_data$imageid  # Set row names to imageid
  
  row_annotations <- merged_data %>% select(PFS_33p) %>% as.data.frame()
  
  
  # Create row annotation with PFS and PFS_progression
  row_annotation <- rowAnnotation(
    PFS = row_annotations$PFS_33p,
    col = list(PFS = PFS_color),
    #annotation_name = c(PFS = "PFS", PFS_progression = "PFS value"),
    annotation_name_rot = 55,
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
    annotation_legend_param = list(
      PFS = list(title = "PFS"))
  )
  
  # Create heatmap
  heatmap <- Heatmap(
    heatmap_matrix,
    name = "% of cell type per patient",
    col = color_mapping,
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),
    clustering_method_rows = "ward.D",
    clustering_distance_rows = "euclidean",
    clustering_method_columns = "ward.D",
    clustering_distance_columns = "euclidean",
    row_dend_width = unit(1.5, "cm"),
    column_dend_height = unit(1.5, "cm"),,
    show_row_names = TRUE,  # Show row names
    row_title_gp = gpar(fontsize = 13, fontface = "bold"),
    row_names_gp = gpar(fontsize = 10, fontface = "bold"),
    left_annotation = row_annotation,
    show_row_dend = TRUE, show_column_dend = TRUE,
    column_title = "Composition of cell types across patients",
    row_title = "Patients",
    column_names_rot = 55,
    #rect_gp = gpar(col = "grey90"), # possibility to add boarders around each cell
    border = "white"
  )
  
  # Draw the heatmap to finalize clustering computation
  hm <- draw(heatmap)
  return(hm)
}
heatmap_global_cell_types(df, PFS_data, PFS_color)

#_______________________________________________________________________________
# Heatmaps and barplots ordered based on heatmap for global cell types and immune cell types

# Global cell types: heatmap with PFS annotations and barplot ordered based on the heatmap, used for FIGURE 1
globalcelltype_combined_heatmap_barplot <- function(df, PFS_data, PFS_color, cell_type_colors) {
  # Aggregate data by imageid and phenotype and calculate percentages
  agg_data <- df %>%
    group_by(imageid, GlobalCellType) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(imageid) %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ungroup()
  
  # Spread data to wide format for heatmap input
  wide_data <- agg_data %>%
    select(imageid, GlobalCellType, percentage) %>%
    pivot_wider(names_from = GlobalCellType, values_from = percentage, values_fill = list(percentage = 0))
  
  # Merge with PFS and PFS_progression data
  merged_data <- left_join(wide_data, PFS_data, by = "imageid")
  color_mapping <- colorRamp2(c(vmin, vmax), c("white", "#D73027"))
  
  # Prepare data for heatmap
  heatmap_matrix <- as.matrix(merged_data %>% select(-imageid, -PFS_33p))
  rownames(heatmap_matrix) <- merged_data$imageid  # Set row names to imageid
  
  row_annotations <- merged_data %>% select(PFS_33p) %>% as.data.frame()
  
  # Create row annotation with PFS and PFS_progression
  row_annotation <- rowAnnotation(
    PFS = row_annotations$PFS_33p,
    col = list(PFS = PFS_color),
    annotation_name_rot = 55,
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
    annotation_legend_param = list(PFS = list(title = "PFS"))
  )
  
  # Create heatmap
  heatmap <- Heatmap(
    heatmap_matrix,
    name = "% of global cell type per patient",
    col = color_mapping,
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),
    clustering_method_rows = "ward.D",
    clustering_distance_rows = "euclidean",
    clustering_method_columns = "ward.D",
    clustering_distance_columns = "euclidean",
    row_dend_width = unit(1.5, "cm"),
    column_dend_height = unit(1.5, "cm"),
    show_row_names = TRUE,  # Show row names
    row_title_gp = gpar(fontsize = 13, fontface = "bold"),
    row_names_gp = gpar(fontsize = 10, fontface = "bold"),
    left_annotation = row_annotation,
    show_row_dend = TRUE, show_column_dend = TRUE,
    column_title = "Global cell type composition across patients",
    row_title = "Patients",
    column_names_rot = 55,
    #rect_gp = gpar(col = "grey90"), # possibility to add boarders around each cell
    border = "white"
  )
  
  # Draw the heatmap to finalize clustering computation
  hm <- draw(heatmap)
  
  # Extract the order of imageids from the heatmap
  row_order <- row_order(hm)
  ordered_image_ids <- rownames(heatmap_matrix)[row_order]
  ordered_image_ids <- rev(ordered_image_ids)
  
  col_order <- column_order(hm)
  ordered_phenotypes <- colnames(heatmap_matrix)[col_order]
  ordered_phenotypes <- rev(ordered_phenotypes)
  
  # Function to draw barplot based on the order from the heatmap, also ordered by tumor %
  barplot_based_on_heatmap <- function(df, ordered_phenotypes, cell_type_colors) {
    # Convert 'phenotype' to factor with the desired order
    df$GlobalCellType <- factor(df$GlobalCellType, levels = ordered_phenotypes)
    
    # Calculate the proportions of each cell type within each imageid
    proportion_data <- df %>%
      group_by(imageid, GlobalCellType) %>%
      summarise(count = n(), .groups = 'drop') %>%
      group_by(imageid) %>%
      mutate(total_count = sum(count)) %>%
      ungroup() %>%
      mutate(Proportion = (count / total_count) * 100)
    
    # Order image ids based on the heatmap
    proportion_data$imageid <- factor(proportion_data$imageid, levels = ordered_image_ids)
    
    # Create a data frame from the color dictionary
    color_order <- cell_type_colors[levels(proportion_data$GlobalCellType)]
    
    # Create the composition plot
    composition_plot <- ggplot(proportion_data, aes(y = factor(imageid, levels = ordered_image_ids), x = Proportion, fill = GlobalCellType)) +
      geom_bar(stat = "identity", position = "stack", 
               width = 0.9,  # Adjust width for spacing between bars
               color = "NA"  # Set bar border color
      ) + 
      labs(x = "Proportion", y = "") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12),
            axis.text.y = element_text(size = 10, face = "bold"),
            axis.title = element_text(size = 16),
            plot.title = element_text(size = 18, face = "bold"),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 12, face = "bold")) +
      ggtitle("") +
      scale_fill_manual(values = color_order, name = "Global Cell Type")
    
    return(composition_plot)
  }
  
  # Draw the barplot based on the order from the heatmap
  barplot <- barplot_based_on_heatmap(df, ordered_phenotypes, cell_type_colors)
  
  return(list(heatmap = hm, barplot = barplot))
}
result <- globalcelltype_combined_heatmap_barplot(df, PFS_data, PFS_color, all_celltype_colors)

# Plot the heatmap and barplot
print(result$heatmap)
print(result$barplot)

##### Globalcelltypes heatmap with endothelial separately

# Global cell types + endothelial
df$Global_endothelial <- ifelse(
  df$GlobalCellType == "Stromal", 
  ifelse(grepl("Endothelial", df$phenotype), "Endothelial", "Stromal"), 
  df$GlobalCellType)

vmin=0
vmax=45

heatmap_global_cell_types <- function(df, PFS_data, PFS_color) {
  # Aggregate data by imageid and phenotype and calculate percentages
  agg_data <- df %>%
    group_by(imageid, Global_endothelial) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(imageid) %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ungroup()
  
  # Spread data to wide format for heatmap input
  wide_data <- agg_data %>%
    select(imageid, Global_endothelial, percentage) %>%
    pivot_wider(names_from = Global_endothelial, values_from = percentage, values_fill = list(percentage = 0))
  
  # Merge with PFS and PFS_progression data
  merged_data <- left_join(wide_data, PFS_data, by = "imageid")
  color_mapping <- colorRamp2(c(vmin, vmax), c("white", "#D73027"))
  
  # Prepare data for heatmap
  heatmap_matrix <- as.matrix(merged_data %>% select(-imageid, -PFS_33p))
  rownames(heatmap_matrix) <- merged_data$imageid  # Set row names to imageid
  
  row_annotations <- merged_data %>% select(PFS_33p) %>% as.data.frame()
  
  
  # Create row annotation with PFS and PFS_progression
  row_annotation <- rowAnnotation(
    PFS = row_annotations$PFS_33p,
    col = list(PFS = PFS_color),
    #annotation_name = c(PFS = "PFS", PFS_progression = "PFS value"),
    annotation_name_rot = 55,
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
    annotation_legend_param = list(
      PFS = list(title = "PFS"))
  )
  
  # Create heatmap
  heatmap <- Heatmap(
    heatmap_matrix,
    name = "% of cell type per patient",
    col = color_mapping,
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),
    clustering_method_rows = "ward.D",
    clustering_distance_rows = "euclidean",
    clustering_method_columns = "ward.D",
    clustering_distance_columns = "euclidean",
    row_dend_width = unit(1.5, "cm"),
    column_dend_height = unit(1.5, "cm"),,
    show_row_names = TRUE,  # Show row names
    row_title_gp = gpar(fontsize = 13, fontface = "bold"),
    row_names_gp = gpar(fontsize = 10, fontface = "bold"),
    left_annotation = row_annotation,
    show_row_dend = TRUE, show_column_dend = TRUE,
    column_title = "Composition of cell types across patients",
    row_title = "Patients",
    column_names_rot = 55,
    #rect_gp = gpar(col = "grey90"), # possibility to add boarders around each cell
    border = "white"
  )
  
  # Draw the heatmap to finalize clustering computation
  hm <- draw(heatmap)
  return(hm)
}
heatmap_global_cell_types(df, PFS_data, PFS_color)


########
# Heatmap and barplot for immune subset(imageids and immune cells ordered based on heatmap) used for FIGURE 1
immune_combined_heatmap_barplot <- function(df, PFS_data, PFS_color, cell_type_colors) {
  # Subset immune cells
  df_immune <- df %>% filter(!phenotype %in% c("Tumor", "Stromal", "Endothelial"))
  # Aggregate data by imageid and phenotype and calculate percentages
  agg_data <- df_immune %>%
    group_by(imageid, phenotype) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(imageid) %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ungroup()
  
  # Spread data to wide format for heatmap input
  wide_data <- agg_data %>%
    select(imageid, phenotype, percentage) %>%
    pivot_wider(names_from = phenotype, values_from = percentage, values_fill = list(percentage = 0))
  
  # Merge with PFS and PFS_progression data
  merged_data <- left_join(wide_data, PFS_data, by = "imageid")
  color_mapping <- colorRamp2(c(vmin, vmax), c("white", "#D73027"))
  
  # Prepare data for heatmap
  heatmap_matrix <- as.matrix(merged_data %>% select(-imageid, -PFS_33p))
  rownames(heatmap_matrix) <- merged_data$imageid  # Set row names to imageid
  
  row_annotations <- merged_data %>% select(PFS_33p) %>% as.data.frame()
  
  
  # Create row annotation with PFS and PFS_progression
  row_annotation <- rowAnnotation(
    PFS = row_annotations$PFS_33p,
    col = list(PFS = PFS_color),
    #annotation_name = c(PFS = "PFS", PFS_progression = "PFS value"),
    annotation_name_rot = 55,
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
    annotation_legend_param = list(PFS = list(title = "PFS"))
    )
  
  # Create heatmap
  heatmap <- Heatmap(
    heatmap_matrix,
    name = "% of immune subtype per patient",
    col = color_mapping,
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    column_title_gp = gpar(fontsize = 13, fontface = "bold"),
    clustering_method_rows = "ward.D",
    clustering_distance_rows = "euclidean",
    clustering_method_columns = "ward.D",
    clustering_distance_columns = "euclidean",
    row_dend_width = unit(1.5, "cm"),
    column_dend_height = unit(1.5, "cm"),
    show_row_names = TRUE,  # Show row names
    row_title_gp = gpar(fontsize = 13, fontface = "bold"),
    row_names_gp = gpar(fontsize = 10, fontface = "bold"),
    left_annotation = row_annotation,
    show_row_dend = TRUE, show_column_dend = TRUE,
    column_title = "Immune subtype composition across patients",
    row_title = "Patients",
    column_names_rot = 55,
    #rect_gp = gpar(col = "grey90"), # possibility to add boarders around each cell
    border = "white"
  )
  
  # Draw the heatmap to finalize clustering computation
  hm <- draw(heatmap)
  
  # Extract the order of imageids from the heatmap
  row_order <- row_order(hm)
  ordered_image_ids <- rownames(heatmap_matrix)[row_order]
  ordered_image_ids <- rev(ordered_image_ids)
  
  col_order <- column_order(hm)
  ordered_phenotypes <- colnames(heatmap_matrix)[col_order]
  ordered_phenotypes <- rev(ordered_phenotypes)
  
  # Function to draw barplot based on the order from the heatmap, also ordered by tumor %
  barplot_based_on_heatmap <- function(df_immune, ordered_phenotypes, cell_type_colors) {
    # Convert 'phenotype' to factor with the desired order
    df_immune$phenotype <- factor(df_immune$phenotype, levels = ordered_phenotypes)
    
    # Calculate the proportions of each cell type within each imageid
    proportion_data <- df_immune %>%
      group_by(imageid, phenotype) %>%
      summarise(count = n(), .groups = 'drop') %>%
      group_by(imageid) %>%
      mutate(total_count = sum(count)) %>%
      ungroup() %>%
      mutate(Proportion = (count / total_count) * 100)
    
    # Order image ids based on the heatmap
    proportion_data$imageid <- factor(proportion_data$imageid, levels = ordered_image_ids)
    
    # Create a data frame from the color dictionary
    color_order <- cell_type_colors[levels(proportion_data$phenotype)]
    
    # Create the composition plot
    composition_plot <- ggplot(proportion_data, aes(y = factor(imageid, levels = ordered_image_ids), x = Proportion, fill = phenotype)) +
      geom_bar(stat = "identity", position = "stack", 
               width = 0.9,  # Adjust width for spacing between bars
               color = "NA"  # Set bar border color
      ) + 
      labs(x = "Proportion", y = "") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 0, hjust = 1, size = 12),
            axis.text.y = element_text(size = 10, face = "bold"),
            axis.title = element_text(size = 16),
            plot.title = element_text(size = 18, face = "bold"),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 12, face = "bold")) +
      ggtitle("") +
      scale_fill_manual(values = color_order, name = "Cell Type")
    
    return(composition_plot)
  }
  
  # Draw the barplot based on the order from the heatmap
  barplot <- barplot_based_on_heatmap(df_immune, ordered_phenotypes, cell_type_colors)
  
  return(list(heatmap = hm, barplot = barplot))
}
result <- immune_combined_heatmap_barplot(df, PFS_data, PFS_color, all_celltype_colors)

# Plot the heatmap and barplot
print(result$heatmap)
print(result$barplot)

