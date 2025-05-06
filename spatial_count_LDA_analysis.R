# Code for analysing RCNs computed with spatial count and LDA

# Author: Alva Grönholm
# Creation date: 12.1.2025
# Modification date: 5.5.2025

# Purpose
# This code generates:
# 1. Barplots of phenotypes in RCNs and RCN cell counts
# 2. Heatmap of phenotype proportions across RCNs
# 3. Z-scored heatmaps (by cluster and sample) showing phenotype distributions per patients
# 4. Boxplots of RCN % in clinical groups (short and long PFS)
# 5. Ranc abundance curves of phenotypes in RCN


#_______________________________________________________________________________
#Load libraries
library(ggplot2)
library(RColorBrewer)
library(ComplexHeatmap)
library(dplyr)
library(ggpubr)


#_______________________________________________________________________________

#Renaming phenotypes
df <- df %>%
  mutate(phenotype = recode(phenotype,
                            "CD163_neg_Macrophages" = "CD163- Macrophages",
                            "CD163_pos_Macrophages" = "CD163+ Macrophages",
                            "CD4_T_cells" = "CD4+ T cells",
                            'CD8_T_cells' = 'CD8+ T cells',
                            'B_cells' = 'B cells',
                            'CD163_pos_DC' = 'CD163+ DC',
                            'CD4_FOXP3_T_cells' = 'FOXP3+ CD4+ T cells'))


# Annotations of spatial_count 13 microns clusters
df <- df %>%
  mutate(spatial_count_40_k15 = recode(spatial_count_40_k15,
                                    "0" = "RCN 04: TSI - DC",
                                    "1" = "RCN 01: Tumor core",
                                    "2" = "RCN 05: Stroma, immune excluded",
                                    '3' = 'RCN 06: Stroma 1',
                                    '4' = 'RCN 08: Endothelial',
                                    '5' = 'RCN 03: TSI - stroma',
                                    '6' = 'RCN 09: DC',
                                    '7' = 'RCN 10: CD163+ DC',
                                    '8' = 'RCN 13: CD8+ T cells',
                                    '9' = 'RCN 02: TSI - tumor',
                                    '10' = 'RCN 11: CD163- Macrophages',
                                    '11' = 'RCN 12: CD163+ Macrophages',
                                    '12' = 'RCN 14: CD4+ T cells',
                                    '13' = 'RCN 07: Stroma 2',
                                    '14' = 'RCN 15: B-cells'))

# Annotations of spatial_count 25 microns clusters
df <- df %>%
  mutate(spatial_count_77_k12 = recode(spatial_count_77_k12,
                                       "0" = "RCN 01: Tumor core",
                                       "1" = "RCN 06: Stroma, immune infiltrated",
                                       "2" = "RCN 03: TSI - stroma",
                                       '3' = 'RCN 09: DC',
                                       '4' = 'RCN 04: TSI - DC',
                                       '5' = 'RCN 10: CD163+ DC',
                                       '6' = 'RCN 11: CD8+ T cells',
                                       '7' = 'RCN 12: B-cells',
                                       '8' = 'RCN 07: Endothelial',
                                       '9' = 'RCN 02: TSI - tumor',
                                       '10' = 'RCN 08: Immune compartment',
                                       '11' = 'RCN 05: Stroma, immune excluded'))

# Annotations of spatial_count 35 microns clusters
df <- df %>%
  mutate(spatial_count_108_k13 = recode(spatial_count_108_k13,
                                       "0" = "RCN 12: CD8+ T cells",
                                       "1" = "RCN 05: TSI - DC",
                                       "2" = "RCN 07: Stroma, immune infiltrated",
                                       '3' = 'RCN 01: Tumor core',
                                       '4' = 'RCN 04: TSI - stroma',
                                       '5' = 'RCN 09: Immune compartment',
                                       '6' = 'RCN 11: CD163+ DC',
                                       '7' = 'RCN 10: DC',
                                       '8' = 'RCN 13: B-cells',
                                       '9' = 'RCN 03: TSI - tumor 2',
                                       '10' = 'RCN 02: TSI - tumor 1',
                                       '11' = 'RCN 08: Endothelial',
                                       '12' = 'RCN 06: Stroma, immune excluded'))

# Annotations of LDA 13 microns clusters
df <- df %>%
  mutate(lda_40_k13 = recode(lda_40_k13,
                                        "0" = "RCN 04: Stroma, immune excluded",
                                        "1" = "RCN 02: TSI - tumor",
                                        "2" = "RCN 08: DC",
                                        '3' = 'RCN 11: B-cells and Macrophages',
                                        '4' = 'RCN 09: CD163+ DC',
                                        '5' = 'RCN 05: Stroma 1',
                                        '6' = 'RCN 12: CD8+ T cells',
                                        '7' = 'RCN 06: Stroma - CD4+ T cells',
                                        '8' = 'RCN 01: Tumor core',
                                        '9' = 'RCN 10: CD163+ Macrophages',
                                        '10' = 'RCN 07: Endothelial',
                                        '11' = 'RCN 13: CD4+ T cells',
                                        '12' = 'RCN 03: TSI - DC'))

# Annotations of LDA 25 microns clusters
df <- df %>%
  mutate(lda_77_k15 = recode(lda_77_k15,
                             "0" = "RCN 13: CD4+ T cells",
                             "1" = "RCN 05: TSI - CD163+ DC",
                             "2" = "RCN 08: Stroma, immune infiltrated",
                             '3' = 'RCN 09: Endothelial',
                             '4' = 'RCN 06: TSI - macrophages',
                             '5' = 'RCN 15: B-cells',
                             '6' = 'RCN 10: DC',
                             '7' = 'RCN 02: TSI - tumor',
                             '8' = 'RCN 14: CD8+ T cells',
                             '9' = 'RCN 03: TSI - stroma',
                             '10' = 'RCN 07: Stroma, immune excluded',
                             '11' = 'RCN 12: Macrophages',
                             '12' = 'RCN 04: TSI - DC',
                             '13' = 'RCN 01: Tumor core',
                             '14' = 'RCN 11: CD163+ DC'))

# Annotations of LDA 35 microns clusters
df <- df %>%
  mutate(lda_108_k13 = recode(lda_108_k13,
                             "0" = "RCN 11: Macrophages",
                             "1" = "RCN 02: TSI 1",
                             "2" = "RCN 08: Endothelial",
                             '3' = 'RCN 06: Stroma, immune excluded',
                             '4' = 'RCN 10: B-cells and CD163+ DC',
                             '5' = 'RCN 03: TSI 2',
                             '6' = 'RCN 12: CD8+ T cells',
                             '7' = 'RCN 05: TSI - myeloid',
                             '8' = 'RCN 01: Tumor core',
                             '9' = 'RCN 09: DC',
                             '10' = 'RCN 07: Stroma, immune infiltrated',
                             '11' = 'RCN 13: CD4+ T cells',
                             '12' = 'RCN 04: TSI - stroma'))

#_______________________________________________________________________________ 

# Colors
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
  "Endothelial" = "#EB881F",
  "Immune" = "#DA3683"
)

RCN_colors <- c("Stroma, immune infiltrated"= "#EB881F", 
                'Tumor, no infiltration'= "#B1CCE3", 
                'Tumor-Stromal Interface'= "#d62728",
                'DC'= "#b2df8a", 
                'Stroma, not infiltrated'= "#fffd83", 
                'Endothelial'= "#C39BD3",
                'CD8+ T cells'= "#EAA09A", 
                'CD163+ DC'= "#2E7629", 
                'CD4+ T cells'= "#DA3683", 
                'B-cell dominated'= "#393b79",
                'CD163- Macrophages'= "#FAE5D3")


PFS_color <- c("Short" = "#A4C8E1", "Middle class" = "#0076A8", "Long" = "#003D5B")


#_______________________________________________________________________________
# Barplots of RCNs

# Barplot function with configurable x-axis and fill variable
RCN_barplot <- function(df, x_axis = "RCN", fill_variable = "phenotype", colors = NULL) {
  # Calculate proportions based on the chosen fill variable
  proportions <- df %>%
    group_by(!!sym(x_axis), !!sym(fill_variable)) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(!!sym(x_axis)) %>%
    mutate(percentage = count / sum(count)) %>%
    ungroup()
  
  # Create the bar plot with configurable x-axis and fill variable
  ggplot(df, aes_string(x = x_axis, fill = fill_variable)) +
    geom_bar(position = "fill", color = 'black') + 
    scale_fill_manual(values = colors) +
    labs(title = "Proportion of phenotypes in clusters, LDA 35 microns",
         x = "Clusters",
         y = "Proportion",
         fill = "Phenotype") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 11, colour = 'black', face = "bold"))
}

g <- RCN_barplot(df, x_axis = "spatial_count_108_k13", fill_variable = "GlobalCellType", colors = all_celltype_colors)


# Countplot for RCNs
# Create a barplot for the total count of each RCN
RCN_count_barplot <- function(df, x_axis = "RCN") {
  # Count occurrences per RCN
  count_data <- df %>%
    group_by(!!sym(x_axis)) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Create the bar plot
  ggplot(count_data, aes_string(x = x_axis, y = "count")) +
    geom_bar(stat = "identity", fill = "gray", color = 'black') +  # Gray bars for count
    labs(x = "Clusters", y = "N (cells)", title = "Cell Counts per RCN") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1),  # Hide x-axis labels to align
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    scale_y_continuous(breaks = scales :: pretty_breaks(n=5))
}

RCN_count_barplot(df, x_axis = 'spatial_count_108_k13')


#_______________________________________________________________________________

# Barplot plotting phenotype - imageid for each cluster

# Filter out Middle class
df_filtered <- df_addedsamples %>% filter(PFS_33p != 'Middle class')

# Calculate proportions
df_proportions <- df_filtered %>% group_by(cluster_kmeans_10, imageid, PFS_33p, phenotype) %>%
  summarize(count = n(), .groups = 'drop') %>%
  group_by(cluster_kmeans_10, imageid) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

# Order imageid based on PFS
df_proportions <- df_proportions %>% mutate(imageid = factor(imageid, levels = unique(imageid[order(PFS_33p == 'Short', decreasing = TRUE)])))

# Plot for each cluster
ggplot(df_proportions, aes(x = percentage, y=imageid, fill = phenotype)) +
  geom_bar(stat = 'identity', color = "black") +
  facet_wrap(~ cluster_kmeans_10, scales = 'free_y') +
  scale_y_discrete(limits = rev) +
  scale_fill_manual(values = all_celltype_colors) +
  labs(
    title = "Percentage of phenotypes per patient in each cluster",
    x = "Percentage",
    y = "Patients",
    fill = "Phenotype"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(size = 10),
        axis.text.y = element_text(size = 8))


# Split the data by cluster
clusters <- split(df_proportions, df_proportions$cluster_kmeans_10)

for (cluster_num in names(clusters)) {cluster_data <- clusters[[cluster_num]]
  
  # Arrange data to ensure Short is at the top and Long at the bottom
  cluster_data <- cluster_data %>%
    arrange(PFS_33p == "Short", desc(imageid)) %>%
    mutate(imageid = factor(imageid, levels = unique(imageid)))
  
  # Plot with black boxes
  p <- ggplot(cluster_data, aes(x = percentage, y = imageid, fill = phenotype)) +
    geom_bar(stat = "identity", color = "black") + 
    scale_y_discrete(limits = rev) + # Reverse y-axis for desired order
    scale_fill_manual(values = all_celltype_colors) + 
    labs(
      title = paste("Cluster", cluster_num),
      x = "Percentage",
      y = "Patients",
      fill = "Phenotype"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 9),
      plot.title = element_text(size = 14, face = "bold")
    )
  
  # Print the plot to view
  print(p) 
  
  ggsave(
    filename = paste0("", cluster_num, "_Barplot.pdf"),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300
  )
}

#_______________________________________________________________________________
# Heatmap of RCN and phenotype

# Heatmap colors
vmin=3
vmax=65

# Heatmap for phneotypes in RCN
heatmap_RCN <- function(df, RCN, x_axis) {
  df <- df %>% mutate(RCN = RCN)
  df <- df %>% mutate(x_axis = x_axis)
  # Aggregate data by RCN and phenotype and calculate percentages
  agg_data <- df %>%
    group_by(x_axis, RCN) %>% 
    summarise(count = n(), .groups = 'drop') %>%
    group_by(RCN) %>%                  
    mutate(percentage = count / sum(count) * 100) %>% 
    ungroup()                          
  
  # Spread data to wide format for heatmap input
  wide_data <- agg_data %>%
    select(RCN, x_axis, percentage) %>%
    pivot_wider(names_from = x_axis, values_from = percentage, values_fill = list(percentage = 0))
  
  color_mapping <- colorRamp2(c(vmin, vmax), c("white", "#cc2127"))
  
  # Prepare data for heatmap
  heatmap_matrix <- as.matrix(wide_data %>% select(-RCN))
  rownames(heatmap_matrix) <- wide_data$RCN  # Set row names to RCN
  
  custom_order <- c("Tumor, no infiltration", 'Tumor-Stromal Interface', 'Stroma, not infiltrated', 'Stroma, immune infiltrated', 'Endothelial', 'B-cell dominated','CD8+ T cells',
                    'DC', 'CD163+ DC', 'CD163- Macrophages')
  
  # Create heatmap
  heatmap <- Heatmap(
    heatmap_matrix,
    name = "% of phenotype per RCN",
    col = color_mapping,
    cluster_rows = FALSE,
    #cluster_columns = FALSE,
    #row_order = custom_order,
    #column_order = c('Tumor', 'Stromal', 'Endothelial', 'B cells', 'CD8+ T cells', 'CD4+ T cells', 'FOXP3+ CD4+ T cells', 'DC', 'CD163+ DC','CD163- Macrophages', 'CD163+ Macrophages'),
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
    show_row_dend = TRUE, show_column_dend = TRUE,
    column_title = "LDA: 35 microns",
    row_title = "RCN",
    column_names_rot = 55,
    border = "white"
  )
  
  # Draw the heatmap to finalize clustering computation
  hm <- draw(heatmap)
  return(hm)
}

heatmap_RCN(df, RCN = df$lda_108_k10, x_axis = df$phenotype)


png("F:/Alva/Images and figures/spatial_analysis/comparison_LDA_vs_spatial_count/HM_spatial_count_phenotype.png", width = 10, height = 8, units = "in", res = 300)
draw(p)
dev.off()

#_______________________________________________________________________________

# Heatmap with Z-scoring - (Z-scored for clusters)
heatmap_zscore_clusters <- function(df, RCN, x_axis, PFS_data, PFS_col, cluster_col) {
  # Ensure RCN and x_axis are part of the dataframe
  df <- df %>% mutate(RCN = RCN, x_axis = x_axis)
  # Exclude the "Middle class" group
  selected_imageids <- PFS_data$imageid[PFS_data[[PFS_col]] != "Middle class"]
  group_data <- df %>% filter(x_axis %in% selected_imageids)
  # Aggregate and calculate percentages
  agg_data <- group_data %>%
    group_by(RCN, x_axis) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(x_axis) %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ungroup()
  # Spread data to wide format
  wide_data <- agg_data %>%
    select(RCN, x_axis, percentage) %>%
    pivot_wider(names_from = x_axis, values_from = percentage, values_fill = list(percentage = 0))
  # Ensure percentages are numeric
  percentages_matrix <- as.matrix(wide_data %>% select(-RCN))
  # Z-score the data
  percentages_zscored <- t(scale(t(percentages_matrix)))
  rownames(percentages_zscored) <- wide_data$RCN  
  # Define annotation for `PFS_33p` column
  pfs_annotations <- PFS_data %>%
    filter(imageid %in% colnames(percentages_zscored)) %>%
    arrange(match(imageid, colnames(percentages_zscored)))  # Ensure ordering matches the heatmap columns
  pfs_colors <- list(PFS_33p = c("Short" = "#A4C8E1", "Long" = "#003D5B"))
  annotations <- HeatmapAnnotation(
    PFS_33p = pfs_annotations[[PFS_col]],
    col = pfs_colors,
    annotation_name_side = "left",
    annotation_legend_param = list(title = "PFS Group", labels_gp = gpar(fontsize = 10))
  )
  # Define color gradient for the heatmap
  color_gradient <- colorRamp2(c(-2,-1, 0, 1, 2), c("#3b4cc0", '#8db0fe', 'white','#f4987a', '#b40426'))
  # Create the heatmap
  hmap <- Heatmap(
    percentages_zscored,
    name = "Cluster % per patient (Z-scored)",
    top_annotation = annotations,
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10),
    column_title = "Patients",
    row_title = "Clusters",
    row_names_side = "left",
    row_dend_side = "right",
    col = color_gradient,
    column_names_rot = 50
  )
  # Draw the heatmap
  draw(hmap)
}

heatmap_zscore_clusters(df = df, RCN = df$cluster_lda_10, x_axis = df$imageid, PFS_data = PFS_data, PFS_col = "PFS_33p", cluster_col = "cluster_kmeans_10")

# Heatmap with Z-scoring - (Z-scored for samples)
heatmap_zscore_samples <- function(df, RCN, x_axis, PFS_data, PFS_col, cluster_col) {
  # Ensure RCN and x_axis are part of the dataframe
  df <- df %>% mutate(RCN = RCN, x_axis = x_axis)
  # Exclude the "Middle class" group
  selected_imageids <- PFS_data$imageid[PFS_data[[PFS_col]] != "Middle class"]
  group_data <- df %>% filter(x_axis %in% selected_imageids)
  # Aggregate and calculate percentages
  agg_data <- group_data %>%
    group_by(RCN, x_axis) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(x_axis) %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ungroup()
  # Spread data to wide format
  wide_data <- agg_data %>%
    select(RCN, x_axis, percentage) %>%
    pivot_wider(names_from = x_axis, values_from = percentage, values_fill = list(percentage = 0))
  # Ensure percentages are numeric
  percentages_matrix <- as.matrix(wide_data %>% select(-RCN))
  # Z-score the data
  percentages_zscored <- (scale((percentages_matrix)))
  rownames(percentages_zscored) <- wide_data$RCN 
  # Define annotation for `PFS_33p` column
  pfs_annotations <- PFS_data %>%
    filter(imageid %in% colnames(percentages_zscored)) %>%
    arrange(match(imageid, colnames(percentages_zscored)))  # Ensure ordering matches the heatmap columns
  pfs_colors <- list(PFS_33p = c("Short" = "#A4C8E1", "Long" = "#003D5B"))
  annotations <- HeatmapAnnotation(
    PFS_33p = pfs_annotations[[PFS_col]],
    col = pfs_colors,
    annotation_name_side = "left",
    annotation_legend_param = list(title = "PFS Group", labels_gp = gpar(fontsize = 10))
  )
  # Define color gradient for the heatmap
  color_gradient <- color_gradient <- colorRamp2(c(-2,-1, 0, 1, 2), c("#3b4cc0", '#8db0fe', 'white','#f4987a', '#b40426'))
  # Create the heatmap
  hmap <- Heatmap(
    percentages_zscored,
    name = "Cluster % per patient (Z-scored), Spatial count, 13",
    cluster_rows = FALSE,
    top_annotation = annotations,
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10),
    column_title = "Patients",
    row_title = "Clusters",
    row_names_side = "left",
    row_dend_side = "right",
    col = color_gradient,
    column_names_rot = 50
  )
  # Draw the heatmap
  draw(hmap)
}

h <- heatmap_zscore_samples(df = df, RCN = df$spatial_count_40_k15, x_axis = df$imageid, PFS_data = PFS_data, PFS_col = "PFS_33p")


pdf("F:/Alva/Final_spatial_analysis_plots/z_scored_Heatmap_spatial_count_13.pdf", width = 12, height = 8)
draw(h)
dev.off()

#_______________________________________________________________________________
# Clinical data

# Heatmap colors
vmin=0
vmax=45 

# Heatmap for all cell types with PFS annotation - removes Middle class
heatmap_cluster <- function(df, PFS_data, PFS_color) {
  # Remove the Middle class from both df and PFS_data
  df <- df %>% filter(PFS_33p != "Middle class")
  PFS_data <- PFS_data %>% filter(PFS_33p != "Middle class")
  # Aggregate data by imageid and cluster and calculate percentages
  agg_data <- df %>%
    group_by(imageid, cluster_kmeans_10) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(imageid) %>%
    mutate(percentage = count / sum(count) * 100) %>%
    ungroup()
  
  # Spread data to wide format for heatmap input
  wide_data <- agg_data %>%
    select(imageid, cluster_kmeans_10, percentage) %>%
    pivot_wider(names_from = cluster_kmeans_10, values_from = percentage, values_fill = list(percentage = 0))
  
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
    border = "white"
  )
  
  # Draw the heatmap to finalize clustering computation
  hm <- draw(heatmap)
  return(hm)
}
heatmap_cluster(df_residual, PFS_data, PFS_color)

#_______________________________________________________________________________
##Boxplots for a single RCN

# Calculate percentage of each cluster per image
percentage_df <- df %>%
  group_by(imageid, PFS_33p, spatial_count_40_k15) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count)) %>%
  ungroup()

# Filter to include only the desired cluster
cluster_of_interest <- "RCN 06: Stroma 1" 

boxplot_df <- percentage_df %>%
  filter(spatial_count_40_k15 == cluster_of_interest)

boxplot_df <- boxplot_df%>%filter(PFS_33p != "Middle class")

# Define the custom theme
custom_theme <- theme(panel.border = element_rect(colour = "black", linewidth = 1, fill = NA), 
                      panel.background = element_blank(), 
                      plot.title = element_text(hjust = 0.5))

# Create the boxplot
a <- ggplot(boxplot_df, aes(x = factor(PFS_33p, levels = c("Short", "Long")), y = percentage, fill = PFS_33p)) + 
  geom_boxplot(lwd = 1, alpha = 0.9) + 
  geom_jitter(position = position_jitter(0.2), size = 1, alpha = 0.9) + 
  scale_fill_manual(values = PFS_color) + 
  scale_y_continuous(limits = c(0, 0.3)) +
  labs(title = paste("Percentage of", cluster_of_interest, "by PFS Status"), 
       x = "PFS Status", 
       y = "% cluster") + 
  stat_compare_means(method = "wilcox.test", label.x = 1.65, label.y.npc = "top", size = 5) + 
  custom_theme + 
  theme(axis.text = element_text(size = 16), 
        axis.title = element_text(size = 16), 
        aspect.ratio = 1, 
        legend.position = "none")


###################
# Boxplots for all RCNs

# Calculate percentage of each cluster per image
percentage_df <- df %>%
  group_by(imageid, PFS_33p, lda_108_k13) %>%
  summarise(count = n()) %>%
  mutate(percentage = count / sum(count)) %>%
  ungroup()

# Filter out undesired classes and images
percentage_df <- percentage_df %>%
  filter(PFS_33p != "Middle class") 

# Define the custom theme
custom_theme <- theme(
  panel.border = element_rect(colour = "black", linewidth = 1, fill = NA),
  panel.background = element_blank(),
  plot.title = element_text(hjust = 0.5),
  axis.text = element_text(size = 12),
  axis.title = element_text(size = 14),
  strip.text = element_text(size = 14),
  legend.position = "none"
)

# Perform Wilcoxon tests for each cluster
p_values <- percentage_df %>%
  group_by(lda_108_k13) %>%
  summarise(p_value = wilcox.test(percentage ~ PFS_33p)$p.value) %>%
  mutate(p_label = ifelse(p_value < 0.001, "***", 
                          ifelse(p_value < 0.01, "**", 
                                 ifelse(p_value < 0.05, "*", "ns"))))

# Merge p-values into the main dataframe
percentage_df <- percentage_df %>%
  left_join(p_values, by = "lda_108_k13")

# Create the boxplot with facets for each cluster
ggplot(percentage_df, aes(x = factor(PFS_33p, levels = c("Short", "Long")), y = percentage, fill = PFS_33p)) +
  geom_boxplot(lwd = 1, alpha = 0.9) +
  geom_jitter(position = position_jitter(0.2), size = 1, alpha = 0.9) +
  scale_fill_manual(values = PFS_color) +
  labs(title = "Percentage of Clusters by PFS Status (LDA: radius 35)",
       x = "PFS Status",
       y = "Percentage of Cluster") +
  custom_theme +
  theme(aspect.ratio = 1) +
  facet_wrap(~lda_108_k13, scales = "free_y") +
  stat_compare_means(
    aes(label = paste0("p = ", signif(p_value, 3))),
    method = "wilcox.test",
    label.x.npc = 0.25,
    label.y.npc = 'top',
    size = 4
  ) 
  #+ ylim(0, max(percentage_df$percentage)*1.1)


#_______________________________________________________________________________
# Compare immune cell abundances in stroma no infiltration 

# Filter for immune cells
immune_df <- df %>%
  filter(GlobalCellType == "Immune")

# Create a dataset for LDA RCN 04
lda_df <- immune_df %>%
  filter(lda_40_k13 == "RCN 04: Stroma, immune excluded") %>%
  mutate(RCN = "LDA")

# Create a dataset for Spatial Count RCN 05
spatial_df <- immune_df %>%
  filter(spatial_count_40_k15 == "RCN 05: Stroma, immune excluded") %>%
  mutate(RCN = "Spatial Count")

# Combine datasets
comparison_df <- bind_rows(lda_df, spatial_df)

# Calculate percentage of immune cells per image
percentage_df <- comparison_df %>%
  group_by(imageid, RCN) %>%
  summarise(count = n(), .groups = "drop") %>%
  mutate(percentage = count / sum(count))

# Define the custom theme
custom_theme <- theme(
  panel.border = element_rect(colour = "black", linewidth = 1, fill = NA),
  panel.background = element_blank(),
  plot.title = element_text(hjust = 0.5)
)

# Custom hex color codes
lda_color <- "#eba5d3"  
spatial_color <- "#AeCCe4"  

# Create the boxplot with custom colors
ggplot(percentage_df, aes(x = RCN, y = percentage, fill = RCN)) + 
  geom_boxplot(lwd = 1, alpha = 0.9) + 
  geom_jitter(position = position_jitter(0.2), size = 1, alpha = 0.9) + 
  labs(
    y = "Immune cells (%)"
  ) + 
  stat_compare_means(method = "wilcox.test",  label.x = 1.5, label.y.npc = "top", size = 5) + 
  custom_theme + 
  theme(
    axis.text = element_text(size = 16), 
    axis.title = element_text(size = 16), 
    aspect.ratio = 1, 
    legend.position = "none"
  ) + 
  scale_fill_manual(values = c("LDA" = lda_color, "Spatial Count" = spatial_color))  

#_______________________________________________________________________________
# Rank abundance curve

# Count cluster abundances for both methods
lda_counts <- df %>%
  count(cluster_lda_10) %>%
  rename(Cluster = cluster_lda_10, Abundance = n) %>%
  mutate(Method = "LDA")

spatial_counts <- df %>%
  count(cluster_kmeans_10) %>%
  rename(Cluster = cluster_kmeans_10, Abundance = n) %>%
  mutate(Method = "Spatial count")

# Combine data for ggplot
rac_data <- bind_rows(lda_counts, spatial_counts) %>%
  arrange(Method, desc(Abundance)) %>%  
  group_by(Method) %>%
  mutate(Rank = row_number(),  
         Relative_Abundance = Abundance / sum(Abundance))  

# Convert Cluster to a factor to ensure consistent colors
rac_data$Cluster <- as.factor(rac_data$Cluster)

# Plot the Rank Abundance Curve with consistent cluster colors
ggplot(rac_data, aes(x = Rank, y = Relative_Abundance, color = Cluster, shape = Method)) +
  geom_line(aes(linetype = Method), color = 'black', size =1, alpha = 0.8) +  
  scale_color_manual(values = RCN_colors) +
  geom_point(size = 3) +  
  scale_y_log10() +  
  theme_minimal() +
  labs(title = "Rank Abundance Curve (RAC)",
       x = "Rank (Most Abundant → Least Abundant)",
       y = "Relative Abundance (log scale)",
       color = "Cluster") +
  theme(legend.title = element_text(size = 12),
        legend.position = "right") 


###RAC for a specific RCN

# Set the specific RCN you want to analyze
specific_RCN <- 'Stroma, immune infiltrated'

calculate_counts <- function(df, cluster_cols, specific_RCN) {
  map_dfr(cluster_cols, ~ {
    df %>%
      filter(.data[[.x]] == specific_RCN) %>%  
      count(phenotype, .data[[.x]]) %>% 
      rename(Phenotype = phenotype, Abundance = n) %>%
      mutate(Method = .x, Cluster = .data[[.x]])  
  })
}

# Example usage
cluster_columns <- c("cluster_lda_10", "spatial_count_40_k10", "spatial_count_77_k10", "spatial_count_108_k10", 'lda_77_k10', 'lda_108_k10')
result_counts <- calculate_counts(df, cluster_columns, specific_RCN)

# Combine, sort, and calculate relative abundance dynamically
rac_data <- result_counts %>%
  group_by(Method) %>%
  arrange(Method, desc(Abundance)) %>%
  mutate(Rank = row_number(),  
         Relative_Abundance = Abundance / sum(Abundance))  


# Convert Phenotype to factor for proper x-axis ordering
rac_data$Phenotype <- factor(rac_data$Phenotype, levels = unique(rac_data$Phenotype))

ggplot(rac_data, aes(x = Rank, y = Relative_Abundance, color = Phenotype, shape = Method)) +
  geom_line(aes(group = Method), color = 'black', size = 1, alpha = 0.8) +  
  geom_point(size = 5) +  
  scale_color_manual(values = all_celltype_colors) +  
  scale_y_log10() +  
  theme_minimal() +
  labs(title = paste("Rank Abundance Curve (RAC) for RCN", specific_RCN),
       x = "Rank (Most Abundant → Least Abundant)",
       y = "Relative Abundance (log scale)",
       color = "Phenotype") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  
        legend.title = element_text(size = 16),
        legend.position = "right") +
  scale_shape_manual(values = 18:length(unique(rac_data$Method)))


################################################################################

counts <- df %>% 
  group_by(cluster_kmeans_10, phenotype) %>%
  summarise(count = n(), .groups = 'drop')

counts

################################################################################

# Calculating the number of cells within a radius of 40 pixels

count_radius_40 <- read.csv("Z:/afarkkila_microscopy/Data/CellCycle/analysis/spatial_analysis/spatial_count_40.csv")

mean(count_radius_40$Tumor_spatial_40)
median(count_radius_40$Tumor_spatial_40)



