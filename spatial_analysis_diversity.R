# Code for RCN Diversity Analysis using Shannon and Simpson Indices

# Author: Alva Grönholm
# Creation date: 20.2.2025
# Modification date: 3.5.2025

# Purpose
# This analysis evaluates phenotypic diversity within spatial and LDA-derived clusters across multiple spatial 
# scales using Shannon and Simpson diversity indices. 
# The results are visualized with a scatter plot to compare RCN diversity across methods and radii

#_______________________________________________________________________________

library(vegan)
library(fossil)
library(tidyr)


#_______________________________________________________________________________

# Annotations of spatial_count 35 microns clusters
df <- df %>%
  mutate(spatial_count_108_k13 = recode(spatial_count_108_k13,
                                        "RCN 08: Endothelial" = 'Endothelial',
                                        'RCN 06: Stroma, immune excluded' = 'Stroma, immune excluded',
                                        'RCN 01: Tumor core' = 'Tumor core',
                                        'RCN 10: DC' = "DC"))

# Annotations of LDA 35 microns clusters
df <- df %>%
  mutate(lda_108_k13 = recode(lda_108_k13,
                              "RCN 08: Endothelial" = 'Endothelial',
                              'RCN 06: Stroma, immune excluded' = 'Stroma, immune excluded',
                              'RCN 01: Tumor core' = 'Tumor core',
                              "RCN 09: DC" = 'DC'))

df <- df %>%
  mutate(spatial_count_40_k15 = recode(spatial_count_40_k15,
                                       "RCN 08: Endothelial" = 'Endothelial',
                                       'RCN 05: Stroma, immune excluded' = 'Stroma, immune excluded',
                                       'RCN 01: Tumor core' = 'Tumor core',
                                       'RCN 09: DC' = 'DC'))

# Annotations of LDA 13 microns clusters
df <- df %>%
  mutate(lda_40_k13 = recode(lda_40_k13,
                             "RCN 07: Endothelial" = 'Endothelial',
                             'RCN 04: Stroma, immune excluded' = 'Stroma, immune excluded',
                             'RCN 01: Tumor core' = 'Tumor core',
                             'RCN 08: DC' = "DC"))

df <- df %>%
  mutate(lda_77_k15 = recode(lda_77_k15,
                              "RCN 09: Endothelial" = 'Endothelial',
                              'RCN 07: Stroma, immune excluded' = 'Stroma, immune excluded',
                              'RCN 01: Tumor core' = 'Tumor core',
                              'RCN 10: DC' = "DC"))

df <- df %>%
  mutate(spatial_count_77_k12 = recode(spatial_count_77_k12,
                                       "RCN 07: Endothelial" = 'Endothelial',
                                       'RCN 05: Stroma, immune excluded' = 'Stroma, immune excluded',
                                       'RCN 01: Tumor core' = 'Tumor core',
                                       'RCN 09: DC' = "DC"))

#_______________________________________________________________________________

# Define the function to calculate proportions
calculate_proportions <- function(df, cluster_col, phenotype_col) {
  
  # Count the phenotypes within each cluster
  phenotype_counts <- df %>% 
    count(!!sym(cluster_col), !!sym(phenotype_col))
  
  # Spread the counts into a wide format (each phenotype gets its own column)
  phenotype_counts_wide <- phenotype_counts %>%
    spread(key = !!sym(phenotype_col), value = n, fill = 0)
  
  # Normalize the counts to proportions by row
  phenotype_data <- phenotype_counts_wide[ , -1]  # Exclude the first column (cluster_lda_10)
  row_sums <- rowSums(phenotype_data, na.rm = TRUE)  # Compute row sums
  
  # Normalize each column by the row sums
  proportions <- phenotype_data %>%
    mutate(across(everything(), ~ . / row_sums))
  
  # Add the sample names (the first column from phenotype_counts_wide)
  proportions$sample <- phenotype_counts_wide[[1]]  # Add sample names back
  
  return(proportions)
}

# Usage:
cluster_col <- "lda_108_k13"  # Specify your cluster column
phenotype_col <- "phenotype"    # Specify your phenotype column

# Assuming your data frame is named `df`, call the function
proportions_kmeans <- calculate_proportions(df, cluster_col, phenotype_col)

# Simpson's Diversity Index (1 - Simpson's Index)
simpsons_index <- apply(proportions_kmeans[ , -ncol(proportions_kmeans)], 1, function(row) {
  diversity(row, index = "simpson")  # Simpson's index from vegan
})

# Shannon's Entropy Index
shannons_index <- apply(proportions_kmeans[ , -ncol(proportions_kmeans)], 1, function(row) {
  diversity(row, index = "shannon")  # Shannon's index from vegan
})

# Combine the sample names with the diversity indices
result_lda_35 <- data.frame(
  sample = proportions_kmeans$sample,  # Sample names
  simpsons_index = simpsons_index,
  shannons_index = shannons_index
)


#_______________________________________________________________________________
# Plotting diversity indexes

result_count_13$method <- "Spatial count"
result_lda_13$method <- "LDA"
result_count_35$method <- "Spatial count"
result_lda_35$method <- "LDA"
result_count_25$method <- "Spatial count"
result_lda_25$method <- "LDA"

result_count_13$radius <- "13"
result_lda_13$radius <- "13"
result_count_25$radius <- "25"
result_lda_25$radius <- "25"
result_count_35$radius <- "35"
result_lda_35$radius <- "35"

shannon_df <- bind_rows(result_count_13, result_lda_13, result_count_25, result_lda_25, result_count_35, result_lda_35)

shannon_df <- shannon_df %>% filter((sample %in% c('Stroma, immune excluded')))

# Shannon
ggplot(shannon_df, aes(x=sample, y=shannons_index, color=shannons_index)) +
  geom_point(aes(shape = method), size = 4) +
  scale_color_gradient2(low = "#0b71b0", mid = "gray90", high = "#cc2127", midpoint = 1) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cluster", y = 'Shannon index', color = 'Shannon index', shape = 'Method') +
  scale_shape_manual(values = c(16, 17))

# Simpson
ggplot(shannon_df, aes(x=sample, y=simpsons_index, color=simpsons_index)) +
  geom_point(aes(shape = method), size = 4) +
  scale_color_gradient2(low = "#0b71b0", mid = "gray90", high = "#cc2127", midpoint = 0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "Cluster", y = 'Simpson index', color = 'Simpson index', shape = 'Method') +
  scale_shape_manual(values = c(16, 17))

# Combined plot

# Combined plot: Shannon index on the y-axis, but color scale based on Simpson index
gg <- ggplot(shannon_df, aes(x = radius, y = shannons_index, color = simpsons_index)) +
  geom_point(aes(shape = method), size = 5) +
  scale_color_gradient2(low = "#0b71b0", mid = "gray90", high = "#cc2127", midpoint = 0.5, limits = c(0.0, 0.9)) +  # Use Simpson index for the color scale
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1)) +
  labs(x = "Radius", y = 'Shannon index', color = 'Simpson index', shape = 'Method') +
  scale_shape_manual(values = c(16, 17))  # Dot for LDA, Triangle for Spatial count

