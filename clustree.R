# Code for clustree analysis (spatial count and LDA)

# Author: Alva Grönholm
# Creation date: 1.2.2025
# Modification date: 3.5.2025

# Purpose
# This code computes clustree plots for spatial count and LDA across 13, 25 and 35 micron radii.

#_______________________________________________________________________________
library(clustree)
library(ggplot2)
library(mclust)  
library(dplyr)   

#_______________________________________________________________________________

# SPATIAL COUNT

rcn_data_spatial <- df %>% select(spatial_count_40_k15, spatial_count_77_k12, spatial_count_108_k13, imageid, CellID)

rcn_data_spatial <- rcn_data_spatial %>% dplyr::rename(cluster_1 = spatial_count_40_k15, cluster_2 = spatial_count_77_k12, cluster_3 = spatial_count_108_k13)

rcn_data_spatial <- rcn_data_spatial %>% mutate(cluster_1 = recode(cluster_1,
                                                                   "RCN 01: Tumor core" = 1,
                                                                   "RCN 02: TSI - tumor" = 2,
                                                                   "RCN 03: TSI - stroma" = 3,
                                                                   "RCN 04: TSI - DC" = 4,
                                                                   "RCN 05: Stroma, immune excluded" = 5,
                                                                   "RCN 06: Stroma 1" = 6,
                                                                   "RCN 07: Stroma 2" = 7,
                                                                   "RCN 08: Endothelial" = 8,
                                                                   "RCN 09: DC" = 9,
                                                                   "RCN 10: CD163+ DC" = 10,
                                                                   "RCN 11: CD163- Macrophages" = 11,
                                                                   "RCN 12: CD163+ Macrophages" = 12,
                                                                   "RCN 13: CD8+ T cells" = 13,
                                                                   "RCN 14: CD4+ T cells" = 14,
                                                                   "RCN 15: B-cells" = 15))


rcn_data_spatial <- rcn_data_spatial %>% mutate(cluster_2 = recode(cluster_2,
                                                                   "RCN 01: Tumor core" = 1,
                                                                   "RCN 02: TSI - tumor" = 2,
                                                                   "RCN 03: TSI - stroma" = 3,
                                                                   "RCN 04: TSI - DC" = 4,
                                                                   "RCN 05: Stroma, immune excluded" = 5,
                                                                   "RCN 06: Stroma, immune infiltrated" = 6,
                                                                   "RCN 07: Endothelial" = 7,
                                                                   "RCN 08: Immune compartment" = 8,
                                                                   "RCN 09: DC" = 9,
                                                                   "RCN 10: CD163+ DC" = 10,
                                                                   "RCN 11: CD8+ T cells" = 11,
                                                                   "RCN 12: B-cells" = 12))

rcn_data_spatial <- rcn_data_spatial %>% mutate(cluster_3 = recode(cluster_3,
                                                                   "RCN 01: Tumor core" = 1,
                                                                   "RCN 02: TSI - tumor 1" = 2,
                                                                   "RCN 03: TSI - tumor 2" = 3,
                                                                   "RCN 04: TSI - stroma" = 4,
                                                                   "RCN 05: TSI - DC" = 5,
                                                                   "RCN 06: Stroma, immune excluded" = 6,
                                                                   "RCN 07: Stroma, immune infiltrated" = 7,
                                                                   "RCN 08: Endothelial" = 8,
                                                                   "RCN 09: Immune compartment" = 9,
                                                                   "RCN 10: DC" = 10,
                                                                   "RCN 11: CD163+ DC" = 11,
                                                                   "RCN 12: CD8+ T cells" = 12,
                                                                   "RCN 13: B-cells" = 13))



# Run Clustree using the correct prefix
a <- clustree(rcn_data_spatial, prefix = "cluster_", prop_filter = 0.1) +
  scale_edge_color_continuous(low = '#0076A8', high = '#d62728') +
  theme_minimal() +
  labs(title = "Spatial count radius 13 microns -> 25 microns -> 35 microns") 

# Save Clustree plot directly
ggsave("F:/Alva/Final_spatial_analysis_plots/clustree_plot_spatial_count.pdf", plot = a, width = 9, height = 7, device = "pdf")

#_______________________________________________________________________________

# LDA

rcn_data_lda <- df %>% select(lda_40_k13, lda_77_k15, lda_108_k13, imageid, CellID)

rcn_data_lda <- rcn_data_lda %>% dplyr::rename(cluster_1 = lda_40_k13, cluster_2 = lda_77_k15, cluster_3 = lda_108_k13)

rcn_data_lda <- rcn_data_lda %>% mutate(cluster_1 = recode(cluster_1,
                                                   "RCN 01: Tumor core" = 1,
                                                   "RCN 02: TSI - tumor" = 2,
                                                   "RCN 03: TSI - DC" = 3,
                                                   "RCN 04: Stroma, immune excluded" = 4,
                                                   "RCN 05: Stroma 1" = 5,
                                                   "RCN 06: Stroma - CD4+ T cells" = 6,
                                                   "RCN 07: Endothelial" = 7,
                                                   "RCN 08: DC" = 8,
                                                   "RCN 09: CD163+ DC" = 9,
                                                   "RCN 10: CD163+ Macrophages" = 10,
                                                   "RCN 11: B-cells and Macrophages" = 11,
                                                   "RCN 12: CD8+ T cells" = 12,
                                                   "RCN 13: CD4+ T cells" = 13))


rcn_data_lda <- rcn_data_lda %>% mutate(cluster_2 = recode(cluster_2,
                                                           "RCN 01: Tumor core" = 1,
                                                           "RCN 02: TSI - tumor" = 2,
                                                           "RCN 03: TSI - stroma" = 3,
                                                           "RCN 04: TSI - DC" = 4,
                                                           "RCN 05: TSI - CD163+ DC" = 5,
                                                           "RCN 06: TSI - macrophages" = 6,
                                                           "RCN 07: Stroma, immune excluded" = 7,
                                                           "RCN 08: Stroma, immune infiltrated" = 8,
                                                           "RCN 09: Endothelial" = 9,
                                                           "RCN 10: DC" = 10,
                                                           "RCN 11: CD163+ DC" = 11,
                                                           "RCN 12: Macrophages" = 12,
                                                           "RCN 13: CD4+ T cells" = 13,
                                                           "RCN 14: CD8+ T cells" = 14,
                                                           "RCN 15: B-cells" = 15))

rcn_data_lda <- rcn_data_lda %>% mutate(cluster_3 = recode(cluster_3,
                                                           "RCN 01: Tumor core" = 1,
                                                           "RCN 02: TSI 1" = 2,
                                                           "RCN 03: TSI 2" = 3,
                                                           "RCN 04: TSI - stroma" = 4,
                                                           "RCN 05: TSI - myeloid" = 5,
                                                           "RCN 06: Stroma, immune excluded" = 6,
                                                           "RCN 07: Stroma, immune infiltrated" = 7,
                                                           "RCN 08: Endothelial" = 8,
                                                           "RCN 09: DC" = 9,
                                                           "RCN 10: B-cells and CD163+ DC" = 10,
                                                           "RCN 11: Macrophages" = 11,
                                                           "RCN 12: CD8+ T cells" = 12,
                                                           "RCN 13: CD4+ T cells" = 13))



# Run Clustree using the correct prefix
b <- clustree(rcn_data_lda, prefix = "cluster_", prop_filter = 0.1) +
  scale_edge_color_continuous(low = '#0076A8', high = '#d62728') +
  theme_minimal() +
  labs(title = "LDA radius 13 microns -> 25 microns -> 35 microns")
