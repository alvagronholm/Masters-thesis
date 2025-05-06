# Code for Dotplots of Log2 Fold Changes in RCNs by PFS Group (Spatial count & LDA)

# Author: Alva Grönholm
# Creation date: 1.3.2025
# Modification date: 3.5.2025

# Purpose
# This code computes dotplots of the relative abundance of annotated RCNs between short and long PFS groups 
# across spatial scales, using log2 fold change based on proportional representation within each group.
# Spatial count and LDA is plotted separately 

#_______________________________________________________________________________

library(tidyverse)

#_______________________________________________________________________________
# Spatial count

# Annotations of spatial_count 40 microns clusters
df <- df %>%
  mutate(spatial_count_40_k15 = recode(spatial_count_40_k15,
                                       "RCN 04: TSI - DC" = "TSI - DC",
                                       "RCN 01: Tumor core" = "Tumor core",
                                       "RCN 05: Stroma, immune excluded" = "Stroma, immune excluded",
                                       'RCN 06: Stroma 1' = 'Stroma 1',
                                       'RCN 08: Endothelial' = 'Endothelial',
                                       'RCN 03: TSI - stroma' = 'TSI - stroma',
                                       'RCN 09: DC' = 'DC',
                                       'RCN 10: CD163+ DC' = 'CD163+ DC',
                                       'RCN 13: CD8+ T cells' = 'CD8+ T cells',
                                       'RCN 02: TSI - tumor' = 'TSI - tumor',
                                       'RCN 11: CD163- Macrophages' = 'CD163- Macrophages',
                                       'RCN 12: CD163+ Macrophages' = 'CD163+ Macrophages',
                                       'RCN 14: CD4+ T cells' = 'CD4+ T cells',
                                       'RCN 07: Stroma 2' = 'Stroma 2',
                                       'RCN 15: B-cells' = 'B-cells'))

# Annotations of spatial_count 25 microns clusters
df <- df %>%
  mutate(spatial_count_77_k12 = recode(spatial_count_77_k12,
                                       "RCN 01: Tumor core" = "Tumor core",
                                       "RCN 06: Stroma, immune infiltrated" = "Stroma, immune infiltrated",
                                       "RCN 03: TSI - stroma" = "TSI - stroma",
                                       'RCN 09: DC' = 'DC',
                                       'RCN 04: TSI - DC' = 'TSI - DC',
                                       'RCN 10: CD163+ DC' = 'CD163+ DC',
                                       'RCN 11: CD8+ T cells' = 'CD8+ T cells',
                                       'RCN 12: B-cells' = 'B-cells',
                                       'RCN 07: Endothelial' = 'Endothelial',
                                       'RCN 02: TSI - tumor' = 'TSI - tumor',
                                       'RCN 08: Immune compartment' = 'Immune compartment',
                                       'RCN 05: Stroma, immune excluded' = 'Stroma, immune excluded'))

# Annotations of spatial_count 35 microns clusters
df <- df %>%
  mutate(spatial_count_108_k13 = recode(spatial_count_108_k13,
                                        "RCN 12: CD8+ T cells" = "CD8+ T cells",
                                        "RCN 05: TSI - DC" = "TSI - DC",
                                        "RCN 07: Stroma, immune infiltrated" = 'Stroma, immune infiltrated',
                                        'RCN 01: Tumor core' = 'Tumor core',
                                        'RCN 04: TSI - stroma' = 'TSI - stroma',
                                        'RCN 09: Immune compartment' = 'Immune compartment',
                                        'RCN 11: CD163+ DC' = 'CD163+ DC',
                                        'RCN 10: DC' = 'DC',
                                        'RCN 13: B-cells' = 'B-cells',
                                        'RCN 03: TSI - tumor 2' = 'TSI - tumor 2',
                                        'RCN 02: TSI - tumor 1' = 'TSI - tumor 1',
                                        'RCN 08: Endothelial' = 'Endothelial',
                                        'RCN 06: Stroma, immune excluded' = 'Stroma, immune excluded'))

# Filter and reshape data
df_filtered <- df %>%
  filter(PFS_33p != "Middle class") %>% # Exclude "Middle class"
  select(CellID, imageid, PFS_33p, spatial_count_40_k15, spatial_count_77_k12, spatial_count_108_k13) %>%
  pivot_longer(cols = starts_with("spatial_count"), names_to = "radius", values_to = "RCN")

# Calculate proportions
proportions <- df_filtered %>%
  group_by(PFS_33p, RCN, radius) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(PFS_33p, radius) %>%
  mutate(total = sum(count), proportion = count / total)

# Reshape wide to calculate fold change (FC)
proportions_wide <- proportions %>%
  select(PFS_33p, RCN, radius, proportion) %>%
  pivot_wider(names_from = PFS_33p, values_from = proportion)

# Compute log2 FC (Long/Short)
fc_results <- proportions_wide %>%
  mutate(FC = Long / Short, log2_FC = log2(FC))

dotplot_spatialcount <- ggplot(fc_results, aes(x = log2_FC, y = RCN, color = radius)) +
  geom_point(size = 5) +
  scale_color_manual(values = c("#C39BD3", "#AED581", "#FFB877")) +  # Customize colors for radii
  labs(x = "Fold Change", color = "Radius") +
  theme_minimal() +
  scale_x_continuous(limits = c(-4.6, 4), breaks = seq(-4, 4, by = 1)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme(
    panel.grid.major = element_line(color = "grey50"),  # Change major grid color
    panel.grid.minor = element_line(color = "grey50")   # Change minor grid color
  )


#_______________________________________________________________________________
# LDA

# Annotations of LDA 13 microns clusters
df <- df %>%
  mutate(lda_40_k13 = recode(lda_40_k13,
                             "RCN 04: Stroma, immune excluded" = "Stroma, immune excluded",
                             "RCN 02: TSI - tumor" = "TSI - tumor",
                             "RCN 08: DC" = "DC",
                             'RCN 11: B-cells and Macrophages' = 'B-cells and Macrophages',
                             'RCN 09: CD163+ DC' = 'CD163+ DC',
                             'RCN 05: Stroma 1' = 'Stroma 1',
                             'RCN 12: CD8+ T cells' = 'CD8+ T cells',
                             'RCN 06: Stroma - CD4+ T cells' = 'Stroma - CD4+ T cells',
                             'RCN 01: Tumor core' = 'Tumor core',
                             'RCN 10: CD163+ Macrophages' = 'CD163+ Macrophages',
                             'RCN 07: Endothelial' = 'Endothelial',
                             'RCN 13: CD4+ T cells' = 'CD4+ T cells',
                             'RCN 03: TSI - DC' = 'TSI - DC'))

# Annotations of LDA 25 microns clusters
df <- df %>%
  mutate(lda_77_k15 = recode(lda_77_k15,
                             "RCN 13: CD4+ T cells" = "CD4+ T cells",
                             "RCN 05: TSI - CD163+ DC" = "TSI - CD163+ DC",
                             "RCN 08: Stroma, immune infiltrated" = "Stroma, immune infiltrated",
                             'RCN 09: Endothelial' = 'Endothelial',
                             'RCN 06: TSI - macrophages' = 'TSI - macrophages',
                             'RCN 15: B-cells' = 'B-cells',
                             'RCN 10: DC' = 'DC',
                             'RCN 02: TSI - tumor' = 'TSI - tumor',
                             'RCN 14: CD8+ T cells' = 'CD8+ T cells',
                             'RCN 03: TSI - stroma' = 'TSI - stroma',
                             'RCN 07: Stroma, immune excluded' = 'Stroma, immune excluded',
                             'RCN 12: Macrophages' = 'Macrophages',
                             'RCN 04: TSI - DC' = 'TSI - DC',
                             'RCN 01: Tumor core' = 'Tumor core',
                             'RCN 11: CD163+ DC' = 'CD163+ DC'))

# Annotations of LDA 35 microns clusters
df <- df %>%
  mutate(lda_108_k13 = recode(lda_108_k13,
                              "RCN 11: Macrophages" = "Macrophages",
                              "RCN 02: TSI 1" = "TSI 1",
                              "RCN 08: Endothelial" = "Endothelial",
                              'RCN 06: Stroma, immune excluded' = 'Stroma, immune excluded',
                              'RCN 10: B-cells and CD163+ DC' = 'B-cells and CD163+ DC',
                              'RCN 03: TSI 2' = 'TSI 2',
                              'RCN 12: CD8+ T cells' = 'CD8+ T cells',
                              'RCN 05: TSI - myeloid' = 'TSI - myeloid',
                              'RCN 01: Tumor core' = 'Tumor core',
                              'RCN 09: DC' = 'DC',
                              'RCN 07: Stroma, immune infiltrated' = 'Stroma, immune infiltrated',
                              'RCN 13: CD4+ T cells' = 'CD4+ T cells',
                              'RCN 04: TSI - stroma' = 'TSI - stroma'))



# Filter and reshape data
df_filtered <- df %>%
  filter(PFS_33p != "Middle class") %>% # Exclude "Middle class"
  select(CellID, imageid, PFS_33p, lda_40_k13, lda_77_k15, lda_108_k13) %>%
  pivot_longer(cols = starts_with("lda"), names_to = "radius", values_to = "RCN")

# Calculate proportions
proportions <- df_filtered %>%
  group_by(PFS_33p, RCN, radius) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(PFS_33p, radius) %>%
  mutate(total = sum(count), proportion = count / total)

# Reshape wide to calculate fold change (FC)
proportions_wide <- proportions %>%
  select(PFS_33p, RCN, radius, proportion) %>%
  pivot_wider(names_from = PFS_33p, values_from = proportion)

# Compute log2 FC (Long/Short)
fc_results <- proportions_wide %>%
  mutate(FC = Long / Short, log2_FC = log2(FC))


dotplot <- ggplot(fc_results, aes(x = log2_FC, y = RCN, color = radius)) +
  geom_point(size = 5) +
  scale_color_manual(values = c("#C39BD3", "#AED581", "#FFB877")) +  # Customize colors for radii
  labs(x = "Fold Change", color = "Radius") +
  theme_minimal() +
  scale_x_continuous(limits = c(-4.6, 4), breaks = seq(-4, 4, by = 1)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme(
    panel.grid.major = element_line(color = "grey50"),  # Change major grid color
    panel.grid.minor = element_line(color = "grey50")   # Change minor grid color
  )
