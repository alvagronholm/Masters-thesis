Code for master's thesis: Profiling the Spatial Tumor Microenvironment Architecture in Highly Multiplexed Images of High-Grade Serous Ovarian Carcinoma

1. Composition Analysis

Files: composition_analysis.R

Purpose: Visualizes global and immune cell type compositions across patients.

Output: Heatmaps and barplots annotated with PFS status.

2. Spatial Count & LDA + Clustering Pipeline

Files: compute_spatial_count_and_LDA.ipynb

Purpose:
- Compute spatial counts and LDA.
- Optimize KMeans clusters via elbow curves and KneeLocator.
- Visualize spatial RCNs and export clustered data.
- Compare RCNs across spatial scales and modalities.

Output:
- Clustered RCN CSV files.
- RCN visualizations across methods and radii.

3.  RCN Analysis

Files: spatial_count_LDA_analysis.R

Purpose:
- Generate barplots and heatmaps of phenotype distributions in RCNs.
- Perform group-wise comparisons using boxplots and rank-abundance curves.

Output:
- Proportional phenotype heatmaps.
- Z-scored patient-level heatmaps.
- Clinical group comparisons.

4. RCN Diversity Analysis

Files: spatial_analysis_diversity.R

Purpose: Evaluates phenotypic diversity (Shannon & Simpson indices) within RCNs.

Output: Scatter plots comparing diversity across spatial scales and analysis methods.

5. Clustree Analysis

Files: clustree.R

Purpose: Computes clustree plots for RCNs derived from spatial count and LDA methods across 13, 25, and 35 µm radii.

Output: Clustree visualizations to assess RCN changes across spatial scales.

6. RCN Abundances across Clinical Groups (Dotplots)

Files: dotplot_RCN_clinical_groups.R

Purpose: Computes log2 fold change in RCN proportions between short and long PFS groups.

Output: Dotplots for spatial count and LDA-derived RCNs across scales.
