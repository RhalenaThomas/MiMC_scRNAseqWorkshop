# MiMC_scRNAseqWorkshop
An overview of single cell RNA sequencing: from cell ranger to annotation in Seurat

# Workbook for MiMC workshop April 26th, 2023
Data is in a separate source.
If you wish to follow this workshop and are not part of the workshop you can find the data used here:
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE186780
Use the CellRanger outputs

# Workbook contents
1. Creating a Seurat object from CellRanger output.
2. Preparing and cleaning the data
  a) Visualize QC
  b) Filter out unwanted cells
  c) Identify and remove doublets
  d) Normalization and scale
  e) Select Variable features
3. Merging and Harmonizing samples
  a) Merge samples
  b) Use Seurat find anchors to integrate
  c) Compare merged vs integrated
4. Dimensional reduction clustering and visualization
  a) PCA and component selection
  b) UMAP
  c) Clustering and visualization
5. Cluster annotation
  a) Visualize expression of known cell type markers
  b) Find cluster markers and look them up in reference cell type library
  c) Manual cluster annotation 
  d) Decisions on merging clusters
6. Automated cluster annotation
  a) Seurat label transfer
  b) scClassify
