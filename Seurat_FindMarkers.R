# Load libraries
library(Seurat)
library(tidyverse)

#Read the data
harmony_data <- readRDS("C:/Users/anish/OneDrive/Desktop/Masters/sc-RNA/harmony_ifnb")
View(harmony_data@meta.data)

# Visualize the read data by grouping them according to condition and clusters

clusters <- DimPlot(harmony_data,reduction = "umap",group.by = "seurat_clusters",label = TRUE)
conditions <- DimPlot(harmony_data,reduction = "umap",group.by = "stim")

conditions|clusters

# Set Default Assay to RNA

DefaultAssay(harmony_data) <- "RNA"
# Finding differentially expressed markers in different conditions
markers3 <- FindConservedMarkers(harmony_data,ident.1 = 3, grouping.var = "stim")
head(markers3)

# Visualize the markers by plotting them
FeaturePlot(harmony_data,features = c("FCGR3A"),min.cutoff = "q10")

# Renaming the cells in cluster 3 (It is known in this case because of the dataset)
Idents(harmony_data)
harmony_data <- RenameIdents(harmony_data,"CD16 Monocytes" = "CD16 Mono")

DimPlot(harmony_data,reduction = "umap",label = T)

# Setting the idents of the data to seurat annotations column
Idents(harmony_data) <- harmony_data@meta.data$seurat_annotations
Idents(harmony_data)

# Visualize the plot with labels
DimPlot(harmony_data,reduction = "umap",label = T)

## To find markers before and after IFNB treatment
# Create a new column 

harmony_data$celltype_cond <- paste0(harmony_data$seurat_annotations,"_",harmony_data$stim)
View(harmony_data@meta.data)
Idents(harmony_data) <- harmony_data$celltype_cond

DimPlot(harmony_data,reduction = "umap",label = T)

# Comparison of CD16 Monocytes control group and stimulated group
interferon_response <- FindMarkers(harmony_data,ident.1 = "CD16 Mono_CTRL", ident.2 = "CD16 Mono_STIM")
head(interferon_response)

FeaturePlot(harmony_data,features = c("FCGR3A","AIF1","IFIT1"),split.by = "stim",min.cutoff = "q10")
