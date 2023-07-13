library(Seurat)
library(tidyverse)

harmony_data <- readRDS("C:/Users/anish/OneDrive/Desktop/Masters/sc-RNA/harmony_ifnb")
View(harmony_data@meta.data)


clusters <- DimPlot(harmony_data,reduction = "umap",group.by = "seurat_clusters",label = TRUE)
conditions <- DimPlot(harmony_data,reduction = "umap",group.by = "stim")

conditions|clusters


markers3 <- FindConservedMarkers(harmony_data,ident.1 = 3, grouping.var = "stim")
head(markers3)

FeaturePlot(harmony_data,features = c("FCGR3A"),min.cutoff = "q10")

Idents(harmony_data)
harmony_data <- RenameIdents(harmony_data,"CD16 Monocytes" = "CD16 Mono")

DimPlot(harmony_data,reduction = "umap",label = T)

Idents(harmony_data) <- harmony_data@meta.data$seurat_annotations
Idents(harmony_data)

DimPlot(harmony_data,reduction = "umap",label = T)

harmony_data$celltype_cond <- paste0(harmony_data$seurat_annotations,"_",harmony_data$stim)
View(harmony_data@meta.data)
Idents(harmony_data) <- harmony_data$celltype_cond

DimPlot(harmony_data,reduction = "umap",label = T)


interferon_response <- FindMarkers(harmony_data,ident.1 = "CD16 Mono_CTRL", ident.2 = "CD16 Mono_STIM")
head(interferon_response)

FeaturePlot(harmony_data,features = c("FCGR3A","AIF1","IFIT1"),split.by = "stim",min.cutoff = "q10")
