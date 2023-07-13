library(harmony)
library(Seurat)
library(SeuratData)
library(tidyverse)
library(ggplot2)


AvailableData()
InstallData("ifnb")
LoadData("ifnb")

ifnb$mitoPercent <- PercentageFeatureSet(ifnb, pattern = "^-MT")
View(ifnb@meta.data)

filtered_ifnb <- subset(ifnb, subset = nCount_RNA > 800 & nFeature_RNA >200
                        & mitoPercent < 5)
filtered_ifnb <- NormalizeData(filtered_ifnb)
filtered_ifnb <- FindVariableFeatures(filtered_ifnb)
filtered_ifnb <- ScaleData(filtered_ifnb)
filtered_ifnb <- RunPCA(filtered_ifnb)
ElbowPlot(filtered_ifnb)

filtered_ifnb <- RunUMAP(filtered_ifnb,dims = 1:20, reduction = "pca")
before<-DimPlot(filtered_ifnb,group.by = "stim", reduction = "umap")

harmony_ifnb <- RunHarmony(filtered_ifnb,group.by.vars = "stim", 
                           plot_convergence = FALSE)

harmony_embed <- Embeddings(harmony_ifnb, reduction = "harmony")

harmony_ifnb <- harmony_ifnb %>% RunUMAP(reduction = "harmony", dims = 1:20) %>%
  FindNeighbors(reduction = "harmony",dims = 1:20) %>%
  FindClusters(resolution = 0.5)

after <- DimPlot(harmony_ifnb,reduction = "umap",group.by = "stim")

before|after
#setwd("C:/Users/anish/OneDrive/Desktop/Masters/sc-RNA")

write_rds(harmony_ifnb,file = "harmony_ifnb")
