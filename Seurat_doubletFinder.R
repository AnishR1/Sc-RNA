library(DoubletFinder)
library(Seurat)
library(tidyverse)
library(ggplot2)

counts <- ReadMtx(mtx = "C:/Users/anish/OneDrive/Desktop/Masters/sc-RNA/doublet_data/10k_PBMC_3p_nextgem_Chromium_Controller_raw_feature_bc_matrix/raw_feature_bc_matrix/matrix.mtx.gz",
        features = "C:/Users/anish/OneDrive/Desktop/Masters/sc-RNA/doublet_data/10k_PBMC_3p_nextgem_Chromium_Controller_raw_feature_bc_matrix/raw_feature_bc_matrix/features.tsv.gz",
        cells = "C:/Users/anish/OneDrive/Desktop/Masters/sc-RNA/doublet_data/10k_PBMC_3p_nextgem_Chromium_Controller_raw_feature_bc_matrix/raw_feature_bc_matrix/barcodes.tsv.gz")

ct_seurat <- CreateSeuratObject(counts = counts)
ct_seurat$mitoPercent <- PercentageFeatureSet(ct_seurat,pattern = "^MT-")
filtered_seurat <- subset(ct_seurat,subset = nFeature_RNA > 500 & nCount_RNA > 800
                          & mitoPercent < 10)
filtered_seurat <- NormalizeData(object = filtered_seurat)
filtered_seurat <- FindVariableFeatures(object = filtered_seurat)
filtered_seurat <- ScaleData(object = filtered_seurat)
filtered_seurat <- RunPCA(object = filtered_seurat)
ElbowPlot(filtered_seurat)
filtered_seurat <- FindNeighbors(object = filtered_seurat, dims = 1:20)
filtered_seurat <- FindClusters(object = filtered_seurat)
filtered_seurat <- RunUMAP(object = filtered_seurat,dims = 1:20)

sweep.Res <- paramSweep_v3(filtered_seurat,PCs = 1:20,sct = FALSE)
sweep.stats <- summarizeSweep(sweep.Res,GT = FALSE)
pk_data <- find.pK(sweep.stats = sweep.stats)

ggplot(pk_data,aes(pK,BCmetric,group = 1))+ geom_point() + geom_line()

pK <- pk_data %>% filter(BCmetric == max(BCmetric)) %>% select(pK)
pK <- as.numeric(as.character(pK[[1]]))

annotations <- filtered_seurat@meta.data$seurat_clusters
homotypic_prop <- modelHomotypic(annotations = annotations)
nExp_poi <- round(0.076*nrow(filtered_seurat@meta.data))
nExp_adj <- round(nExp_poi*(1-homotypic_prop))

filtered_seurat <- doubletFinder_v3(filtered_seurat,
                                    PCs = 1:20,
                                    pK = pK, pN = 0.25,
                                    nExp = nExp_adj,
                                    reuse.pANN = FALSE, sct = FALSE)

DimPlot(filtered_seurat,reduction = "umap", group.by = "DF.classifications_0.25_0.2_691")
