# Load required libraries
library(Seurat)
library(tidyverse)

#Read the file and store count matrix
sparseMatrix <- Read10X_h5(filename = "C:/Users/anish/Downloads/20k_NSCLC_DTC_3p_nextgem_Multiplex_count_raw_feature_bc_matrix.h5")
#str(sparseMatrix)
counts <- sparseMatrix$`Gene Expression`

#Create seurat object
counts_seurat <- CreateSeuratObject(counts = counts, min.cells = 3, min.features = 200, project = "NSCLC")
# str(counts_seurat)
#Calculate mitochondrial gene percent 
counts_seurat[["percent_mt"]] <- PercentageFeatureSet(counts_seurat, pattern = "^MT")
View(counts_seurat@meta.data)

# Plot count matrix to check distribution of data
VlnPlot(counts_seurat,features = c("nFeature_RNA","nCount_RNA","percent_mt"),ncol = 3)
FeatureScatter(counts_seurat,feature1 = "nFeature_RNA", feature2 = "nCount_RNA") + geom_smooth(method = "lm")

# Filter the genes with by setting threshold for copy numbers and mitochondrial gene percent
counts_seurat <- subset(counts_seurat, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                          percent_mt < 5)
counts_seurat

# Run the Normalization
counts_seurat <- NormalizeData(counts_seurat)
# Find variability in data
counts_seurat <- FindVariableFeatures(counts_seurat,selection.method = "vst",nfeatures = 2000)
top10 <- head(VariableFeatures(counts_seurat),10)
top10
plot1 <- VariableFeaturePlot(counts_seurat)
LabelPoints(plot1,points = top10, repel = TRUE)

all_genes <- rownames(counts_seurat)

# Scale data
counts_seurat <- ScaleData(counts_seurat,features = all_genes)
# Run PCA
counts_seurat <- RunPCA(counts_seurat, features = VariableFeatures(counts_seurat))
print(counts_seurat[["pca"]],dims = 1:5, nfeatures = 10)
DimHeatmap(counts_seurat,dims = 1, cells = 5000, balanced = TRUE)
ElbowPlot(counts_seurat)

# Find Neighbours and assign clusters
counts_seurat <- FindNeighbors(counts_seurat, dims = 1:15)
counts_seurat <- FindClusters(counts_seurat, resolution = c(0.1,0.3,0.5,0.7,1))

DimPlot(counts_seurat,group.by = "RNA_snn_res.0.1", label = TRUE)

# Set idents of count matrix
Idents(counts_seurat)
Idents(counts_seurat) <- "RNA_snn_res.0.1"
Idents(counts_seurat)

# Run UMAP
counts_seurat <- RunUMAP(counts_seurat,dims = 1:15)
DimPlot(counts_seurat, reduction = "umap", label = TRUE)
