#Load required libraries
library("Seurat")
library("SeuratDisk")

#Reading different sc-RNA file formats to SeuratData to process them
rds <- readRDS("C:/Users/anish/Downloads/ependymal_cells.rds")
H5 <- Read10X_h5("C:/Users/anish/Downloads/20k_PBMC_3p_HT_nextgem_Chromium_X_filtered_feature_bc_matrix.h5",use.names = TRUE, unique.features = TRUE)
H5[1:10,1:10]
x <- CreateSeuratObject(counts = H5)
mat_obj <- ReadMtx(features = "C:/Users/anish/Downloads/raw_feature_bc_matrix/features.tsv.gz",
        mtx = "C:/Users/anish/Downloads/raw_feature_bc_matrix/matrix.mtx.gz",
        cells = "C:/Users/anish/Downloads/raw_feature_bc_matrix/barcodes.tsv.gz")
mat_obj[1:10,1:10]

Convert("C:/Users/anish/Downloads/adata_SS2_for_download.h5ad",dest = "h5Seurat",overwrite = TRUE)
H5_Seurat_data <- LoadH5Seurat("C:/Users/anish/Downloads/adata_SS2_for_download.h5Seurat")
