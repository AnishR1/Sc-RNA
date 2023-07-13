library(Seurat)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(future)
library(future.apply)

#options(future.globals.maxSize = 1000 * 1024)
#plan(multicore)
setwd("C:/Users/anish/OneDrive/Desktop/Masters/sc-RNA")
dirs <- list.dirs(path = "C:/Users/anish/OneDrive/Desktop/Masters/sc-RNA/data",recursive = F, full.names = T)
for (i in dirs){
  
  x <- ReadMtx(mtx = paste0(i, "/matrix.mtx.gz"),
          features = paste0(i,"/features.tsv.gz"),
          cells = paste0(i,"/barcodes.tsv.gz"))
  
  name1 <- gsub("_filtered_feature_bc_matrix","",i)
  name <- gsub("C:/Users/anish/OneDrive/Desktop/Masters/sc-RNA/data/","", name1)
  assign(name,CreateSeuratObject(counts = x))
}

final_matrix <- merge(HB17_background,y = c(HB17_PDX,HB17_tumor,HB30_PDX,HB30_tumor,HB53_background,HB53_tumor), 
                      add.cell.ids = ls()[2:8],project = "HB")

final_matrix
View(final_matrix@meta.data)

final_matrix$sample <- rownames(final_matrix@meta.data)
View(final_matrix@meta.data)
final_matrix@meta.data <- separate(final_matrix@meta.data,col = "sample", into = c("Patient", "Type","Barcode"),
                         sep = "_")
View(final_matrix@meta.data)

final_matrix$mitoPercent <- PercentageFeatureSet(final_matrix,pattern = "^MT-")
View(final_matrix@meta.data)

filtered_matrix <- subset(final_matrix, subset = nCount_RNA > 800 &
                            nFeature_RNA > 500 &
                            mitoPercent < 10)

filtered_matrix <- NormalizeData(filtered_matrix)

