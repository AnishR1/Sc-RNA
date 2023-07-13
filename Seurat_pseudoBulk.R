library(Seurat)
library(ExperimentHub)
library(tidyverse)
library(DESeq2)

eh <- ExperimentHub()
query(eh,"Kang")
sce <- eh[["EH2259"]]
sce_seurat <- as.Seurat(sce,data = NULL)

sce_seurat$mitoPercent <- PercentageFeatureSet(sce_seurat,pattern = "^MT-")
filtered_seurat <- subset(sce_seurat, subset = nFeature_originalexp >200 &
         nFeature_originalexp < 2500 &
         nCount_originalexp > 800 & mitoPercent < 5 & multiplets == "singlet")

filtered_seurat <- NormalizeData(filtered_seurat)
filtered_seurat <- FindVariableFeatures(filtered_seurat)
filtered_seurat <- ScaleData(filtered_seurat)
filtered_seurat <- RunPCA(filtered_seurat)
ElbowPlot(filtered_seurat)

filtered_seurat <- RunUMAP(filtered_seurat,dims = 1:20)

cell_plot <- DimPlot(filtered_seurat,reduction = "umap", group.by = "cell",label = T)
condition <- DimPlot(filtered_seurat,reduction = "umap", group.by = "stim")

cell_plot | condition

filtered_seurat$id_condtion <- paste0(filtered_seurat$stim,filtered_seurat$ind)

DefaultAssay(filtered_seurat)
cts <- AggregateExpression(filtered_seurat,
                           group.by = c("cell","id_condtion"),
                           assays = "originalexp",
                           slot = "counts",
                           return.seurat = F)
cts <- cts$originalexp
cts.t <- as.data.frame(t(cts))

splitrows <- gsub("_.*","",rownames(cts.t))
cts.split <- split.data.frame(cts.t, f = factor(splitrows))

modified_split <- lapply(cts.split, function(x)
  {
  rownames(x) <- gsub('.*_(.*)',"\\1", rownames(x))
  t(x)
})

counts_bcell <- modified_split$`B cells`

colData <- data.frame(samples = colnames(counts_bcell) )
colData <- colData %>% 
  mutate(condition = ifelse(grepl("stim",samples),"Stimulated","Control")) %>%
  column_to_rownames(var = "samples")

dds <-  DESeqDataSetFromMatrix(countData = counts_bcell,
                               colData = colData,
                               design = ~condition)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds <- DESeq(dds)

resultsNames(dds)
res <- results(dds, name = "condition_Stimulated_vs_Control")
res
