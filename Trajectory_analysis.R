library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(tidyverse)
library(ggplot2)

markers <- read.delim("C:/Users/anish/Downloads/ABC_Marker.txt", header = T)
metaData <- read.delim("C:/Users/anish/Downloads/ABC_Meta.txt", header = T)
expr <- read.delim("C:/Users/anish/Downloads/ABC_umi_matrix_7551_cells/ABC_umi_matrix_7551_cells.csv",
                   header = T, sep = ",")

expr_t <- t(expr)

Seu_expr <- CreateSeuratObject(counts = expr_t)

Seu_expr@meta.data<- merge(Seu_expr@meta.data,metaData,by.x = "row.names", by.y = "cell_id")
Seu_expr@meta.data <- Seu_expr@meta.data %>%
  column_to_rownames(var = "Row.names")
Seu_expr$MitoPercent <- PercentageFeatureSet(Seu_expr,pattern = "^MT-")
Seu_expr_filtered <- subset(Seu_expr,subset = nCount_RNA > 800 & nFeature_RNA > 500
                   & MitoPercent < 10)

unique(Seu_expr_filtered@meta.data$population)
Idents(Seu_expr_filtered) <- Seu_expr_filtered$population
b_cells_seurat <- subset(Seu_expr_filtered,idents = "b")

unique(b_cells_seurat@meta.data$redefined_cluster)


b_cells_seurat <- NormalizeData(b_cells_seurat)
b_cells_seurat <- FindVariableFeatures(b_cells_seurat)
b_cells_seurat <- ScaleData(b_cells_seurat)
b_cells_seurat <- RunPCA(b_cells_seurat)

b_cells_seurat <- FindNeighbors(b_cells_seurat,dims = 1:30)
b_cells_seurat <- FindClusters(b_cells_seurat, resolution = 0.9)

b_cells_seurat <- RunUMAP(b_cells_seurat, dims = 1:30, n.neighbors = 50)

p1 <- DimPlot(b_cells_seurat, reduction = "umap", group.by = "redefined_cluster",
              label = T)
p2 <- DimPlot(b_cells_seurat, reduction = "umap", group.by = "seurat_clusters",
              label = T)

p1|p2


cds <- as.cell_data_set(b_cells_seurat)
colData(cds)
fData(cds)$gene_names <- rownames(fData(cds))
counts(cds)

recreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(recreate.partition) <- cds@colData@rownames
recreate.partition <- as.factor(recreate.partition)
cds@clusters$UMAP$partitions <- recreate.partition


cluster_list <- b_cells_seurat@active.ident
cds@clusters$UMAP$clusters <- cluster_list


cds@int_colData@listData$reducedDims$UMAP <- b_cells_seurat@reductions$umap@cell.embeddings

before_traj <- plot_cells(cds,
           color_cells_by = "cluster", label_groups_by_cluster = F
           group_label_size = 5) +
  theme(legend.position = "right")

cluster_names <- plot_cells(cds,
                            color_cells_by = "redefined_cluster", label_groups_by_cluster = F,
                            group_label_size = 5) +
  scale_color_manual(values = c("red","blue","green","maroon","yellow","grey","cyan"))+
  theme(legend.position = "right")
before_traj | cluster_names

cds <- learn_graph(cds,use_partition = FALSE)

plot_cells(cds,
           color_cells_by = "redefined_cluster",
           label_cell_groups = T,
           label_roots = F,
           label_leaves = F,
           label_branch_points = F,
           group_label_size = 5)

cds <- order_cells(cds, reduction_method = "UMAP",root_cells = colnames(cds[,clusters(cds)==5]))
plot_cells(cds,
           color_cells_by = "pseudotime", label_cell_groups = T,
           label_leaves = F,
           label_branch_points = F,
           group_label_size = 5)


pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data_pseudo <- as.data.frame(colData(cds))
ggplot(data = data_pseudo,aes(monocle3_pseudotime,reorder(redefined_cluster,monocle3_pseudotime, median),fill ="redefined_cluster"))+
  geom_boxplot()

deg_bcells <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)

deg_bcells %>% arrange(q_value) %>%
  filter(status == "OK") %>%
  head()

FeaturePlot(b_cells_seurat, features = c("STMN1","CD52","HMGN2"))


b_cells_seurat$pseudotime <- pseudotime(cds)
Idents(b_cells_seurat) <- b_cells_seurat$redefined_cluster
FeaturePlot(b_cells_seurat,features = "pseudotime",label = T)
