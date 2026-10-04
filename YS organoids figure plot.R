####pre-processing single-cell RNA-seq data for YS organoids ----
###R version=3.6
library(Seurat)
library(monocle)
###construction of seurat object
BMP4_Dx_count <- Read10X(data.dir = "10x_matrix_fold_path")
BMP4_Dx <- CreateSeuratObject(counts = BMP4_Dx_count, project = "BMP4_Dx", min.cells = 3, min.features = 200)
BMP4_Dx[["percent.mt"]] <- PercentageFeatureSet(BMP4_Dx, pattern = "^MT")
VlnPlot(BMP4_Dx, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(BMP4_Dx, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(BMP4_Dx, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1
plot2
BMP4_Dx <- subset(BMP4_Dx, subset = nFeature_RNA > 200 & nFeature_RNA < 10000 & percent.mt < 5)
BMP4_Dx <- NormalizeData(BMP4_Dx)
BMP4_Dx <- FindVariableFeatures(BMP4_Dx, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(BMP4_Dx)
all.genes <- rownames(BMP4_Dx)
BMP4_Dx <- ScaleData(BMP4_Dx, features = all.genes) ##select all.genes can be visualized by dotplot 
BMP4_Dx <- RunPCA(BMP4_Dx, features = VariableFeatures(object = BMP4_Dx))
print(BMP4_Dx[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(BMP4_Dx, dims = 1:2, reduction = "pca")
DimPlot(BMP4_Dx, reduction = "pca")
BMP4_Dx <- JackStraw(BMP4_Dx, num.replicate = 100)
BMP4_Dx <- ScoreJackStraw(BMP4_Dx, dims = 1:20)
JackStrawPlot(BMP4_Dx, dims = 1:20)
ElbowPlot(BMP4_Dx)
BMP4_Dx <- FindNeighbors(BMP4_Dx, dims = 1:15)
BMP4_Dx <- FindClusters(BMP4_Dx, resolution = 1)
BMP4_Dx <- RunUMAP(BMP4_Dx, dims = 1:15,return.model = T)
BMP4_Dx <- RunTSNE(BMP4_Dx, dims = 1:15)
DimPlot(BMP4_Dx, reduction = "umap",label = T)
BMP4_Dx.markers <- FindAllMarkers(BMP4_Dx, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(BMP4_Dx.markers,"path/BMP4_Dx.markers.csv")
saveRDS(BMP4_Dx,"path/BMP4_Dx.rds")


#### fig1/figS1. Monocle3 lineage tracing of human embryonic reference ----
library(monocle3)
#cds constraction
expression_matrix <- GetAssayData(ref.endo)
cell_metadata <- ref.endo@meta.data
gene_annotation <- data.frame(rownames(ref.endo), rownames(ref.endo))
rownames(gene_annotation) <- rownames(ref.endo)
colnames(gene_annotation) <- c("GeneSymbol", "gene_short_name")
cds <- new_cell_data_set(expression_matrix,                         
                         cell_metadata = cell_metadata,                         
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds, 
                 alignment_group = "batch",
                 preprocess_method = "PCA")
cds <- reduce_dimension(cds, preprocess_method = "PCA") 
# embedding transfer
cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(ref.endo, reduction = "UMAP")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
# trajectory
cds <- cluster_cells(cds, reduction_method = "UMAP")
plot_cells(cds, color_cells_by = "partition")
cds <- learn_graph(cds,
                   verbose=T,
                   use_partition=F,
                   learn_graph_control = list( 
                     minimal_branch_len = 4 ))
cds <- order_cells(cds)
#reverse pseudo
myselect <- function(cds,select.classify,my_select){  
  cell_ids <- which(colData(cds)[,select.classify] == my_select)  
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex  
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])  
  root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names    (which.max(table(closest_vertex[cell_ids,]))))]  
  root_pr_nodes
}
cds <- order_cells(cds, root_pr_nodes=myselect(cds,select.classify = 'sub_rename_EML',my_select = "DE"))
cds <- order_cells(cds)
plot_cells(cds, color_cells_by = "pseudotime",           
           show_trajectory_graph=T) + 
  plot_cells(cds,           
             color_cells_by = "sub_rename_EML",          
             label_cell_groups=T,           
             label_leaves=FALSE,           
             label_branch_points=FALSE,           
             group_label_size=3,
             alpha = 0.8)
saveRDS(cds,"cds.rds")
saveRDS(ref.endo,"ref.endo.rds")
#reverse pseudotime
original_pseudotime <- cds@principal_graph_aux$UMAP$pseudotime
reversed_pseudotime <- max(original_pseudotime) - original_pseudotime
pData(cds)$reversed_pseudotime <- reversed_pseudotime
plot_cells(cds, 
           color_cells_by = "reversed_pseudotime",
           label_branch_points = F, 
           label_roots = F,
           label_leaves = F,
           cell_size = 0.8
)

#只保留batchCS
cds <- cds[, cds@colData$batch == "batchCS"]
cds <- choose_graph_segments(cds,
                             clear_cds = F) 
saveRDS(cds,"monocle3_batchCS.rds")
plot_cells(cds, 
           color_cells_by = "reversed_pseudotime",
           label_branch_points = F, 
           label_roots = F,
           label_leaves = F,
           cell_size = 0.8)
cds@principal_graph_aux$UMAP$pseudotime <- pData(cds)$reversed_pseudotime
cds_sub <- cds[rowData(cds)$gene_short_name %in% gene_plot, ]
plot_genes_in_pseudotime(cds_sub,
                         color_cells_by="sub_rename_EML",
                         min_expr=0.5,
                         ncol = 4)
table(cds_sub@colData$batch)
#module
pr_graph_test_res <- graph_test(cds, neighbor_graph="knn")
pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=1e-2)
#module tibble
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                cell_group=colData(cds)$sub_rename_EML )
agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Celltype ", colnames(agg_mat))
#module heatmap
pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)
#module umap
plot_cells(cds, 
           genes=gene_module_df %>% filter(module %in% c(6,8, 9, 15, 16, 17, 18, 23,5,2,1,14,4,11)),
           group_cells_by="partition",
           color_cells_by="partition",
           show_trajectory_graph=FALSE)
module_gene <- tibble::tibble(gene = rowData(cds)$gene_short_name , 
                              cell_group=colData(cds)$sub_rename_EML )

### fig3. reference projection ----
ggplot()+geom_point(predict_out$umap %>% filter(pj!="query"),mapping=aes(x=UMAP_1,y=UMAP_2),color="grey")+geom_point(predict_out$umap %>% filter(pj=="query"),mapping=aes(x=UMAP_1,y=UMAP_2,color="red"))+theme_classic()
predict_out_list <- list(predict_out_day2,
                         predict_out_day4,
                         predict_out_day6,
                         predict_out_day8,
                         predict_out_day10,
                         predict_out_day15)




alluvial_data_list <- lapply(predict_out_list, function(x){
  df <- x$full.anno[,1:3]
  
  alluvial_data <- df %>%
    mutate(percentage = n / sum(n) * 100) %>%
    ungroup()
  alluvial_data
})
# cell type percentage pie plot
library(ggplot2)
library(dplyr)
df <- predict_out_day6$full.anno[,1:3]
celltype_counts <- df %>%
  count(sub_pred_EML, name = "count") %>%
  mutate(percentage = count / sum(count) * 100)
actual_celltypes <- unique(celltype_counts$sub_pred_EML)
color_mapping <- ref_panel[names(ref_panel) %in% actual_celltypes]
ggplot(celltype_counts, aes(x = "", y = count, fill = sub_pred_EML)) +
  geom_bar(stat = "identity", width = 1, color = "white") +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = color_mapping) +
  labs(
    title = "Cell Type Proportions",
    fill = "Cell Type"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )
developmental_order <- c(
  "Early_Hypoblast",
  "Late_Hypoblast",
  "Late_EPI",
  "Early_EPI",
  "PriS",
  "Mesoderm",
  "AdvMes",
  "Amnion",
  "ExE_Mes",
  "DE",
  "YSE",
  "HEP",
  "nb_failed",
  "Ambiguous",
  "low_cor"
)
ref_panel <- c(
  "Late_EPI"="#53B539",
  "Early_EPI"="#2E8B57",
  "Early_Hypoblast" ="#8195F8" ,
  "Late_Hypoblast" ="#B899D1" ,
  "ExE_Mes"="#3983E6",
  "HEP"="#6BE8EB",
  "DE"="#ABC8DB",
  "AdvMes"="#DC8941",
  "PriS"="#E76CD9",
  "YSE"="#DA6A8E",
  "Mesoderm"="#56BFA7",
  "Amnion"="#D5382A",
  "nb_failed"="#C0C0C0",
  "Ambiguous"="#EBEBEB",
  "low_cor"="#FFFFFF"
)
celltype_counts <- df %>%
  count(sub_pred_EML, name = "count") %>%
  mutate(percentage = count / sum(count) * 100) %>%
  arrange(desc(percentage)) 
celltype_counts$sub_pred_EML <- factor(
  celltype_counts$sub_pred_EML, 
  levels = celltype_counts$sub_pred_EML
)
celltype_counts <- celltype_counts %>%
  mutate(
    ypos = cumsum(count) - 0.5 * count
  )
ggplot(celltype_counts, aes(x = "", y = count, fill = sub_pred_EML)) +
  geom_bar(stat = "identity", width = 1, color = "white",alpha=0.9) +
  coord_polar("y", start = 0) +
  scale_fill_manual(values = ref_panel) +
  labs(
    title = "Cell Type Proportions",
    fill = "Cell Type"
  ) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  )

### fig5i. flow cytometric analysis of percentages 
## data grouping
# 单核/巨噬
# CD206-CD163+ 单核细胞
# CD206+CD163- 非经典单核
# CD206+CD163+ M2巨噬细胞
# CD42b-CD41+ 早期巨核
# CD42b+CD41- 巨核亚群
# CD42b+CD41+ 成熟巨核
df_summary_marco1 <- df_summary[c(1,2,3,4,5,6,7,8,9),]
df_summary_marco2 <- df_summary[c(16,17,18,19,20,21,22,23,24),]
p1
# 红系
# CD235a+CD71- 成熟红细胞
# CD235a+CD71+ 红系前体细胞
df_summary_ery <- df_summary[c(10,11,12,13,14,15,16,17,18),]
p2
# 造血干
# CD43-CD34+ 原始造血干（e.g.长期）
# CD43+CD34+ 造血祖细胞（e.g.多能祖细胞）
# CD45-CD34+ 内皮祖细胞/原始造血干细胞
# CD45+CD34+ 常见造血祖细胞
df_summary_stem <- df_summary[c(28,29,30,34,35,36),]
p3

merged_df <- df %>%
  mutate(merged_celltype = case_when(
    grepl("CD43[+-]CD34\\+", Label_celltypes) ~ "CD43_CD34+",
    grepl("CD45[+-]CD34\\+", Label_celltypes) ~ "CD45_CD34+",
    TRUE ~ Label_celltypes
  )) %>%
  group_by(merged_celltype, Stages) %>%
  summarise(mean_Value = mean(Value),
            sd_Value = sd(Value), 
            se_Value = sd(Value) / sqrt(n()), 
            n = n(),                  
            .groups = 'drop')
merged_df <- merged_df[c(25,26,27),]
merged_df$Stages <- factor(merged_df$Stages, levels = c("D6", "D10", "D14"))
merged_df$merged_celltype <- "CD34+"

library(readxl)
df <- read_excel("FACS_quantification.xlsx")
library(ggplot2)
library(dplyr)
df_summary <- df %>%
  group_by(Label_celltypes, Stages) %>%
  summarise(
    mean_Value = mean(Value),
    sd_Value = sd(Value),   
    se_Value = sd(Value) / sqrt(n()),  
    n = n(),                 
    .groups = 'drop'
  )

df_summary$Stages <- factor(df_summary$Stages, levels = c("D6", "D10", "D14"))
p3 <- ggplot(df_summary_stem, aes(x = Stages, y = mean_Value, 
                                  color = Label_celltypes, 
                                  group = Label_celltypes)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_Value - se_Value, 
                    ymax = mean_Value + se_Value),
                width = 0.1,  
                size = 0.6) + 
  labs(
    subtitle = "(±SEM)",
    x = "Stage",
    y = "Mean Value",
    color = "Celltype"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray30"),
    legend.position = "right",
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    panel.grid.major = element_line(color = "white", size = 0.3),
    panel.grid.minor = element_blank()
  )
print(p3)


### fig6. PCA of endodermal lineage ----
# data subset and merge
# day2 no endo
# day4 data
obj_day4_DE <- subset(obj_day4,idents = "Definitive endoderm")
obj_day4_DE <- CreateSeuratObject(counts = GetAssayData(obj_day4_DE),
                                  meta.data = data.frame(row.names = rownames(obj_day4_DE@meta.data),
                                                         orig.ident = "Day4",
                                                         nCount_RNA = obj_day4_DE@meta.data$nCount_RNA,
                                                         nFeature_RNA = obj_day4_DE@meta.data$nFeature_RNA,
                                                         celltype = obj_day4_DE@meta.data$cell_types_order,
                                                         sample = "Day4",
                                                         celltype_merge = paste("Day4",obj_day4_DE@meta.data$cell_types_order,sep = "_")
                                  ) )
rm(obj_day4)
# day6 data
DimPlot(obj_day6,group.by = "seurat_clusters",label = T)
Idents(obj_day6) <- obj_day6$seurat_clusters
obj_day6_Endo <- subset(obj_day6,idents = c("4","7","11","13","18","19"))
levels(obj_day6_Endo)
new.cluster.id <- c("YS endo progenitor",
                    "Def Endo",
                    "Endo gut","Endo gut","Endo gut",
                    "YS meso progenitor")
names(new.cluster.id) <- levels(obj_day6_Endo)
obj_day6_Endo <- RenameIdents(obj_day6_Endo, new.cluster.id)
DimPlot(obj_day6_Endo,label = T)
obj_day6_Endo$cell_type <- Idents(obj_day6_Endo)
obj_day6_Endo <- CreateSeuratObject(counts = GetAssayData(obj_day6_Endo),
                                    meta.data = data.frame(row.names = rownames(obj_day6_Endo@meta.data),
                                                           orig.ident = "Day6",
                                                           nCount_RNA = obj_day6_Endo@meta.data$nCount_RNA,
                                                           nFeature_RNA = obj_day6_Endo@meta.data$nFeature_RNA,
                                                           celltype = obj_day6_Endo@meta.data$cell_type,
                                                           sample = "Day6",
                                                           celltype_merge = paste("Day6",obj_day6_Endo@meta.data$cell_type,sep = "_")
                                    ) )
View(obj_day6_Endo@meta.data)
rm(obj_day6)
# day8 data
View(obj_day8@meta.data)
table(obj_day8$cell_types_ordered) 
DimPlot(obj_day8,label = T)
obj_day8_Endo <- subset(obj_day8,idents = c("Yolk-sac mesoderm","Yolk-sac endoderm","Endoderm-gut/DE"))
obj_day8_Endo <- CreateSeuratObject(counts = GetAssayData(obj_day8_Endo),
                                    meta.data = data.frame(row.names = rownames(obj_day8_Endo@meta.data),
                                                           orig.ident = "Day8",
                                                           nCount_RNA = obj_day8_Endo@meta.data$nCount_RNA,
                                                           nFeature_RNA = obj_day8_Endo@meta.data$nFeature_RNA,
                                                           celltype = obj_day8_Endo@meta.data$cell_types_ordered,
                                                           sample = "Day8",
                                                           celltype_merge = paste("Day8",obj_day8_Endo@meta.data$cell_types_ordered,sep = "_")
                                    ) )
View(obj_day8_Endo@meta.data)
rm(obj_day8)
# day10 data
table(obj_day10$cell_types_ordered)
obj_day10_Endo <- subset(obj_day10,idents = c("Yolk-sac mesoderm","Yolk-sac endoderm","Endoderm-gut","Definitive Endoderm"))
obj_day10_Endo <- CreateSeuratObject(counts = GetAssayData(obj_day10_Endo),
                                     meta.data = data.frame(row.names = rownames(obj_day10_Endo@meta.data),
                                                            orig.ident = "Day10",
                                                            nCount_RNA = obj_day10_Endo@meta.data$nCount_RNA,
                                                            nFeature_RNA = obj_day10_Endo@meta.data$nFeature_RNA,
                                                            celltype = obj_day10_Endo@meta.data$cell_types_ordered,
                                                            sample = "Day10",
                                                            celltype_merge = paste("Day10",obj_day10_Endo@meta.data$cell_types_ordered,sep = "_")))
View(obj_day10_Endo@meta.data)
rm(obj_day10)
# day15 data
table(obj_day15$cell_types_ordered)
Idents(obj_day15) <- obj_day15$cell_types_ordered
obj_day15_Endo <- subset(obj_day15,idents = c("Yolk-sac endoderm","Endoderm-gut","Definitive Endoderm","Liver-like","Stomach-like","Pancreas-like"))
obj_day15_Endo <- CreateSeuratObject(counts = GetAssayData(obj_day15_Endo),
                                     meta.data = data.frame(row.names = rownames(obj_day15_Endo@meta.data),
                                                            orig.ident = "Day15",
                                                            nCount_RNA = obj_day15_Endo@meta.data$nCount_RNA,
                                                            nFeature_RNA = obj_day15_Endo@meta.data$nFeature_RNA,
                                                            celltype = obj_day15_Endo@meta.data$cell_types_ordered,
                                                            sample = "Day15",
                                                            celltype_merge = paste("Day15",obj_day15_Endo@meta.data$cell_types_ordered,sep = "_")))
View(obj_day15_Endo@meta.data)   
rm(obj_day15)
#YLQ CS7 data
obj_CS7 <- readRDS('/Users/chenchenyi/Downloads/0000_存档/0002_human\ lineage/CS7_human_embryo.rds')
View(obj_CS7@meta.data)
table(obj_CS7$clusters)
Idents(obj_CS7) <- obj_CS7$clusters
obj_CS7_Endo <- subset(obj_CS7,idents = c("DE/VE","Em/EXE.Meso","YS.DE/VE","YS.EXE.Em/EXE.Meso"))
obj_CS7_Endo <- CreateSeuratObject(counts = GetAssayData(obj_CS7_Endo),
                                   meta.data = data.frame(row.names = rownames(obj_CS7_Endo@meta.data),
                                                          orig.ident = "CS7",
                                                          nCount_RNA = obj_CS7_Endo@meta.data$nCount_RNA,
                                                          nFeature_RNA = obj_CS7_Endo@meta.data$nFeature_RNA,
                                                          celltype = obj_CS7_Endo@meta.data$clusters,
                                                          sample = "CS7",
                                                          celltype_merge = paste("CS7",obj_CS7_Endo@meta.data$clusters,sep = "_")))
View(obj_CS7_Endo@meta.data)
rm(obj_CS7)
#YLQ CS8 data
table(obj_CS8$celltype)
Idents(obj_CS8) <- obj_CS8$celltype
obj_CS8_Endo <- subset(obj_CS8,idents = c("Endo","Visceral Endo","YS Endo","YS EXM meso-A","YS EXM meso-B"))
obj_CS8_Endo <- CreateSeuratObject(counts = GetAssayData(obj_CS8_Endo),
                                   meta.data = data.frame(row.names = rownames(obj_CS8_Endo@meta.data),
                                                          orig.ident = "CS8",
                                                          nCount_RNA = obj_CS8_Endo@meta.data$nCount_RNA,
                                                          nFeature_RNA = obj_CS8_Endo@meta.data$nFeature_RNA,
                                                          celltype = obj_CS8_Endo@meta.data$celltype,
                                                          sample = "CS8",
                                                          celltype_merge = paste("CS8",obj_CS8_Endo@meta.data$celltype,sep = "_")))
View(obj_CS8_Endo@meta.data)
rm(obj_CS8)
# list
obj_day4_Endo <- obj_day4_DE
obj_Endo_list <- list(
  obj_day4_Endo,
  obj_day6_Endo,
  obj_day8_Endo,
  obj_day10_Endo,
  obj_day15_Endo,
  obj_CS7_Endo,
  obj_CS8_Endo
)
names(obj_Endo_list) <- c("Day4",
                          "Day6",
                          "Day8",
                          "Day10",
                          "Day15",
                          "CS7",
                          "CS8")
rm(obj_day4_Endo,
   obj_day6_Endo,
   obj_day8_Endo,
   obj_day10_Endo,
   obj_day15_Endo,
   obj_CS7_Endo,
   obj_CS8_Endo)
saveRDS(obj_Endo_list,"1202_obj_Endo_list_for_PCA.rds")
# merge
library(dplyr)
gene.merge.min <- lapply(obj_Endo_list, function(x){
  rownames(x)
}) %>% Reduce(f = intersect)
obj_Endo_list_sub <- lapply(obj_Endo_list,function(x){
  subset(x,features = gene.merge.min)
})
bulk_ag <- lapply(obj_Endo_list_sub,function(x){
  res <- AggregateExpression(x,group.by = "celltype_merge")$RNA
  if (ncol(res) == 1) {
    celltype <- unique(x$celltype_merge)
    colnames(res) <- celltype
  }
  res
}) %>% purrr::reduce(cbind) %>% as.data.frame
colnames(bulk_ag) <- gsub("-", "_", colnames(bulk_ag))   
saveRDS(bulk_ag,"bulk_ag_Endo_YS_CS78_min.rds")

## pseudobulk and PCA 
library(ggplot2)
library(ggrepel)
library(dplyr)
library(sva) 
set.seed(123)
bulk <- bulk_ag[,c(4,7,11,18,19,21,23,24,25)]
normalized_data <- log2(bulk + 1)
group_info <- c(rep("Group1", 12), rep("Group2", 5))
cell_types <- colnames(bulk)
batch_info
# combat
library(sva)
mod <- model.matrix(~ group_info)
corrected_data <- ComBat(dat = as.matrix(normalized_data), 
                         batch = batch_info, 
                         mod = NULL,
                         par.prior = TRUE)
pca_result <- prcomp(t(corrected_data), scale. = TRUE, center = TRUE)
pca_df <- as.data.frame(pca_result$x)
pca_df$Group <- group_info
pca_df$CellType <- cell_types
pca_df$Batch <- batch_info

variance_explained <- round(100 * pca_result$sdev^2 / sum(pca_result$sdev^2), 2)

library(ggthemes)
p1 <- ggplot(pca_df, aes(x = PC1, y = PC3, color = Group, 
                         # shape = Batch, 
                         label = CellType)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  scale_color_manual(
    values = c("Group1" = "#E41A1C", "Group2" = "#377EB8"),
    labels = c("YS model", "reference")
  ) +
  # scale_shape_manual(values = c(16, 17)) +
  labs(
    # title = "PCA analysis",
    x = paste0("PC1 (", variance_explained[1], "%)"),
    y = paste0("PC3 (", variance_explained[3], "%)"),
    color = "Group"#,
    #shape = "batch"
  ) +
  theme_few() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "bottom"
  )
print(p1)
# save
write.csv(pca_df, "pca_results_batch_corrected.csv", row.names = TRUE)
ggsave("PCA_batch_corrected.pdf", p1, width = 10, height = 8)

### figS4. sankey diagram of cluster distribution ----
library(ggalluvial)
df <- seu.merge@meta.data
alluvial_data <- df %>%
  count(stage, sub_pred_EML) %>%
  group_by(stage) %>%
  mutate(percentage = n / sum(n) * 100) %>%
  ungroup()
alluvial_data <- celltype_percentage
alluvial_data$stage <- gsub("D","Day",alluvial_data$stage)
alluvial_data$stage <- factor(alluvial_data$stage, levels = c("Day2","Day4","Day6","Day8","Day10","Day15"))
alluvial_data$sub_pred_EML <- factor(alluvial_data$sub_pred_EML,
                                     levels = c(
                                       "Amnion",
                                       "YSE",
                                    
                                       "Late_EPI",
                                       "PriS",
                                       "Mesoderm",
                                       "Early_EPI",
                                       
                                       "Early_Hypoblast",
                                       "Late_Hypoblast",
                                       "DE",
                                       
                                       "AdvMes",
                                       "ExE_Mes",
                                       "HEP",
                                       
                                       "nb_failed",
                                       "Ambiguous"
                                     ))
cell_type_colors <- c(
  "Amnion"="#D5382A",
  "YSE"="#DA6A8E",
  
  "Late_EPI"="#53B539",
  "PriS"="#E76CD9",
  "Mesoderm"="#56BFA7",
  
  "Early_EPI"="#2E8B57",
  
  "Early_Hypoblast" ="#8195F8" ,
  "Late_Hypoblast" ="#B899D1" ,
  "DE"="#ABC8DB",
  
  "AdvMes"="#DC8941",
  "ExE_Mes"="#3983E6",
  "HEP"="#6BE8EB",
  
  "nb_failed"="#C0C0C0",
  "Ambiguous"="#EBEBEB"
)
names(cell_type_colors) <- levels(alluvial_data$sub_pred_EML)
ggplot(alluvial_data,
       aes(x = stage, 
           y = percentage, 
           stratum = sub_pred_EML, 
           alluvium = sub_pred_EML,
           fill = sub_pred_EML,
           label = ifelse(percentage > 1.6, paste0(round(percentage, 1), "%"), "")
       )
) +
  geom_flow(alpha = 0.8) +
  geom_stratum(alpha = 0.8) +
  geom_text(stat = "stratum", size = 2) +
  scale_fill_manual(
    values = cell_type_colors,
    name = "Cell Type",
    labels = levels(alluvial_data$sub_pred_EML)  # 可以自定义标签
  ) +
  labs(
    x = "Stage",
    y = "Percentage within Stage (%)"#,
  ) +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),  
    axis.title = element_text(size = 11, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 1, size = 10),
    axis.ticks.x = element_blank()
  )

### fig.s5h percentage of YS endoderm ----
library(dplyr)
library(tidyr)
library(ggplot2)
plot_data <- obj@meta.data %>%
  group_by(stage) %>%
  summarise(
    total_cells = n(),
    ys_cells = sum(cell_type_ordered == "YS endoderm"),
    ys_prop = ys_cells / total_cells * 100,
    .groups = "drop"
  ) %>%
  filter(!stage %in% c("D2"))
plot_data$stage <- factor(
  plot_data$stage,
  levels = c("D4", "D6", "D8", "D10", "D15")
)
plot_data$x <- seq_len(nrow(plot_data))
ymax <- 5
plot_long <- plot_data %>%
  mutate(
    `YS endoderm` = ys_prop,
    `Other cells` = ymax - ys_prop
  ) %>%
  select(x, stage, `YS endoderm`, `Other cells`) %>%
  pivot_longer(
    cols = c(`YS endoderm`, `Other cells`),
    names_to = "cell_group",
    values_to = "proportion"
  )
plot_long$cell_group <- factor(
  plot_long$cell_group,
  levels = c("Other cells", "YS endoderm")
)
xmin_box <- 1
xmax_box <- nrow(plot_data)
xpad <- 0.20
ypad <- 0.15
xmin_axis <- xmin_box - xpad
xmax_axis <- xmax_box + xpad
line_width <- 0.55
ggplot(
  plot_long,
  aes(
    x = x,
    y = proportion,
    fill = cell_group,
    group = cell_group
  )
) +
  geom_area(
    position = "stack",
    alpha = 0.9,
    linewidth = 0
  ) +
  annotate(
    "rect",
    xmin = xmin_box,
    xmax = xmax_box,
    ymin = 0,
    ymax = ymax,
    fill = NA,
    colour = "grey20",
    linewidth = line_width
  ) +
  geom_text(
    data = plot_data %>%
      mutate(label_x = case_when( stage == "D4"  ~ x + 0.26,stage == "D15" ~ x - 0.26, TRUE ~ x ) ),
    aes(x = label_x, y = ys_prop,label = sprintf("%.1f%%", ys_prop)),
    inherit.aes = FALSE,nudge_y = 0.18,hjust = 0.5,size = 3.7) +
  scale_fill_manual(
    values = c( "YS endoderm" = "#B65655","Other cells" = "#EEEEEE"),
    breaks = c("YS endoderm", "Other cells")) +
  scale_x_continuous(limits = c(xmin_axis, xmax_axis),breaks = plot_data$x,labels = as.character(plot_data$stage),expand = c(0, 0)) +
  scale_y_continuous(limits = c(-ypad, ymax + ypad),breaks = 0:ymax,labels = 0:ymax,expand = c(0, 0)) +
  labs(x = NULL,y = "YS endoderm proportion (%)",fill = NULL) +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 13) +
  theme(
    axis.line = element_blank(),
    panel.border = element_rect(colour = "grey20",fill = NA,linewidth = line_width),
    axis.ticks = element_line(colour = "grey20",linewidth = line_width),
    axis.ticks.length = unit(0.16, "cm"),
    axis.text = element_text(colour = "grey15",size = 11),
    axis.title.y = element_text(size = 12,margin = margin(r = 8)),
    panel.grid = element_blank(),
    legend.position = "top",
    legend.text = element_text(size = 11),
    plot.margin = margin(t = 8,r = 10,b = 8,l = 10))

### fig3a/d,fig5f/g/h In vivo YS data reference ----

library(Seurat)
Idents(ys)
DimPlot(ys)
View(ys@meta.data)
ys@meta.data <- meta
DimPlot(ys,group.by = 'annotation1014',label = T)

marker.fig2a <- c(
  #endo
  'SERPINA1','APOA2','АРОСЗ','SPINK1','AHSG',
  #hepatocyte，
  'AFP','APOA1','ALB','TTR','FGB',
  #fibroblast YS
  'LUM','FRZB','COL3A1','SPARC','COL6A2',
  #fibroblast liver
  'PTN','COL1A1','VIM','TPM1','CXCL12',
  #smooth muscle
  'CSRP2','CTHRC1','TPM2','CALD1','ACTA2',
  #mesothelium
  'KRT19','S100A10','KRT8','TMEM98','PDPN','UPK3B')

#YS endoliver endo
Idents(ys) <- ys$annotation1014
ys.endo <- subset(ys,
                  idents = c('Gut','Definitive endoderm','Pancreas-like','YS endoderm','Stomach-like')
)
marker.fig2d <- c(
  #lipid-metabolic
  'SERPINA1','APOA2','APOC3','SPINK1','AHSG','AFP','APOA1','ALB','TTR',
  #hemato-common pathway
  'F2','F5','F10',
  #hemato-extrinstic pathway
  'F3','F7',
  #hemato-intrinsic pathway
  'F12','F11','F8',
  #hemato-fibrinogen & crosslinking
  'FGA','FGB','FGG','F13A1',
  #hemato-anticoagulation
  'SERPINC1','PROC','PROS1',
  #hemato-fibrinolysis
  'PLAT','PLAU','SERPINE1','SERPINB2',
  #growth factor
  'EPO' ,'EGFR' ,'THPO'
)

marker.fig2e <- c(
  #lipid metabolic
  'APOE',
  'APOM',
  'APOA1',
  'AFP',
  'FABP1',
  'FASN',
  'CYP51A1',
  'APOB',
  'GPAM',
  'SREBF1',
  'LSS',
  'ALB',
  #retinoid metabolism
  'RBP4',
  'RBP2',
  'RBP3',
  'RARRES2',
  #fibrinogen and cross linking
  'FGB',
  'FGG',
  'FGA',
  #apoptosis regulation
  'NFKBIA',
  'JUN',
  'ATF3'
)
marker.fig3a <- c(
  #Canonical
  'CD34',
  'MLLT3',
  'SPINK2',
  'HOPX',
  'HLF',
  'RAB27B',
  'MYB',
  #early
  'DDIT4',
  'SLC2A3',
  'RGS16',
  'LIN28A',
  #definitive
  'KIT',
  'ITGA4',
  'CD74',
  'PROCR',
  'EMCN',
  'GBP4',
  'ACE',
  #patterning
  'HOXA7',
  'HOXA9',
  'HOXA10',
  'HOXB7',
  'HOXB9')
marker.fig4c <- c(
  'PLVAP',
  'ESAM',
  'PECAM1',
  'CDH5',
  'FLT1',
  'KDR',
  'KCNK17'
)
marker.fig4d <- c(
  #bottom
  'KITLG',
  'DLL1',
  'DLL1',
  'IGF2',
  'JAG1',
  'FBN1',
  'WNT5A',
  'FN1',
  'FN1',
  'DLK1',
  'EPO',
  'THPO',
  'VTN',
  #top
  'KIT',
  'NOTCH1',
  'NOTCH2',
  'IGF1R',
  'NOTCH2',
  'ITGA5',
  'ITGB1',
  'FZD3',
  'ITGAV',
  'ITGB1',
  'ITGA4',
  'ITGB1',
  'NOTCH4',
  'EPOR',
  'MPL',
  'ITGAV',
  'ITGB1'
)
DotPlot(
  ys.endo,
  features = marker.fig2d
)RotatedAxis()


DotPlot(obj,features = c('CD14','CD64','TREM2'))

# marker gene dotplot
ref <- readRDS("ys_seurat_v4.rds")
obj <- readRDS("BMP4_YS.combined_all.rds")
obj@meta.data <- meta
Idents(obj) <- obj$cell_type_ordered
DimPlot(obj,label = T)
# obj.subset <- subset(obj,idents = c('YS endoderm','YS medoderm','Gut','Stomach-like','Pancreas-like','Definitive endoderm'))
obj.subset <- subset(obj,idents = c('YS endoderm','YS medoderm'))
# table(obj.subset$bulk.label)
# Idents(obj.subset) <- obj.subset$bulk.label

# View(ref@meta.data)
Idents(ref) <- ref$LVL2 
# table(ref.subset$bulk.label)
# Idents(ref.subset) <- ref.subset$bulk.label

bulk.obj.stage.type.av <- as.data.frame((AverageExpression(obj.subset,group.by = c('cell_type_ordered')))$RNA)
bulk.ref.stage.type.av <- as.data.frame((AverageExpression(ref.subset,group.by = c('LVL2')))$RNA)

#correlation by all genes
df <- merge(bulk.obj.stage.type.av,bulk.ref.stage.type.av, by = 0)
rownames(df) <- df$Row.names
df$Row.names <- NULL
df <- df[, c(
  "YS endoderm",
  "YS medoderm",
  
  "ENDODERM",
  "SMOOTH-MUSCLE",
  "FIBROBLAST",
  "MESOTHELIUM"
  
), drop = FALSE]

df_normalized <- log2(1+df)
correlation_matrix <- sapply(3:ncol(df_normalized), function(i) {
  cor(df_normalized[, 1:2], df_normalized[, i], method = "spearman")
})
# correlation_matrix <- cor(df_normalized, method = "spearman")
correlation_matrix <- t(correlation_matrix)
rownames(correlation_matrix) <- colnames(df_normalized)[-c(1:2)]
colnames(correlation_matrix) <- colnames(df_normalized)[1:2]

mat <- t(correlation_matrix)
vals <- as.vector(mat)

minv <- min(vals, na.rm = TRUE)
maxv <- max(vals, na.rm = TRUE)

pheatmap::pheatmap(
  mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("lightblue", "white", "pink"))(100),
  breaks = seq(minv, maxv, length.out = 101),
  display_numbers = F
)


corrplot::corrplot(mat,
                   method = 'shade',
                   col = colorRampPalette(c("lightblue", "white", "pink"))(100),
                   is.corr=T,
                   addCoef.col=F
)

## marker gene dotplot
meta <- readRDS("2026_ys_metadata.rds")
obj@meta.data <- meta
Idents(obj) <- obj$cell_type_ordered
obj.subset <- subset(obj,idents = c('YS endoderm'))
obj.subset$type.stage <- paste(obj.subset$stage,obj.subset$cell_type_ordered)

ref.subset <- subset(ref,idents = c('ENDODERM'))
ref.subset$type.stage <- paste(ref.subset$stage,ref.subset$LVL2)

Idents(obj.subset) <- 'ys'
obj.subset$group <- 'Yolk sac model'
Idents(ref.subset) <- 'ref'
ref.subset$group <- 'Yolk sac reference'
merge <- merge(obj.subset,ref.subset,add.cell.ids = c('ys','ref'))
Idents(merge) <- merge$group

#dotplot
library(ggplot2)
library(grid)
library(dplyr)
library(tibble)
marker.fig2d.list <- list(
  "lipid-metabolic" = c(
    "SERPINA1","APOA2","APOC3","SPINK1","AHSG","AFP","APOA1","ALB","TTR" ),
  "hemato-common pathway" = c(
    "F2","F5","F10"),
  "hemato-extrinsic pathway" = c(
    "F3","F7"),
  "hemato-intrinsic pathway" = c(
    "F12","F11","F8"),
  "hemato-fibrinogen & crosslinking" = c(
    "FGA","FGB","FGG","F13A1" ),
  "hemato-anticoagulation" = c(
    "SERPINC1","PROC","PROS1" ),
  "hemato-fibrinolysis" = c(
    "PLAT","PLAU","SERPINE1","SERPINB2"),
  "growth factor" = c(
    "EPO","EGFR","THPO"))
marker_df <- stack(marker.fig2d.list)
colnames(marker_df) <- c("gene", "cluster")
merge$stage <- factor((merge$stage),levels=c("D6","D8","D10","D15","CS10","CS11","CS14","CS15","CS17","CS18","CS22","CS23"))
library(tibble)
marker.fig2d <- tibble(
  cluster = rev(c(
    rep("Lipid metabolic", 9),
    rep("Haemostasis-common pathway", 3),
    rep("Haemostasis-extrinsic pathway", 2),
    rep("Haemostasis-intrinsic pathway", 3),
    rep("Haemostasis-fibrinogen & crosslinking", 4),
    rep("Haemostasis-anticoagulation", 3),
    rep("Haemostasis-fibrinolysis", 4),
    rep("Growth factors", 3)
  )),
  gene = rev(c(
    "SERPINA1","APOA2","APOC3","SPINK1","AHSG","AFP","APOA1","ALB","TTR",
    "F2","F5","F10",
    "F3","F7",
    "F12","F11","F8",
    "FGA","FGB","FGG","F13A1",
    "SERPINC1","PROC","PROS1",
    "PLAT","PLAU","SERPINE1","SERPINB2",
    "EPO","EGFR","THPO"
  ))
)

jjDotPlot(
  object = obj,
  markerGene = marker.fig2d,
  id = "stage",
  scale = T,
  base_size = 15,
  legend.position = "left",
  ytree = F,
  bar.width = 3,
  plot.margin = c(1, 8, 1, 1),
  anno = TRUE,anno_y_position =12, aesGroName = "cluster",segWidth = 0.8,lwd = 1.5,textRot = 0,textSize = 12,hjust =0)

marker.fig2e <- data.frame(
  gene = c(
    # lipid metabolic
    'APOE',
    'APOM',
    'APOA1',
    'AFP',
    'FABP1',
    'FASN',
    'CYP51A1',
    'APOB',
    'GPAM',
    'SREBF1',
    'LSS',
    'ALB',
    # retinoid metabolism
    'RBP4',
    'RBP2',
    'RBP3',
    'RARRES2',
    # fibrinogen and cross linking
    'FGB',
    'FGG',
    'FGA',
    # apoptosis regulation
    'NFKBIA',
    'JUN',
    'ATF3'
  ),
  cluster = c(
    rep('Lipid metabolic', 12),
    rep('Retinoid metabolism', 4),
    rep('Fibrinogen and cross linking', 3),
    rep('Apoptosis regulation', 3)
  ),
  stringsAsFactors = FALSE
)
Idents(merge) <- merge$stage
merge.sub.n.cs23 <- subset(merge,idents = c('CS23'),invert = T)
jjDotPlot(
  object = merge.sub.n.cs23,markerGene = marker.fig2e,id = "stage",scale = T,
  base_size = 15,legend.position = "left",ytree = F,bar.width = 3,plot.margin = c(5, 2, 1, 1),
  anno = TRUE,anno_y_position = 11.8, aesGroName = "cluster",segWidth = 0.2,lwd = 1.5,textRot = 45,textSize = 12,hjust =0)

marker.fig3a <- tibble(
  cluster = rev(c(
    rep("Patterning", 5),
    rep("Definitive", 7),
    rep("Early", 4),
    rep("Canonical", 7)
  )),
  gene = rev(c(
    'HOXA7','HOXA9','HOXA10','HOXB7','HOXB9',
    'KIT','ITGA4','CD74','PROCR','EMCN','GBP4','ACE',
    'DDIT4','SLC2A3','RGS16','LIN28A',
    'CD34','MLLT3','SPINK2','HOPX','HLF','RAB27B','MYB'
  ))
)

jjDotPlot(
  object = obj.sub,markerGene = marker.fig3a,id = "stage",
  scale = T,base_size = 15,legend.position = "left",ytree = F,bar.width = 3,
  plot.margin = c(1, 4, 1, 1),anno = TRUE,anno_y_position = 4.88,  aesGroName = "cluster",
  segWidth = 0.8,lwd = 1.5,textRot = 0,textSize = 12,hjust =0,bar.width = 6
)+coord_flip(clip = "off")


### figS5 correlation heatmap ----

bulk.hed #merged YS blastoids peseudo-bulk data

###CS7
# CS7.cts <- AggregateExpression(CS7, group.by = "sub_cluster",assays = "RNA")
# CS7.cts <- as.data.frame(CS7.cts$RNA)
# CS7.cts[] <- lapply(CS7.cts, as.integer)
# rm(CS7)
# colnames(CS7.cts) <- paste("CS7",colnames(CS7.cts),sep = "_")
# saveRDS(CS7.cts,"CS7.cts.rds")
bulk.CS7

###CS8
# CS8.cts <- AggregateExpression(CS8, group.by = "celltype",assays = "RNA")
# CS8.cts <- as.data.frame(CS8.cts$RNA)
# colnames(CS8.cts) <- paste("CS8",colnames(CS8.cts),sep = "_")
# rm(CS8)
# saveRDS(CS8.cts,"CS8.cts.rds")
bulk.CS8

###CS10
# CS10.cts <- AggregateExpression(CS10, group.by = "cell_type",assays = "RNA")
# CS10.cts <- as.data.frame(CS10.cts$RNA)
# colnames(CS10.cts) <- paste("CS10",colnames(CS10.cts),sep = "_")
# rm(CS10)
# saveRDS(CS10.cts,"CS10.cts.rds")
bulk.CS10

###CS12 SWY
# CS12 <- subset(human_seurat_CS12_annotated,subset = stage=="CS12")
# CS13_14 <- subset(SWY,subset = stage=="CS13-14")
# CS15_16 <- subset(SWY,subset = stage=="CS15-16")
# rm(human_seurat_CS12_annotated)
# rm(SWY)
# CS12.cts <- AggregateExpression(CS12, group.by = "developmental.system",assays = "RNA")
# CS13_14.cts <- AggregateExpression(CS13_14, group.by = "developmental.system",assays = "RNA")
# CS15_16.cts <- AggregateExpression(CS15_16, group.by = "developmental.system",assays = "RNA")
# CS12.cts <- as.data.frame(CS12.cts$RNA)
# CS13_14.cts <- as.data.frame(CS13_14.cts$RNA)
# CS15_16.cts <- as.data.frame(CS15_16.cts$RNA)
# colnames(CS12.cts) <- paste("CS12",colnames(CS12.cts),sep = "_")
# colnames(CS13_14.cts) <- paste("CS13-14",colnames(CS13_14.cts),sep = "_")
# colnames(CS15_16.cts) <- paste("CS15-16",colnames(CS15_16.cts),sep = "_")
# bulk.SWY <- cbind(bulk.CS12,CS13_14.cts,CS15_16.cts)
bulk.SWY

###pre-implantation
# preimplannation <- AggregateExpression(preimplannation, group.by = "Group",assays = "RNA")
# preimplannation.cts <- as.data.frame(preimplannation$RNA)
# colnames(preimplannation.cts) <- paste("preimplantation",colnames(preimplannation.cts),sep = "_")
# rm(preimplannation)
# saveRDS(preimplannation.cts,"preimplannation.cts.rds")
bulk.ltq

###CS7-all
gene.trans <- genename.list[["hed_CS7_ltq"]]
bulk.CS7.trans <- bulk.CS7[gene.trans$CS7,]
rownames(bulk.CS7.trans) <- gene.trans$hed
bulk.hed.trans <- bulk.hed[gene.trans$hed,]
bulk.hed.CS7 <- cbind(bulk.CS7.trans,bulk.hed.trans)

###CS8-all
# colnames(bulk.CS8) <- paste("CS8",colnames(bulk.CS8),sep = "_")
gene.CS8 <- genename.list[[4]]
gene.hed <- rownames(bulk.hed)
gene.CS8.hed <- intersect(gene.CS8$CS8, gene.hed)
bulk.hed.CS8 <- cbind(bulk.CS8[gene.CS8.hed, ],bulk.hed[gene.CS8.hed,])

###CS10-all
gene.CS10.hed <- genename.list[["CS10.hed"]]
bulk.CS10.trans <- bulk.CS10[gene.CS10.hed$gene.CS10.ori,]
rownames(bulk.CS10.trans) <- gene.CS10.hed$gene.hed.ori
bulk.hed.CS10 <- merge(bulk.CS10.trans,bulk.hed,by = "row.names")
rownames(bulk.hed.CS10) <- bulk.hed.CS10$Row.names
bulk.hed.CS10$Row.names <- NULL

###CS12-all
gene.SWY.hed <- genename.list[["CS12.hed"]]
bulk.SWY <- bulk.SWY[gene.SWY.hed$CS12,]
rownames(bulk.SWY) <- gene.SWY.hed$gene.hed.ori
bulk.hed.SWY <- merge(bulk.SWY,bulk.hed,by = "row.names")
rownames(bulk.hed.SWY) <- bulk.hed.SWY$Row.names
bulk.hed.SWY$Row.names <- NULL
rm(bulk.SWY,bulk.SWY.hed,gene.SWY.hed,gene.SWY)

###ltq-all
gene.hed.ltq <- genename.list[["hed_CS7_ltq"]]
bulk.ltq.trans <- bulk.ltq[gene.hed.ltq$ltq,]
colnames(bulk.ltq.trans) <- paste("preimplantation",colnames(bulk.ltq.trans),sep = "_")
rownames(bulk.ltq.trans) <- gene.hed.ltq$hed
bulk.hed.trans <- bulk.hed[gene.hed.ltq$hed,]
bulk.hed.ltq <- cbind(bulk.ltq.trans,bulk.hed.trans)
rm(gene.hed.ltq,bulk.ltq,bulk.ltq.trans,bulk.hed.trans)

df #merged
# 计算每个基因的极差
# gene_variability <- apply(df, 1, function(x) max(x) - min(x))
gene_variability <- apply(df, 1, sd)
top_genes <- sort(gene_variability, decreasing = TRUE)[1:50]
top_genes_df <- df[names(top_genes), ]
df <- top_genes_df

#top30 DEG
top_genes_per_cluster_list <- list()
data <- marker.all
top_genes_per_cluster <- data %>%
  group_by(cluster) %>%
  top_n(30, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))
top_genes_per_cluster_1 <- unique(top_genes_per_cluster$gene)
top_genes_per_cluster_list[["top_hed"]] <- top_genes_per_cluster_1
rm(data,top_genes_per_cluster)

#CS7
data <- CS7_markers
top_genes_per_cluster <- data %>%
  group_by(cluster) %>%
  top_n(30, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))
top_genes_per_cluster_2 <- unique(top_genes_per_cluster$gene)
top_genes_per_cluster_list[["top_CS7"]] <- top_genes_per_cluster_2
rm(data,top_genes_per_cluster,top_genes_per_cluster_2,CS7_markers)

top_genes_per_cluster <- c(top_genes_per_cluster_1,top_genes_per_cluster_2)
top_genes_per_cluster <- unique(top_genes_per_cluster)
top_genes_per_cluster_list[["top_hed_CS7"]] <- top_genes_per_cluster
rm(top_genes_per_cluster_2)

#CS8
data <- CS8_markers
top_genes_per_cluster <- data %>%
  group_by(cluster) %>%
  top_n(30, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))
top_genes_per_cluster_2 <- unique(top_genes_per_cluster$gene)
top_genes_per_cluster_list[["top_CS8"]] <- top_genes_per_cluster_2
rm(data,top_genes_per_cluster,CS8_markers)

top_genes_per_cluster <- c(top_genes_per_cluster_1,top_genes_per_cluster_2)
top_genes_per_cluster <- unique(top_genes_per_cluster)
top_genes_per_cluster_list[["top_hed_CS8"]] <- top_genes_per_cluster
rm(top_genes_per_cluster_2,top_genes_per_cluster)

#ltq
data <- marker_ltq
top_genes_per_cluster <- data %>%
  group_by(cluster) %>%
  top_n(30, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))
top_genes_per_cluster_2 <- unique(top_genes_per_cluster$gene)
top_genes_per_cluster_list[["top_ltq"]] <- top_genes_per_cluster_2
rm(data,top_genes_per_cluster,marker_ltq)

top_genes_per_cluster <- c(top_genes_per_cluster_1,top_genes_per_cluster_2)
top_genes_per_cluster <- unique(top_genes_per_cluster)
top_genes_per_cluster_list[["top_hed_ltq"]] <- top_genes_per_cluster
rm(top_genes_per_cluster_2,top_genes_per_cluster)

#CS10
data <- CS10_markers
top_genes_per_cluster <- data %>%
  group_by(cluster) %>%
  top_n(30, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))
top_genes_per_cluster_2 <- unique(top_genes_per_cluster$gene)
top_genes_per_cluster_list[["top_CS10"]] <- top_genes_per_cluster_2
rm(data,top_genes_per_cluster,CS10_markers)

top_genes_per_cluster <- c(top_genes_per_cluster_1,top_genes_per_cluster_2)
top_genes_per_cluster <- unique(top_genes_per_cluster)
top_genes_per_cluster_list[["top_hed_CS10"]] <- top_genes_per_cluster
rm(top_genes_per_cluster_2,top_genes_per_cluster)

#SWY
data <- SWY_markers
top_genes_per_cluster <- data %>%
  group_by(cluster) %>%
  top_n(30, avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))
top_genes_per_cluster_2 <- unique(top_genes_per_cluster$gene)
top_genes_per_cluster_list[["top_SWY"]] <- top_genes_per_cluster_2
rm(data,top_genes_per_cluster,SWY_markers)

top_genes_per_cluster <- c(top_genes_per_cluster_1,top_genes_per_cluster_2)
top_genes_per_cluster <- unique(top_genes_per_cluster)
top_genes_per_cluster_list[["top_hed_SWY"]] <- top_genes_per_cluster
rm(top_genes_per_cluster_2,top_genes_per_cluster)

hed_order <- c("d2_early PGC-like","d4_early PGC-like","d6_early PGC-like","d8_early PGC-like","d10_early PGC-like","d15_early PGC-like",
               "d4_late PGC-like","d6_late PGC-like","d8_late PGC-like","d10_late PGC-like","d15_late PGC-like",
               "d2_Epiblast","d4_Epiblast","d6_Epiblast","d8_Epiblast","d10_Epiblast","d15_Epiblast",
               "d2_Transiting epiblast","d4_Transiting epiblast","d6_Transiting epiblast","d8_Transiting epiblast","d10_Transiting epiblast","d15_Transiting epiblast",
               "d2_Primitive Streak","d4_Primitive Streak","d6_Primitive Streak","d8_Primitive Streak","d10_Primitive Streak","d15_Primitive Streak",
               "d4_Nascent mesoderm","d6_Nascent mesoderm","d8_Nascent mesoderm","d10_Nascent mesoderm","d15_Nascent mesoderm",
               "d4_Advanced mesoderm","d6_Advanced mesoderm","d8_Advanced mesoderm","d10_Advanced mesoderm","d15_Advanced mesoderm",                                          
               "d2_Cardiomyocyte","d4_Cardiomyocyte","d6_Cardiomyocyte","d8_Cardiomyocyte","d10_Cardiomyocyte","d15_Cardiomyocyte",
               "d2_Mesenchyme","d4_Mesenchyme","d6_Mesenchyme","d8_Mesenchyme","d10_Mesenchyme","d15_Mesenchyme",
               "d4_Limb-like","d6_Limb-like","d8_Limb-like","d10_Limb-like","d15_Limb-like",
               "d6_Somite LPM","d8_Somite LPM","d10_Somite LPM","d15_Somite LPM",                   
               "d4_Erythrocyte-like","d6_Erythrocyte-like","d8_Erythrocyte-like","d10_Erythrocyte-like","d15_Erythrocyte-like",
               "d4_Blood vessel","d6_Blood vessel","d8_Blood vessel","d10_Blood vessel","d15_Blood vessel",
               "d6_Blood progenitor","d8_Blood progenitor","d10_Blood progenitor","d15_Blood progenitor",                                                     
               "d4_Definitive endoderm","d6_Definitive endoderm","d8_Definitive endoderm","d10_Definitive endoderm","d15_Definitive endoderm",
               "d4_Gut","d6_Gut","d8_Gut","d10_Gut","d15_Gut",
               "d6_Pancreas-like","d8_Pancreas-like","d10_Pancreas-like","d15_Pancreas-like",
               "d4_Stomach-like","d6_Stomach-like","d8_Stomach-like","d10_Stomach-like","d15_Stomach-like",  
               "d6_YS endoderm","d8_YS endoderm","d10_YS endoderm","d15_YS endoderm",
               "d4_YS mesoderm","d6_YS mesoderm","d8_YS mesoderm","d10_YS mesoderm","d15_YS mesoderm",
               "d2_Amnion","d4_Amnion","d6_Amnion","d8_Amnion","d10_Amnion","d15_Amnion"
)
###ltq
bulk <- bulk.hed.merge.list[["bulk.hed.ltq"]]
gene <- top_genes_per_cluster_list[["top_hed_ltq"]]
top_genes <- intersect(gene,rownames(bulk))
bulk <- bulk[top_genes,]
bulk.hed.merge.list[["bulk.hed.ltq.top"]] <- bulk
rm(top_genes,gene)
df <- bulk
group1 <- df[, 1:4]
group2 <- df[, 5:112]
cor_matrix <- cor(group1, group2, use = "pairwise.complete.obs",method = "spearman")

cor_matrix <- cor_matrix[,hed_order] 
ltq_order <- c("preimplantation_ICM","preimplantation_Epiblast","preimplantation_Hypoblast",  
               "preimplantation_Trophoblast")
cor_matrix <- cor_matrix[ltq_order,]
pheatmap::pheatmap(cor_matrix,
                   cluster_rows = FALSE,  
                   cluster_cols = FALSE,  
                   display_numbers = FALSE, 
                   #gaps_row = c(4,23,36,55),
                   angle_col = 315,
                   cellwidth = 10,
                   cellheight = 10,
                   filename = "1024.heatmap.pre.png"
)
cor_matrix_list[["ltq"]] <- cor_matrix
rm(df,group1,group2,cor_matrix)

###CS7
bulk <- bulk.hed.merge.list[["bulk.hed.CS7"]]
gene <- top_genes_per_cluster_list[["top_hed_CS7"]]
top_genes <- intersect(gene,rownames(bulk))
bulk <- bulk[top_genes,]
bulk.hed.merge.list[["bulk.hed.CS7.top"]] <- bulk
rm(top_genes,gene)
df<- bulk
group1 <- df[, 1:19]
group2 <- df[, 20:127]
cor_matrix <- cor(group1, group2, use = "pairwise.complete.obs",method = "spearman")
cor_matrix <- cor_matrix[,hed_order] 
CS7_order <- c("CS7_PGC","CS7_Epiblast","CS7_Primitive Streak","CS7_Nascent Mesoderm",
               "CS7_Emergent Mesoderm","CS7_Advanced Mesoderm","CS7_Myeloid Progenitors",
               "CS7_Axial Mesoderm","CS7_Erythroblasts","CS7_Blood Progenitors",
               "CS7_Erythro-Myeloid Progenitors","CS7_Hemogenic Endothelium",
               "CS7_DE(NP)","CS7_DE(P)","CS7_YS Endoderm","CS7_YS Mesoderm",
               "CS7_Hypoblast","CS7_NNE","CS7_Amnion")
setdiff(CS7_order,rownames(cor_matrix))
cor_matrix <- cor_matrix[CS7_order,]
pheatmap::pheatmap(cor_matrix,
                   cluster_rows = FALSE,  
                   cluster_cols = FALSE,  
                   display_numbers = FALSE, 
                   #gaps_row = c(4,23,36,55),
                   angle_col = 315,
                   cellwidth = 10,
                   cellheight = 10,
                   filename = "1024.heatmap.CS7.png"
)

cor_matrix_list[["CS7"]] <- cor_matrix
rm(df,group1,group2,cor_matrix)

###CS8
bulk <- bulk.hed.merge.list[["bulk.hed.CS8"]]
gene <- top_genes_per_cluster_list[["top_hed_CS8"]]
top_genes <- intersect(gene,rownames(bulk))
bulk <- bulk[top_genes,]
bulk.hed.merge.list[["bulk.hed.CS8.top"]] <- bulk
rm(top_genes,gene)
df<- bulk
group1 <- df[, 1:13]
group2 <- df[, 14:121]
cor_matrix <- cor(group1, group2, use = "pairwise.complete.obs",method = "spearman")
cor_matrix <- cor_matrix[,hed_order] 
CS8_order <- c("CS8_Epi/Ecto","CS8_Gast/PS","CS8_Meso",
               "CS8_Ery","CS8_HEP",
               "CS8_Endo","CS8_Visceral Endo","CS8_YS Endo","CS8_YS EXM meso-A",
               "CS8_YS EXM meso-B","CS8_AM EXM Meso","CS8_AM","CS8_Noto")
cor_matrix <- cor_matrix[CS8_order,]

pheatmap::pheatmap(cor_matrix,
                   cluster_rows = FALSE,  
                   cluster_cols = FALSE,  
                   display_numbers = FALSE, 
                   #gaps_row = c(4,23,36,55),
                   angle_col = 315,
                   cellwidth = 10,
                   cellheight = 10,
                   filename = "1024.heatmap.CS8.png"
)

cor_matrix_list[["CS8"]] <- cor_matrix
rm(df,group1,group2,cor_matrix)

###CS10
bulk <- bulk.hed.merge.list[["bulk.hed.CS10"]]
gene <- top_genes_per_cluster_list[["top_hed_CS10"]]
top_genes <- intersect(gene,rownames(bulk))
bulk <- bulk[top_genes,]
bulk.hed.merge.list[["bulk.hed.CS10.top"]] <- bulk
rm(top_genes,gene)
df<- bulk
group1 <- df[, 1:19]
group2 <- df[, 20:127]
cor_matrix <- cor(group1, group2, use = "pairwise.complete.obs",method = "spearman")
cor_matrix <- cor_matrix[,hed_order] 

pheatmap::pheatmap(cor_matrix,
                   cluster_rows = FALSE,  
                   cluster_cols = FALSE,  
                   display_numbers = FALSE, 
                   #gaps_row = c(4,23,36,55),
                   angle_col = 315,
                   cellwidth = 10,
                   cellheight = 10,
                   filename = "1024.heatmap.CS10.png"
)
cor_matrix_list[["CS10"]] <- cor_matrix
rm(df,group1,group2,cor_matrix)

###CS12
bulk <- bulk.hed.merge.list[["bulk.hed.SWY"]]
gene <- top_genes_per_cluster_list[["top_hed_SWY"]]
top_genes <- intersect(gene,rownames(bulk))
bulk <- bulk[top_genes,]
bulk.hed.merge.list[["bulk.hed.SWY.top"]] <- bulk
rm(top_genes,gene)
df<- bulk
group1 <- df[, 1:56]
group2 <- df[, 57:164]
cor_matrix <- cor(group1, group2, use = "pairwise.complete.obs",method = "spearman")
cor_matrix <- cor_matrix[,hed_order] 
SWY_order <- c("CS12_PGC","CS12_somatic LPM","CS12_splanchnic LPM","CS12_somite","CS12_endothelium",
               "CS12_IM","CS12_limb","CS12_head mesoderm", "CS12_craniofacial",
               "CS12_miscellaneous","CS12_fibroblast",
               "CS12_blood",
               "CS12_epithelium",
               "CS12_endoderm","CS12_epidermis",       
               "CS12_neural progenitor","CS12_schwann", "CS12_neuron","CS12_sensory neuron",
               
               "CS13-14_PGC","CS13-14_somatic LPM","CS13-14_splanchnic LPM","CS13-14_somite","CS13-14_endothelium",
               "CS13-14_IM","CS13-14_limb","CS13-14_head mesoderm", "CS13-14_craniofacial",
               "CS13-14_miscellaneous","CS13-14_fibroblast",
               "CS13-14_blood",
               "CS13-14_epithelium",
               "CS13-14_endoderm","CS13-14_epidermis",       
               "CS13-14_neural progenitor","CS13-14_schwann", "CS13-14_neuron","CS13-14_sensory neuron",
               
               "CS15-16_somatic LPM","CS15-16_splanchnic LPM","CS15-16_somite","CS15-16_endothelium",
               "CS15-16_IM","CS15-16_limb","CS15-16_head mesoderm", "CS15-16_craniofacial",
               "CS15-16_miscellaneous","CS15-16_fibroblast",
               "CS15-16_blood",
               "CS15-16_epithelium",
               "CS15-16_endoderm","CS15-16_epidermis",       
               "CS15-16_neural progenitor","CS15-16_schwann", "CS15-16_neuron","CS15-16_sensory neuron"
)
cor_matrix <- cor_matrix[SWY_order,]

pheatmap::pheatmap(cor_matrix,
                   cluster_rows = FALSE,  
                   cluster_cols = FALSE,  
                   display_numbers = FALSE, 
                   #gaps_row = c(4,23,36,55),
                   angle_col = 315,
                   cellwidth = 10,
                   cellheight = 10,
                   filename = "1024.heatmap.SWY.png"
)
cor_matrix_list[["SWY"]] <- cor_matrix
rm(df,group1,group2,cor_matrix)

mtx.all <- rbind(cor.matrix.ltq,cor.matrix.CS7,cor.matrix.CS8,cor.matrix.CS10,cor.matrix.CS12)

## CS7/8 correlation
#处理bulk文件
bulk.CS7
bulk.CS8
genename.list

gene.CS7 <- genename.list[["CS7_ltq"]]
gene.CS8 <- genename.list[["gene.CS8"]]

gene.CS7.CS8 <- merge(gene.CS7,gene.CS8,by = "x")
genename.list[["CS7_CS8_ltq"]] <- gene.CS7.CS8
saveRDS(genename.list,)

bulk.CS7.cut <- bulk.CS7[gene.CS7.CS8$CS7,] 
rownames(bulk.CS7.cut) <-  gene.CS7.CS8$CS8

bulk.CS8.cut <- bulk.CS8[gene.CS7.CS8$CS8,]

colnames(bulk.CS7.cut) <- paste("CS7",colnames(bulk.CS7.cut),sep = "_")
colnames(bulk.CS8.cut) <- paste("CS8",colnames(bulk.CS8.cut),sep = "_")

#top_gene
library(dplyr)
top_genes_per_cluster_list
top_CS7 <- top_genes_per_cluster_list$top_CS7
top_CS8 <- top_genes_per_cluster_list$top_CS8
top_genes_per_cluster <- c(top_CS7,top_CS8)
top_genes_per_cluster <- unique(top_genes_per_cluster)
top_genes_per_cluster_list[["top_CS7_CS8"]] <- top_genes_per_cluster
rm(top_genes_per_cluster_1,top_genes_per_cluster_2,top_genes_per_cluster)
saveRDS(top_genes_per_cluster_list,"1226top_genes_per_cluster_list.rds")

#heatmap
bulk <- cbind(bulk.CS7.cut,bulk.CS8.cut)
gene <- top_genes_per_cluster_list[["top_CS7_CS8"]]
top_genes <- intersect(gene,rownames(bulk))
bulk <- bulk[top_genes,]
df<- bulk
group1 <- df[, 1:19]
group2 <- df[, 20:32]
cor_matrix <- cor(group1, group2, use = "pairwise.complete.obs",method = "spearman")
CS8_order <- c("CS8_Epi/Ecto","CS8_Gast/PS","CS8_Meso",
               "CS8_Ery","CS8_HEP",
               "CS8_Endo","CS8_Visceral Endo","CS8_YS Endo","CS8_YS EXM meso-A",
               "CS8_YS EXM meso-B","CS8_AM EXM Meso","CS8_AM","CS8_Noto")
cor_matrix <- cor_matrix[,CS8_order]
CS7_order <- c("CS7_PGC","CS7_Epiblast","CS7_Primitive Streak","CS7_Nascent Mesoderm","CS7_Emergent Mesoderm",
               "CS7_Advanced Mesoderm","CS7_Myeloid Progenitors","CS7_Axial Mesoderm","CS7_Erythroblasts",
               "CS7_Blood Progenitors","CS7_Erythro-Myeloid Progenitors","CS7_Hemogenic Endothelium",
               "CS7_DE(NP)","CS7_DE(P)","CS7_YS Endoderm","CS7_YS Mesoderm","CS7_Hypoblast","CS7_NNE","CS7_Amnion" )
cor_matrix <- t(cor_matrix)
cor_matrix <- cor_matrix[,CS7_order]
cor_matrix <- t(cor_matrix)
pheatmap::pheatmap(cor_matrix,
                   cluster_rows = FALSE,  
                   cluster_cols = FALSE,  
                   display_numbers = FALSE, 
                   #gaps_row = c(4,23,36,55),
                   angle_col = 315,
                   cellwidth = 20,
                   cellheight = 20,
                   filename = "1126.heatmap.CS7.CS8.topDEG.pdf"
)

cor_matrix_list[["CS8"]] <- cor_matrix
rm(df,group1,group2,cor_matrix)

bulk.CS7.CS8.all <- cbind(bulk.CS7.cut,bulk.CS8.cut)

df <- bulk.CS7.CS8.all
group1 <- df[, 1:19]
group2 <- df[, 20:32]
cor_matrix <- cor(group1, group2, use = "pairwise.complete.obs",method = "spearman")
CS8_order <- c("CS8_Epi/Ecto","CS8_Gast/PS","CS8_Meso",
               "CS8_Ery","CS8_HEP",
               "CS8_Endo","CS8_Visceral Endo","CS8_YS Endo","CS8_YS EXM meso-A",
               "CS8_YS EXM meso-B","CS8_AM EXM Meso","CS8_AM","CS8_Noto")
cor_matrix <- cor_matrix[,CS8_order]
CS7_order <- c("CS7_PGC","CS7_Epiblast","CS7_Primitive Streak","CS7_Nascent Mesoderm","CS7_Emergent Mesoderm",
               "CS7_Advanced Mesoderm","CS7_Myeloid Progenitors","CS7_Axial Mesoderm","CS7_Erythroblasts",
               "CS7_Blood Progenitors","CS7_Erythro-Myeloid Progenitors","CS7_Hemogenic Endothelium",
               "CS7_DE(NP)","CS7_DE(P)","CS7_YS Endoderm","CS7_YS Mesoderm","CS7_Hypoblast","CS7_NNE","CS7_Amnion" )
cor_matrix <- t(cor_matrix)
cor_matrix <- cor_matrix[,CS7_order]

cor_matrix <- t(cor_matrix)
pheatmap::pheatmap(cor_matrix,
                   cluster_rows = FALSE,  
                   cluster_cols = FALSE,  
                   display_numbers = FALSE, 
                   #gaps_row = c(4,23,36,55),
                   angle_col = 315,
                   cellwidth = 20,
                   cellheight = 20,
                   filename = "1126.heatmap.CS7.CS8.all.pdf"
)


