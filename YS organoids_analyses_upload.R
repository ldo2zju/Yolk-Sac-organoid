####pre-processing single-cell RNA-seq data for YS organoids
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

###RNA velocity analyses
library(stringr)
Bmp4_YS_D6_velocyto <- read.loom.matrices("/path/Bmp4_YS_D6.loom", engine = "hdf5r")
Bmp4_YS_D6 <- readRDS("/path/BMP4_YS_D6_updated_annotation.rds")

cell_annotation_abbre <- c("Epi","T-Epi","PS","Nas Mes",
                           "Eme Mes","Adv Mes/Cardio",
                           "BV","T-Def Endo",
                           "Def Endo","Endo-gut",
                           "YS Meso prog","YS Endo prog",
                           "Amnion","PGC","Mig PGC")


Bmp4_YS_D6 <- SetIdent(Bmp4_YS_D6,value = "update_celltypes_ordered")
names(cell_annotation_abbre) <- levels(Bmp4_YS_D6)
Bmp4_YS_D6 <- RenameIdents(Bmp4_YS_D6,cell_annotation_abbre)
Bmp4_YS_D6@meta.data$cell_annotation_abbre <- as.character(Idents(Bmp4_YS_D6))
saveRDS(Bmp4_YS_D6,"/path/BMP4_YS_D6_updated_annotation.rds")

# spliced
Bmp4_YS_D6_spliced <- Bmp4_YS_D6_velocyto[["spliced"]]
colnames(Bmp4_YS_D6_spliced) <- colnames(Bmp4_YS_D6_spliced) %>% str_replace("Bmp4_YS_D6:", "")  %>% str_replace("x", "-1")

colnames_intersect <- intersect(rownames(Bmp4_YS_D6@meta.data), colnames(Bmp4_YS_D6_spliced))
Bmp4_YS_D6_spliced <- Bmp4_YS_D6_spliced[, which(colnames(Bmp4_YS_D6_spliced) %in% colnames_intersect)]
## integrate into seuratobj
Bmp4_YS_D6[["spliced"]] <- CreateAssayObject(counts = as.matrix(Bmp4_YS_D6_spliced))

# check rownames
rownames_intersect <- intersect(rownames(Bmp4_YS_D6_spliced), rownames(Bmp4_YS_D6@assays$RNA@counts))
Bmp4_YS_D6_intersect <- Bmp4_YS_D6[rownames_intersect, ]

# unspliced ###
Bmp4_YS_D6_unspliced <- Bmp4_YS_D6_velocyto[["unspliced"]]
colnames(Bmp4_YS_D6_unspliced) <- colnames(Bmp4_YS_D6_unspliced) %>% str_replace("Bmp4_YS_D6:", "")  %>% str_replace("x", "-1")

colnames_intersect <- intersect(rownames(Bmp4_YS_D6@meta.data), colnames(Bmp4_YS_D6_unspliced))
Bmp4_YS_D6_unspliced <- Bmp4_YS_D6_unspliced[, which(colnames(Bmp4_YS_D6_unspliced) %in% colnames_intersect)]
## integrate into seuratobj
Bmp4_YS_D6[["unspliced"]] <- CreateAssayObject(counts = as.matrix(Bmp4_YS_D6_unspliced))

# check rownames
rownames_intersect <- intersect(rownames(Bmp4_YS_D6_unspliced), rownames(Bmp4_YS_D6@assays$RNA@counts))
Bmp4_YS_D6_intersect <- Bmp4_YS_D6[rownames_intersect, ]


# cell type as character
Bmp4_YS_D6_intersect$update_celltypes <- as.character(Bmp4_YS_D6_intersect$update_celltypes)
Bmp4_YS_D6_intersect$update_celltypes_ordered <- as.character(Bmp4_YS_D6_intersect$update_celltypes_ordered)
Bmp4_YS_D6_intersect$cell_annotation_abbre <- as.character(Bmp4_YS_D6_intersect$cell_annotation_abbre)
# convert
colnames(Bmp4_YS_D6_intersect@meta.data)  <- colnames(Bmp4_YS_D6_intersect@meta.data) %>% str_replace_all("\\.", "_")
colnames(Bmp4_YS_D6_intersect@assays$RNA@meta.features) <- colnames(Bmp4_YS_D6_intersect@assays$RNA@meta.features) %>% str_replace_all("\\.", "_")
Bmp4_YS_D6_intersect@assays$RNA@meta.features[VariableFeatures(Bmp4_YS_D6_intersect), 'highly_variable_genes'] <- 1
Bmp4_YS_D6_intersect@assays$RNA@meta.features[is.na(Bmp4_YS_D6_intersect@assays$RNA@meta.features$highly_variable_genes), 'highly_variable_genes'] <- 0

library(reticulate)
conda_list() 
reticulate::use_condaenv("r-reticulate", required = TRUE)

.regularise_df <- function(df, drop_single_values = TRUE) {
  if (ncol(df) == 0) df[['name']] <- rownames(df)
  if (drop_single_values) {
    k_singular <- sapply(df, function(x) length(unique(x)) == 1)
    if (sum(k_singular) > 0)
      warning(paste('Dropping single category variables:'),
              paste(colnames(df)[k_singular], collapse=', '))
    df <- df[, !k_singular, drop=F]
    if (ncol(df) == 0) df[['name']] <- rownames(df)
  }
  return(df)
}

seurat2anndata <- function(
    obj, outFile = NULL, slot = 'counts', main_layer = 'RNA', transfer_layers = c("spliced", "unspliced"), drop_single_values = TRUE
) {
  
  if (compareVersion(as.character(obj@version), '3.0.0') < 0)
    obj <- Seurat::UpdateSeuratObject(object = obj)
  
  X <- Seurat::GetAssayData(object = obj, assay = main_layer, slot = slot)
  
  obs <- .regularise_df(obj@meta.data, drop_single_values = drop_single_values)
  
  var <- .regularise_df(Seurat::GetAssay(obj, assay = main_layer)@meta.features, drop_single_values = drop_single_values)
  
  obsm <- NULL
  reductions <- names(obj@reductions)
  if (length(reductions) > 0) {
    obsm <- sapply(
      reductions,
      function(name) as.matrix(Seurat::Embeddings(obj, reduction=name)),
      simplify = FALSE
    )
    names(obsm) <- paste0('X_', tolower(names(obj@reductions)))
  }
  
  layers <- list()
  for (layer in transfer_layers) {
    mat <- Seurat::GetAssayData(object = obj, assay = layer, slot = slot)
    layers[[layer]] <- Matrix::t(mat)
  }
  
  anndata <- reticulate::import('anndata', convert = FALSE)
  
  adata <- anndata$AnnData(
    X = Matrix::t(X),
    obs = obs,
    var = var,
    obsm = obsm,
    layers = layers
  )
  
  if (!is.null(outFile))
    adata$write(outFile, compression = 'gzip')
  
  adata
}

# use seurat2ann script
seurat2anndata(obj = Bmp4_YS_D6_intersect, outFile = "/home/chengtao/single_cell_result/BMP4_YS/Resubmission/Bmp4_YS_D6_velocyto.h5ad")
##analyses on python

















