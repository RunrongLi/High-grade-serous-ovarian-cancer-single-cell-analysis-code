suppressMessages(library(Seurat))

set.seed(123)

# ---------------------------------------------------------------
# 1. Load QC-filtered object
# ---------------------------------------------------------------
GSE165897 <- readRDS("output/01_GSE165897_qc.rds")

# ---------------------------------------------------------------
# 2. Split by sample
# ---------------------------------------------------------------
pt.list <- SplitObject(GSE165897, split.by = "sample")

# ---------------------------------------------------------------
# 3. Per-sample normalization and variable feature selection
# ---------------------------------------------------------------
pt.list <- lapply(pt.list, function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(
    x,
    selection.method = "vst",
    nfeatures        = 2000,
    verbose          = FALSE
  )
  x
})

# ---------------------------------------------------------------
# 4. Select integration features
# ---------------------------------------------------------------
features <- SelectIntegrationFeatures(
  object.list = pt.list,
  nfeatures   = 2000
)

# ---------------------------------------------------------------
# 5. Per-sample scaling and PCA
# ---------------------------------------------------------------
pt.list <- lapply(pt.list, function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
  x
})

# ---------------------------------------------------------------
# 6. Find integration anchors (RPCA)
# ---------------------------------------------------------------
pt.anchors <- FindIntegrationAnchors(
  object.list     = pt.list,
  anchor.features = features,
  reduction       = "rpca"
)

# ---------------------------------------------------------------
# 7. Integrate
# ---------------------------------------------------------------
SeuObj.merge <- IntegrateData(anchorset = pt.anchors)
DefaultAssay(SeuObj.merge) <- "integrated"

# ---------------------------------------------------------------
# 8. Downstream dimensionality reduction and clustering
# ---------------------------------------------------------------
SeuObj.merge <- ScaleData(SeuObj.merge, verbose = FALSE)
SeuObj.merge <- RunPCA(SeuObj.merge, npcs = 30, verbose = FALSE)
SeuObj.merge <- RunUMAP(SeuObj.merge, reduction = "pca", dims = 1:30)
SeuObj.merge <- RunTSNE(SeuObj.merge, reduction = "pca", dims = 1:30)
SeuObj.merge <- FindNeighbors(SeuObj.merge, reduction = "pca", dims = 1:30)
SeuObj.merge <- FindClusters(SeuObj.merge, resolution = 0.5)
DefaultAssay(SeuObj.merge) <- "RNA"
# ---------------------------------------------------------------
# 9. Save
# ---------------------------------------------------------------
saveRDS(SeuObj.merge, "output/02_GSE165897_integrated.rds")