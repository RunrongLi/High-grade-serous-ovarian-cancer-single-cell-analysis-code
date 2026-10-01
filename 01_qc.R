# 01_qc.R
# Quality control and cell filtering for GSE165897 HGSOC scRNA-seq data
# Input : data/GSE165897_UMIcounts_HGSOC.tsv.gz
#         data/GSE165897_cellInfo_HGSOC.tsv.gz
# Output: output/01_GSE165897_qc.rds

suppressMessages(library(Seurat))
suppressMessages(library(Matrix))
# ---------------------------------------------------------------
# 1. Load count matrix and cell metadata
# ---------------------------------------------------------------
mtx <- read.delim("data/GSE165897_UMIcounts_HGSOC.tsv.gz")
cl  <- read.delim("data/GSE165897_cellInfo_HGSOC.tsv.gz")
rownames(mtx) <- mtx[, 1]
mtx1 <- mtx[, -1]

# ---------------------------------------------------------------
# 2. Create Seurat object
# ---------------------------------------------------------------
GSE165897 <- CreateSeuratObject(
  counts       = mtx1,
  project      = "GSE165897"
)
GSE165897$sample<- cl$sample
# ---------------------------------------------------------------
# 3. Filter cells
# ---------------------------------------------------------------
GSE165897[["percent.mt"]] <- PercentageFeatureSet(GSE165897, pattern = "^MT-")
GSE165897 <- subset(
  GSE165897,
  subset = nCount_RNA   > 500 &
    nFeature_RNA > 1000 &
    percent.mt   < 20
)

# ---------------------------------------------------------------
# 4. Save
# ---------------------------------------------------------------
dir.create("output", showWarnings = FALSE)
saveRDS(GSE165897, "output/01_GSE165897_qc.rds")