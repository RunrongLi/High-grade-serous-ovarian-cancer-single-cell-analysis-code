suppressMessages(library(Seurat))
set.seed(123)
DefaultAssay(SeuObj.merge) <- "RNA"
Idents(SeuObj.merge) <- "celltype"
caf.markers <- FindMarkers(
  object           = SeuObj.merge,
  ident.1          = "CAFs",
  only.pos         = FALSE,
  min.pct          = 0.25,
  logfc.threshold  = 0.25,
  test.use = "wilcox"
)
