# 03_visualization.R
# Marker gene visualization for GSE165897 integrated object
# Input : output/02_GSE165897_integrated.rds
# Output: output/03_marker_vlnplot.pdf

suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(purrr))
suppressMessages(library(grid))
# ---------------------------------------------------------------
# 1. Load integrated object
# ---------------------------------------------------------------
SeuObj.merge <- readRDS("output/02_GSE165897_integrated.rds")
DefaultAssay(SeuObj.merge) <- "RNA"

# ---------------------------------------------------------------
# 2. Helper functions
# ---------------------------------------------------------------
modify_vlnplot <- function(obj, feature, pt.size = 0,
                           plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  p <- VlnPlot(obj, features = feature, pt.size = pt.size, ...) +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(
      legend.position  = "none",
      axis.text.x      = element_blank(),
      axis.text.y      = element_blank(),
      axis.ticks.x     = element_blank(),
      axis.ticks.y     = element_line(),
      axis.title.y     = element_text(size = 45, angle = 0, vjust = 0.5),
      plot.margin      = plot.margin
    )
  return(p)
}

StackedVlnPlot <- function(obj, features, pt.size = 0,
                           plot.margin = unit(c(-0.75, 0, -0.75, 0), "cm"), ...) {
  plot_list <- purrr::map(features, function(x) {
    modify_vlnplot(obj = obj, feature = x, pt.size = pt.size,
                   plot.margin = plot.margin, ...)
  })
  
  plot_list[[length(plot_list)]] <- plot_list[[length(plot_list)]] +
    theme(
      axis.text.x  = element_text(size = 45, angle = 60, hjust = 1, vjust = 1),
      axis.ticks.x = element_line()
    )
  
  p <- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}

# ---------------------------------------------------------------
# 3. Plot marker genes
# ---------------------------------------------------------------
markers <- c(
  "EPCAM", "KRT8", "KRT19",     
  "CD2", "CD3D", "CD3E",        
  "KLRF1", "KLRC1", "KLRD1",    
  "CD19", "CD79A", "IGHG1",     
  "CD14", "AIF1",               
  "PECAM1", "VWF",              
  "COL1A1", "ACTA2", "DCN"      
)

p <- StackedVlnPlot(SeuObj.merge, markers, pt.size = 0)

# ---------------------------------------------------------------
# 4. Save
# ---------------------------------------------------------------
dir.create("output", showWarnings = FALSE)
ggsave("output/03_marker_vlnplot.pdf", p, width=28,height = 18)

#Cell type annotation is performed manually based on marker expression.