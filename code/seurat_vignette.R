.libPaths(c('/n/app/bcbio/R4.3.1_singlecell_dev/', '/n/app/bcbio/R4.3.1_singlecell'))
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(bcbioR)
options(future.globals.maxSize= 891289600)

# localdir <- "../final/10x_example_data/"
localdir <- "final/outs_test/"
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))

Assays(object)
DefaultAssay(object) <- "Spatial.008um"


vln.plot <- VlnPlot(object, features = "nCount_Spatial.008um", pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
count.plot <- SpatialFeaturePlot(object, features = "nCount_Spatial.008um") + theme(legend.position = "right")

# note that many spots have very few counts, in-part
# due to low cellular density in certain tissue regions
vln.plot | count.plot


DefaultAssay(object) <- "Spatial.008um"
object <- NormalizeData(object)

DefaultAssay(object) <- "Spatial.016um"
object <- NormalizeData(object)


DefaultAssay(object) <- "Spatial.008um"
object <- FindVariableFeatures(object)
object <- ScaleData(object)
# we select 50,0000 cells and create a new 'sketch' assay
object <- SketchData(
  object = object,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(object) <- "sketch"
# perform clustering workflow
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, assay = "sketch", reduction.name = "pca.sketch")
object <- FindNeighbors(object, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
object <- FindClusters(object, cluster.name = "seurat_cluster.sketched", resolution = .5)
object <- RunUMAP(object, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)

object <- ProjectData(
  object = object,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)

DefaultAssay(object) <- "sketch"
Idents(object) <- "seurat_cluster.sketched"
p1 <- DimPlot(object, reduction = "umap.sketch", label = F) + ggtitle("Sketched clustering (50,000 cells)") + theme(legend.position = "bottom")

# switch to full dataset
DefaultAssay(object) <- "Spatial.008um"
Idents(object) <- "seurat_cluster.projected"
p2 <- DimPlot(object, reduction = "full.umap.sketch", label = F) + ggtitle("Projected clustering (full dataset)") + theme(legend.position = "bottom")

p1 | p2

SpatialDimPlot(object) & scale_fill_cb_friendly()#, label = T, repel = T, label.size = 4)
