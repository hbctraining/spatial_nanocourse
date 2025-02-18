library(Seurat)
library(qs)

options(future.globals.maxSize= 2000000000)

localdir <- '../path/to/spaceranger/outs/'

# to load raw feature matrix
object <- Load10X_Spatial(data.dir = localdir,
                          filename = 'raw_feature_bc_matrix.h5',
                          bin.size = 16)

# use plot with interactive=T to determine what coordinates define the region of interest
# SpatialDimPlot(object, interactive = T)

cropped.coords <- Crop(object[["slice1.016um"]], x = c(500, 1200), y = c(800, 1600), coords = "plot")
object[["zoom"]] <- cropped.coords

# create new object containing only cells in subset
object_subset <- subset(object, cells = object@images$zoom$centroids@cells)
object_subset[["zoom"]] <- NULL

qsave(object_subset, 'subsetted_seurat_object.qs')