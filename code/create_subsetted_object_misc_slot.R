library(Seurat)
library(qs)

localdir <- "outs/"

# Load full spatial object
object <- Load10X_Spatial(
  data.dir = localdir,
  filename = "raw_feature_bc_matrix.h5",
  bin.size = 16
)

# Crop the image (this keeps all slots, including misc)
cropped <- Crop(
  object[["slice1.016um"]],
  x = c(500, 1200),
  y = c(800, 1600),
  coords = "plot"
)

# Attach cropped image BEFORE subsetting
object[["zoom"]] <- cropped

# Determine cells inside cropped area
cells_to_keep <- object[["zoom"]]$centroids@cells

# Subset the Seurat object
object_subset <- subset(object, cells = cells_to_keep)

# Make sure the subset retains the cropped image
object_subset[["zoom"]] <- cropped
object_subset
# An object of class Seurat 
# 32285 features across 43167 samples within 1 assay 
# Active assay: Spatial.016um (32285 features, 0 variable features)
# 1 layer present: counts
# 2 spatial fields of view present: slice1.016um zoom

object_subset@images[["slice1.016um"]] <- object_subset@images[["zoom"]]
object_subset@images[["zoom"]] <- NULL
object_subset
# An object of class Seurat 
# 32285 features across 43167 samples within 1 assay 
# Active assay: Spatial.016um (32285 features, 0 variable features)
# 1 layer present: counts
# 1 spatial field of view present: slice1.016um

SpatialFeaturePlot(object_subset, 
                   feature = 'nCount_Spatial.016um', 
                   pt.size.factor = 8)

# Save
qsave(object_subset, "MsBrain_FF-A1_subset_misc.qs")
