library(Seurat)
library(qs)

options(future.globals.maxSize= 2000000000)

localdir <- '../path/to/spaceranger/outs/'

# list all samples available in spaceranger output (may need to adjust string-parsing)
samples <- list.files(localdir)[grepl('LIB', list.files(localdir))]
names(samples) <- unlist(transpose(strsplit(samples, '_'))[3])

# one by one, load samples in to seurat objects
# note that to run this you will need to have the hdf5r package installed
for (i in 1:length(samples)){
  object_sample <- Load10X_Spatial(data.dir = paste0(localdir, '/', samples[i], '/', 'spaceranger/outs/'),
                                   slice = names(samples)[i],
                                   filename = 'raw_feature_bc_matrix.h5',
                                   bin.size = 16)  # need to read in one bin at a time for multiple samples
  object_sample$orig.ident <- names(samples)[i]
  assign(paste0(samples[i], "_seurat"),
         object_sample) # stores Seurat object in variable of corresponding sample name
}
seurat_ID <- paste0(samples, "_seurat") # get names of all objects

###### for 2 samples #######

# merge 2 seurat objects into 1 seurat object
object <- merge(x = get(seurat_ID[1]),
                y = get(seurat_ID[2]),
                add.cell.ids = samples,
                project = "mouse_brain_visiumhd")

###### for more than 2 samples #######

rest_of_samples <- get(seurat_ID[2])
for (i in 3:length(seurat_ID)) {
  rest_of_samples <- c(rest_of_samples, get(seurat_ID[i]))
} ## makes a list of all seurat objects

# merge more than 2 seurat objects into 1 seurat object
object <- merge(x = get(seurat_ID[1]),
                y = rest_of_samples,
                add.cell.ids = samples,
                project = "mouse_brain_visiumhd")

qsave(object, 'seurat_object.qs')