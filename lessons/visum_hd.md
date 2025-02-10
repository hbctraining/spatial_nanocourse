---
title: "Visium HD Analysis"
author: "Alex Bartlett, Meeta Mistry and Will Gammerdinger"
date: "February 4th, 2025"
---

Contributors: Alex Bartlett, Meeta Mistry and Will Gammerdinger

Approximate time: XX minutes

# Learning Objectives 

In this lesson, we will:
- Describe the elements of the Seurat object that are unique to spatial technologies (Learning Objective 1)
- Visually inspect and compare spatial scRNA-seq data before and after filtering (Learning Objective 2)
- Interact with the spatial seurat object to superimpose cluster onto the image (Learning Objective 3)

# NGS-based Spatial Transcriptomics Data Analysis

## Visium HD 

Each Visium HD slide has the same 6.5 x 6.5mm capture area as previous Visium products but is covered with about 10 million uniquely-barcoded oligonucleotide squares. These 2 micron tiles are arrayed in a continuous lawn across the entire capture area.

full coding transcriptome probes placed on a gap-less grid, barcoded in 2 μm square regions (bins) which are then grouped into 8 μm square bins for default analysis or annotation. 

The data is accompanied by a matching high resolution bright field (e.g. hematoxylin & eosin, H&E) or fluorescent (e.g. immuno-fluorescence, IF) morphology image. While the 8 μm resolution is a big improvement over original Visium’s ∼55 μm spots, having access to 2 μm bins along with matching morphology information makes it tempting to reconstruct single cells from the data.


### Preprocessing Data with Spaceranger

In the Visium HD assay, the barcodes are patterned in a continuous grid of 2x2 µm squares. By default, the Space Ranger pipeline creates 8x8 µm and 16x16 µm bins of gene expression data. For the purposes of this lesson, we will use the 16µm binning.


Sequencing facilities often output scRNAseq data, including spatial scRNAseq data, in FASTQ format. Because this is VisiumHD data from 10X genomics, we use their proprietary preprocessing software [Space Ranger](https://www.10xgenomics.com/support/software/space-ranger/latest) to process the FASTQ files into a count matrix and other images.

**TODO link to spaceranger report**

**Insert figure of workflow here**

# Setting up our environment

Before we start processing our data, we first need to set-up our R environment...

## Loading Libraries

We load the libraries necessary for processing scRNAseq data.

Run this `.libPaths()` command in order to...

```
# .libPaths(c('/n/app/bcbio/R4.3.1_singlecell_dev/', '/n/app/bcbio/R4.3.1_singlecell'))
.libPaths('../libs/')
```

Next, we will need to load these libraries using:

```
library(tidyverse)
library(Seurat)
library(patchwork)
library(qs)
```

> Note: While not necessary for this lesson, if you are running the full workflow, you may want to consider adding these commands to your R Script because...:
>
> ```
>options(future.globals.maxSize= 891289600)
>library(Azimuth)
>```




# Creating the Seurat Object

The Seurat object is a custom list-like object that has well-defined spaces to store specific information/data for single cell experiments, including spatial experiments and Visium HD.

**Insert picture of Surat object with spatial slots**

The Seurat package provides a function `Load10X_Spatial` to easily create a Seurat object from the output of Space Ranger.

In the Visium HD assay, the barcodes are patterned in a continuous grid of 2x2 µm squares. By default, the Space Ranger pipeline creates 8x8 µm and 16x16 µm bins of gene expression data. Our Seurat object will have data from both of these binnings, but for the purposes of this lesson, we will use the 8µm binning.

**What is this code exactly and why are we running it? Can we put it in a dropdown?**
```
# DO NOT RUN
localdir <- "../final/outs_test/"
object <- Load10X_Spatial(data.dir = localdir,
                          filename = 'raw_feature_bc_matrix.h5')
object <- Load10X_Spatial(data.dir = localdir)
DefaultAssay(object) <- "Spatial.008um"
```

We will read in a...

```
object <- qread('../data_processed/visiumhd_intestine_clustered.qs')
DefaultAssay(object) <- "Spatial.008um"
```

***

**Exercise**

This feels like a nice break for an exercise or even a question. I am not sure what yet because most of this feels like setting up the environment, but maybe have people explore the Seurat object a smidge? I know we are making an image for the scRNA-seq course on the Seurat object, maybe we could provide that and ask people what different slots look like in terms of data? This exercise should be reflective of Learning Objective 1, which might be like "Describe the elements of the Seurat object that are unique to spatial technologies"

***

# Quality Control

The main objective of quality control is to filter the data so that we include only true cells that are of high quality. This makes it so that when we cluster our cells, it is easier to identify distinct cell type populations.

Challenges include:

- **Delineating cells that are poor quality from less complex cells** - Elaborate
- **Choosing appropriate thresholds for filtering, so as to keep high quality cells without removing biologically relevant cell types** - Elaborate

**Maybe insert an image or two for these challenges demonstrating the challenge if possible**

Various metrics can be used to filter low-quality cells from high-quality ones, including:

- **UMI counts per bin** - Elaborate
- **Genes detected per bin** - Elaborate
- **Complexity (novelty score)** - Elaborate
- **Mitochondrial counts ratio** - Elaborate

Space Ranger applies filtering by default, and for this lesson, we will be working with a Seurat object loaded from Space Ranger's filtered data. However, we can compare plots of UMI counts per cell and genes detected per cell before and after Space Ranger's filtering. 

## Pre-filtering

In the figure below, we can see **XYZ** before filtering on...

<p align="center">
<img src="../img/preQC.png" width="600">
</p>

## Post-Filtering

Let's compare that with the filtered output from Space Ranger. First, we will need to create a metadata object using this command:

```
object_meta <- object@meta.data
```

Next, we will use this metadata object to create some plots to help us compare the filtered data to the pre-filtered data.

### Genes per Bin

In this first plot we can observe the number of genes per bin. We can use the following code to visualize this:

```
p1 <- object_meta %>%
   ggplot(aes(x=nCount_Spatial.008um)) +
   geom_density(alpha = 0.2) +
   scale_x_log10() +
   theme_classic() +
   ylab("Cell density") +
  ggtitle('PostQC Genes/Bin')
```

**Let's make a more descriptive name than p1**

This plot should look like:

<p align="center">
<img src="../img/genes_per_bin.png" width="600">
</p>

When we compare this to the pre-filtered data we can see **XYZ**. We expect to see **XYZ** and we can see that this filtering is capturing that. 

### UMIs per Bin

In this next plot we can observe the number of UMIs per bin. We can use the following code to visualize this:

```
p2 <- object_meta %>%
   ggplot(aes(x=nFeature_Spatial.008um)) +
   geom_density(alpha = 0.2) +
   scale_x_log10() +
   theme_classic() +
   ylab("Cell density") +
  ggtitle('PostQC UMIs/Bin')
```

**Let's make a more descriptive name than p2**

This plot should looke like:

```
<p align="center">
<img src="../img/UMIs_per_bin.png" width="600">
</p>
```

**Explain this code**

```
p3 <- p1 | p2
p3
```

**Let's make a more descriptive name than p3**


### Visualize Counts Data

**What is the scope of this section? Still post-filtering or is it general QC?**

We can visualize the number of counts per bin, both as a distribution and layered on top of the tissue image. Note that many spots have very few counts, in part due to low cellular density or cell types with low complexity in certain tissue regions.

First, we will plot... We are hoping to see...

```
vln.plot <- VlnPlot(object, features = "nCount_Spatial.008um", pt.size = 0) + NoLegend()
```

Next, we will plot... We are hoping to see... 

```
count.plot <- SpatialFeaturePlot(object, features = "nCount_Spatial.008um") + theme(legend.position = "right")
```

**What is this code doing?**

```
vln.plot | count.plot
```

***

**Exercise**

Perhaps have participants carry out the same QC plots for complexity and mitochondrial counts and get their impressions on the filtering. You may need to provide the code for creating these plots. This should be reflective of Learning Objective 2 and I am thinking it might be something like "Visually inspect and compare spatial scRNA-seq data before and after filtering"

***

# Normalize Data

Normalization is important in order to make expression counts comparable across genes and/or sample. Here we use a standard log-normalization for spatial data. We note that the best normalization methods for spatial data are still being developed and evaluated. Below we provide the code for doing this, but do not run this because XYZ...

```
# DO NOT RUN
object <- NormalizeData(object)
```

In the interest of computational time and memory resources, ee will provide you with the output of this normalization step below...

# Unsupervised Clustering

The authors of the Seurat package recommend the Seurat v5 sketch clustering workflow exhibits improved performance, especially for identifying rare and spatially restricted groups. Sketch-based analyses aim to ‘subsample’ large datasets in a way that preserves rare populations. Here, we sketch the Visium HD dataset, perform clustering on the subsampled cells, and then project the cluster labels back to the full dataset. Below we provide the code for doing this, but do not run this because it may take about 5 minutes to run.

## Sketch the Visium HD Dataset

This first step is doing... because...

**Can we get a figure of what this is doing?**

```
# DO NOT RUN
object <- FindVariableFeatures(object)
```

Next, we are doing... because...

```
# DO NOT RUN
object <- ScaleData(object)
```

**Can we get a figure or a sample table of what this is doing?**

Lastly, we are selecting 50,0000 cells and creating a new 'sketch' assay in order to... 

```
# DO NOT RUN
object <- SketchData(
  object = object,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)
```

**Break down each argument in the above command**

**Can we give participants an intermediate object here to explore?**

## Perform clustering on subsampled cells

This next step will take ~10 minutes to run, so we will provide you with this object at the end...

First we need to make out default assay in the Seurat object be `sketch` using:
```
# DO NOT RUN
# switch analysis to sketched cells
DefaultAssay(object) <- "sketch"
```

Similiarly to earlier, we will need to find variable features in order to...

```
# DO NOT RUN
object <- FindVariableFeatures(object)
```

We will also need to scale our data because...

```
# DO NOT RUN
object <- ScaleData(object)
```

**Could we provide the object here and have particpants do the next steps?**

Now we would like to perform a PC analysis in order to determine...

```
# DO NOT RUN?
object <- RunPCA(object, assay = "sketch", reduction.name = "pca.sketch")
```

**Break down each argument in the above command**

**Can we plot the PCA?**

Next, we need to find our nearest neighbors in order to...

```
# DO NOT RUN?
object <- FindNeighbors(object, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
```

**Break down each argument in the above command**

Now, we can find clusters...

```
# DO NOT RUN?
object <- FindClusters(object, cluster.name = "seurat_cluster.sketched", resolution = .5)
```

**Break down each argument in the above command**

```
# DO NOT RUN?
object <- RunUMAP(object, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)
```

**Break down each argument in the above command**

**Can we plot the UMAP?**


## Project cluster labels back to the full dataset

Now that we have our clusters from our subsampled dataset, we need to project these onto the full dataset. **Talk more about why?** This next step may take a few so, we will prvoide you with the completed object. But the code to do this would be:

```
# DO NOT RUN
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
```
**Break down each argument in the above command**

Because this is a large object, we wil use `qsave` from the `qs` package...

```
# DO NOT RUN
qsave(object, '../data_processed/visiumhd_intestine_clustered.qs')
```

### Visualizing the clusters

We can now visualize our clusters...

First we will need to set our default assay for the Seurat object to `sketch` using:

```
DefaultAssay(object) <- "sketch"
```

Next, we will assign the string of ...

```
Idents(object) <- "seurat_cluster.sketched"
```

Now we can create a `DimPlot` in order to show...

```
p1 <- DimPlot(object, reduction = "umap.sketch", label = F) +
        ggtitle("Sketched clustering") +
        theme(legend.position = "bottom")
```

**I think p1 already exists in this workflow and also can we call this something more informative?**

This figure should look like:

<p align="center">
<img src="../img/DimPlot.png" width="600">
</p>

From this `Dimplot` we can see that...

Now we will swtich our default assay to `patial.008um` switch to full dataset because...

```
DefaultAssay(object) <- "Spatial.008um"
```

We will need to change the string of... because...

```
Idents(object) <- "seurat_cluster.projected"
```

```
p2 <- DimPlot(object, reduction = "full.umap.sketch", label = F, raster = F)+
        ggtitle("Projected clustering") +
        theme(legend.position = "bottom")
```
**I think p2 already exists in this workflow and also can we call this something more informative?**

**Break down each argument in the above command**

**What is this code below doing?**

```
p1 | p2
```

### Visualizing Clusters on Image

In order to see the clusters superimposed on our image we need to...

```
p1 <- ImageDimPlot(object, cols = 'polychrome') 
```

**I think p1 already exists in this workflow and also can we call this something more informative?**

This figure should look like:

<p align="center">
<img src="../img/Image_DimPlot_1.png" width="600">
</p>

#### Zooming in on clusters of interest

In order to zoom in on an area of interest, we first need to crop the coordiantes using the `Crop` function in the `XYZ` package like:

```
cropped.coords <- Crop(object[["slice1.008um"]],
                       x = c(1550, 1750),
                       y = c(1250, 1450),
                       coords = "plot")
```

**Break down each argument in the above command**

Next, we need to assign these cropped coordiantes to the object canmed `zoom` in ...

```
object[["zoom"]] <- cropped.coords
```

Next, we need to ... in order to....

```
p2 <- ImageDimPlot(object,
                   fov = 'zoom',
                   cols = 'polychrome') 
```

**Break down each argument in the above command**

This figure should look like:

<p align="center">
<img src="../img/Image_DimPlot_2.png" width="600">
</p>

**What is this code doing?**

```
p1 | p2
```

***

**Exercise**

Zoom in on a different area of interest and this should reflect the goals of Learning Objective 3, which might be something like "Interact with the spatial seurat object to superimpose cluster onto the image"

***

# Cell Type Identification

Now that we have identified our desired clusters, we can move on to cell type identification, which will allow us to verify the identity of the cells contained in our various clusters.

[Azimuth](https://azimuth.hubmapconsortium.org/) is a web application that uses an annotated reference dataset to automate the processing, analysis, and interpretation of a new single-cell RNA-seq experiment. While we won't run it on this dataset because **XYZ**, you would be able to run it on your data using the following command:

```
# DO NOT RUN
object <- RunAzimuth(object, reference = "something")
```

**Can we show an image of these annotations on the figre that we've created?**

***

You have now completed your visualization of Visium HD ...

***

[Back to Schedule](../schedule/README.md)

***

*These materials have been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/). These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
