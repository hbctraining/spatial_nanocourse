## Loupe Browser

[Loupe Browser](https://www.10xgenomics.com/support/software/loupe-browser/latest) is a  visualization software from 10X genomics that allows you to explore and analyze your 10x Genomics Chromium and Visium data. You can also convert your Seurat objects into Loupe Browser files using the LoupeR package.

### Load in the data
We will use the cloupe file created by spaceranger, as you will notice this the full image and not the cropped image we have been using so far in class. By default, Space Ranger count produces a cloupe file for Visium HD Spatial Gene Expression data using 8 µm x 8 µm bins (can change this using the `--custom-bin-size` option).

### View panel

The workspace is centered around the **View Panel**, showing an H&E brightfield image of the tissue section overlaid with spots representing Visium Spatial Gene Expression data. Zoom in so the image fills the View space.
* Note the scale bar in the bottom right; set to 2mm by default
* Spot opacity is controlled with slider, making individual spots more or less visible.
* Projection settings allow the user to change individual settings, and then reset it back to original settings
* Different projection types (UMAP, Feature plot), but we will stick to Spatial

### Tool selector
On the left hand side you have the **Tool selector**:

**Clusters**

* Default is set to Clusters - these are clusters determined by spaceranger (similar to what we observed in the spaceranger report)
  * The graph-based clustering algorithm consists of building a sparse nearest-neighbor graph, followed by Louvain Modularity Optimization - most similar to what is performed in Seurat
  * k-means clustering, standard k-means where you need to specify the number k clusters you want
  * Use the "Create new group" to upload a csv of your own cluster labels. First column, cellular barcodes and second column is cell annotation. The first column's contents must match at least a subset of the barcodes in the .cloupe file

1. Hover over image, see the different clusters
2. De-select all and pick and choose select clusters (e.g. Cluster 16)
3. See the Differential expression output table below and genes that have the highest FC

**Features**

* Use the search box to look for the Prox1 gene; from the Allen Brain Atlas is observed to increased expression in Hippocampus. Defaults to log2 scale. Can change to linear scale. Turn on "Filter barcodes" and set the minimum to 1. This will omit bins with zero counts and help see where the gene is expressed exclusively on the tissue slide
* Try the same with Crlf1
* Try the same with Pkp2. See that different parts of the hippocampus are illuminated with each gene.
* Save this list of genes by "Edit list name" - Hippocampus
* Click combine all below
* Go back and select all clusters, then take a look at expression distribution

**Advanced selection**

* Allows you to create rules. Do this using the same three hippocampal genes.
* See how mnay barcodes you are left with when you use AND versus OR

 **Reanalyze**

 * A new box will pop up to guide you.
 * Threshold by features, UMIs and then choose whether you want a new spatial plot and/or UMAP, can even tweak some parameters ("advanced")
