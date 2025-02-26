## Loupe Browser

Loupe Browser is a  visualization software from 10X genomics that allows you to explore and analyze your 10x Genomics Chromium and Visium data. You can also convert your Seurat objects into Loupe Browser files using the LoupeR package.

### Load in the data
We will use the cloupe file created by spaceranger, as you will notice this the full image and not the cropped image we have been using so far in class.


Let's start with the **user interface**:

The workspace is centered around the **View Panel**, showing an H&E brightfield image of the tissue section overlaid with spots representing Visium Spatial Gene Expression data. 
* Spot opacity is controlled with slider, making individual spots more or less visible.
* Projection settings allow the user to change individual settings, and then reset it back to original settings
* Differnt projection types (UMAP, Feature plot), but we will stick to Spatial
  
On the left hand side you have the **Tool selector**.
* Default is set to Clusters - these are clusters determined by spaceranger (simialry
