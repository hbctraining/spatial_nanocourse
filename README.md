## Spatial transcriptomics Nanocourse

| Audience | Computational skills required| Duration |
:----------|:----------|:----------|
| Biologists | [Introduction to R](https://hbctraining.github.io/Intro-to-R-flipped/) | 1-session in-person (~1h 45 min of trainer-led time)|


### Description
This repository contains materials for a module which is part of a Nanocourse organized by the [Single Cell Core at HMS](https://singlecellcore.hms.harvard.edu/). The nanocourse titled "Spatial Transcriptomics: Key Technologies, Experimental Considerations, and Data Analysis" introduces the fundamentals of data analysis for spatial transcriptomics, including common techniques and tools. **In this module we walk through the analysis workflow for Visium HD data analysis.**

### Learning Objectives
- LO1 ...
- ...

### Dataset
The dataset which is used in the pre-reading activities can be found at the link below.

* [Visium HD Dataset](https://github.com/hbctraining/BMI713_Intro_to_HPC/blob/master/data/unix_lesson.zip?raw=true)

### Lessons

[Workshop schedule (trainer-led learning)]()


### Installation Requirements

#### Applications
Download the most recent versions of R and RStudio for your laptop:

 - [R](http://lib.stat.cmu.edu/R/CRAN/) **(version 4.0.0 or above)**
 - [RStudio](https://www.rstudio.com/products/rstudio/download/#download)

#### Packages for R

> **Note 1: Install the packages in the order listed below.**

> **Note 2:  All the package names listed below are case sensitive!**

> **Note 3**: If you have a Mac with an M1 chip, download and install this tool before intalling your packages: https://mac.r-project.org/tools/gfortran-12.2-universal.pkg

> **Note 4**: At any point (especially if you’ve used R/Bioconductor in the past), in the console **R may ask you if you want to update any old packages by asking Update all/some/none? [a/s/n]:**. If you see this, **type "a" at the prompt and hit Enter** to update any old packages. _Updating packages can sometimes take quite a bit of time to run, so please account for that before you start with these installations._  

> **Note 5:** If you see a message in your console along the lines of “binary version available but the source version is later”, followed by a question, **“Do you want to install from sources the package which needs compilation? y/n”, type n for no, and hit enter**.


**(1)** Install the 8 packages listed below from **CRAN** using the `install.packages()` function. 

1. `tidyverse`
1. `BiocManager`
1. `Seurat`
1. `metap`

**Please install them one-by-one as follows:**

```r
install.packages("tidyverse")
install.packages("BiocManager")
install.packages("Seurat")
& so on ...
```

**(2)** Install the 4 packages listed below from **Bioconductor** using the the `BiocManager::install()` function.

1. `BiocPackage_name`


**Please install them one-by-one as follows:**

```r
BiocManager::install("BiocPackage_name")
BiocManager::install("BiocPackage_name")
& so on ...
```

**(3)** Finally, please check that all the packages were installed successfully by **loading them one at a time** using the `library()` function.  

```r
library(Seurat)
library(tidyverse)
...
```

**(4)** Once all packages have been loaded, run sessionInfo().  

```r
sessionInfo()
```



---

*These materials have been developed by members of the teaching team at the [Harvard Chan Bioinformatics Core (HBC)](http://bioinformatics.sph.harvard.edu/) RRID:SCR_025373. These are open access materials distributed under the terms of the [Creative Commons Attribution license](https://creativecommons.org/licenses/by/4.0/) (CC BY 4.0), which permits unrestricted use, distribution, and reproduction in any medium, provided the original author and source are credited.*
