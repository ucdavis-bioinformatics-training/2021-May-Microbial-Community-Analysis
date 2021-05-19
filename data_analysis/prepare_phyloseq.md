# Prepare R and copy over the Phyloseq object for Data Analysis

## Create a new RStudio project

Open RStudio and create a new project, for more info on project see [this page](<https://support.rstudio.com/hc/en-us/articles/200526207-Using-Projects>)

* File > New Project > New Directory > New Project (name the new directory, Ex. mca_analysis)

We first need to make sure we have the necessary packages, phyloseq, ggplot2, gridExtra, gridR, ape, and edgeR are installed (if not install it), and then load each package to veryify they installed correctly

To install the packages, In the R console run the following commands

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")

if (!any(rownames(installed.packages()) == "phyloseq")){
  BiocManager::install("phyloseq")
}
library(phyloseq); packageVersion("phyloseq")

if (!any(rownames(installed.packages()) == "biomformat")){
  BiocManager::install("biomformat")
}
library(biomformat); packageVersion("biomformat")

if (!any(rownames(installed.packages()) == "ggplot2")){
  BiocManager::install("ggplot2")
}
library(ggplot2); packageVersion("ggplot2")

if (!any(rownames(installed.packages()) == "gridExtra")){
  BiocManager::install("gridExtra")
}
library(gridExtra); packageVersion("gridExtra")

if (!any(rownames(installed.packages()) == "vegan")){
  BiocManager::install("vegan")
}
library(vegan); packageVersion("vegan")

if (!any(rownames(installed.packages()) == "edgeR")){
  BiocManager::install("edgeR")
}
library(edgeR); packageVersion("edgeR")

if (!any(rownames(installed.packages()) == "dada2")){
  BiocManager::install("dada2")
}
library(dada2); packageVersion("dada2")

if (!any(rownames(installed.packages()) == "phangorn")){
  BiocManager::install("phangorn")
}
library(phangorn); packageVersion("phangorn")

```

### Download the template Markdown workshop document and open it

In the R console run the following command

```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2021-May-Microbial-Community-Analysis/master/data_analysis/mca_part1.Rmd", "mca_part1.Rmd")
```

### Download the data file for the workshop document and preview/open it


```r
download.file("https://raw.githubusercontent.com/ucdavis-bioinformatics-training/2021-May-Microbial-Community-Analysis/master/data/phyloseq_nochim_silva.RData", "phyloseq_nochim_silva2.RData")
```

### Edit the file YAML portion

The top YAML (YAML ain't markup language) portion of the doc tells RStudio how to parse the document.

<pre><code>---
title: "Microbial Community Analysis in R"
author: your_name
date: current_date
output:
    html_notebook: default
    html_document: default
