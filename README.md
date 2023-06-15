# scAnalysisR 
A package to perform single cell analysis using various tools. The package centers around using `Seurat` objects and has functions for preprocessing, plotting, running differential expression and gene ontology, and naming popultions.

This package works with `Seurat` version 3 objects.

If you use this package for your own analysis, please acknowledge this repository. You can also consider citing my eLife paper. The package I wrote for that paper was updated to work with `Seurat` version 3 objects here.

```
Wells, K. L., Miller, C. N., Gschwind, A. R., Wei, W., Phipps, J. D., Anderson, M. S., & Steinmetz, L. M. (2020). Combined transient ablation and single cell rna sequencing reveals the development of medullary thymic epithelial cells. ELife, 9, 1–80. https://doi.org/10.7554/eLife.60188
```

## Installation

Install `clustifyr`

```R
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("clustifyr")
```

Install devtools

```R
install.packages("devtools")
```

Install `scAnalysisR`

```R
library(devtools)
install_github("CUAnschutzBDC/scAnalysisR")
```

## Usage

This package is meant as a wrapper for many functions associated with single cell Analysis. It performs functions from loading in the data, processing the data, running differential expression, running gene ontology, and plotting.

### Loading in the data
Data that has been processed with `cellranger` can be read in using the `create_seurat_object` function. This has a few important arguments

* `sample` - The name of the sample. This should match the sample name that you provided to cellranger as it is used to identify the path of the cellranger output. This name will also be added as the `orig.ident` in the seurat object.
* `count_path` - If you ran `cellranger`, the path to the output directory you gave cellranger. Otherwise, this path can be more specific. See below for details on the path
* `ADT` - Were ADTs included in this run? If so, set this to `TRUE`, if not set this to `FALSE`. If ADTs were included you must set this to `TRUE` because the input data from cellranger looks different if multiple assays were included.
* `hashtag` - Were HTOs included? If so, set this to `TRUE`.
* `min_features` - What minimum features should be passed to `CreateSeuratObject`. Default is 200
* `min_cells` - What minimum cells should be passed to `CreateSeuratObject`. Default is 3.
* `hashtag_ident` - What are the identities of the hashtags that you passed in your antibody file to cellranger? This should be a list of hashtags ie `c("HTO1", "HTO2", "HTO3")`. This is used to pull out the hashtags and put them into a separate assay called `HTO`.
* `tenx_structure` - What cellranger was run? Options are `count`, `multi`, `multi7`, or `none`. This tells how to find the input data.

**Path to data**
The path to the data is defined by:
* `sample`
* `count_path`
* `tenx_structure`

If you ran `cellranger count` on your data and gave it `results` as your output directory and `sample1` as your sample name, the files we need are in `results/sample1/outs/filtered_feature_bc_matrix`. To tell this funtion how to find this data, you just need:

```R
so <- create_seurat_object(sample = "sample1", count_path = "results", tenx_structure = "count")
```

If you ran `cellranger multi` on your data and gave it `results` as your output directory and `sample1` as your sample name, the files we need are in `results/sample1/outs/count/filtered_feature_bc_matrix`. To tell this funtion how to find this data, you just need:

```R
so <- create_seurat_object(sample = "sample1", count_path = "results", tenx_structure = "multi")
```

If you ran `cellranger multi` with cellranger 7 or above on your data and gave it `results` as your output directory and `sample1` as your sample name, the files we need are in `results/sample1/outs/per_sample_outs/sample1/count/filtered_feature_bc_matrix`. To tell this funtion how to find this data, you just need:

```R
so <- create_seurat_object(sample = "sample1", count_path = "results", tenx_structure = "multi7")
```

The goal of this function is to very easily find your files with minimal input from you. If there is a problem finding the files, the error message will print the path of the directory where it is looking for the files so you can more easily debug. In the case that your output files are not in the above configuration, you can provide the path to the directory that contains the `matrix.mtx`, `features.tsv`, and `barcodes.tsv` files and have `tenx_structure = "none"`


```R
so <- create_seurat_object(sample = "sample1", count_path = "full/path/to/filtered_feature_bc_matrix",  
                           tenx_structure = "none")
```

**HTOs and ADTs**
If you included any antibody capture elements in your cellranger run, you will need to set either `HTO` or `ADT` or both to `TRUE`. This is because the data will be read in as a list and the function will fail if it does not know to look for this.

* If `ADT` is `TRUE` and `HTO` is `FALSE`, your Seurat object will be returned with an `ADT` assay with all of the features from your antibody file as the features. Your seurat object will also have an `RNA` assay that contains the gene expression data.
* If `ADT` is `TRUE` and `HTO` is `TRUE`, your Seurat object will be returned with an `ADT` assay with all of the features from your antibody file except for those listed in `hashtag_idents` as the features. It will also have a `HTO` assay will all of the features from the `hashtag_idents` as your features. If you don't set `hashtag_idents`, it will search for any features that have `hashtag` in the name. Your seurat object will also have an `RNA` assay that contains the gene expression data.
* If `ADT` is `FALSE` and `HTO` is `TRUE`, your Seurat object will be returned with `HTO` assay will all of the features from the `hashtag_idents` as your features. If you don't set `hashtag_idents`, it will search for any features that have `hashtag` in the name. Your seurat object will also have an `RNA` assay that contains the gene expression data.
