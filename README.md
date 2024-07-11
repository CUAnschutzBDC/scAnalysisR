# scAnalysisR 
Please cite using the zenodo [DOI: 10.5281/zenodo.12726713](https://doi.org/10.5281/zenodo.12726713)


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

There are a few steps that you will need to run outside of this package, for example filtering data based on quality, adding mitochondrial percents, and demultiplexing hashtags. If you want any functionality added here, feel free to submit an issue.

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

* If `ADT = TRUE` and `HTO = FALSE`, your Seurat object will be returned with an `ADT` assay with all of the features from your antibody file as the features. Your seurat object will also have an `RNA` assay that contains the gene expression data.
* If `ADT = TRUE` and `HTO = TRUE`, your Seurat object will be returned with an `ADT` assay with all of the features from your antibody file except for those listed in `hashtag_idents` as the features. It will also have a `HTO` assay will all of the features from the `hashtag_idents` as your features. If you don't set `hashtag_idents`, it will search for any features that have `hashtag` in the name. Your seurat object will also have an `RNA` assay that contains the gene expression data.
* If `ADT = FALSE` and `HTO = TRUE`, your Seurat object will be returned with `HTO` assay will all of the features from the `hashtag_idents` as your features. If you don't set `hashtag_idents`, it will search for any features that have `hashtag` in the name. Your seurat object will also have an `RNA` assay that contains the gene expression data.

### PCA and UMAP

#### PCA

To run PCA, use the function `PCA_dimRed`.

```R
so <- PCA_dimRed(so, assay = "RNA", reduction_name = NULL, vars_to_regress = NULL
```

This has only a handful of arguments

* `sample_object` - The seurat object loaded in by `create_seurat_object`
* `assay` - The assay to run `PCA` on. This can be any assay in your object, for example `RNA`, `ADT`, or `SCT`. If you are running on the `ADT` assay, if you have 30 or fewer ADTs you likely should just run `UMAP` directly on the ADT expression data and not perform PCA first.
* `reduction_name` - What the name of the reduction should be. These will be the names of the dimensional redution saved in the seurat object, so you can run this function on all assays without overriding the existing dimensional reductions. Many are provided but can be overridden by providing a avlue to `reduction_name`. The default values are:
  * If `assay = "RNA"`: `reduction_name = pca`
  * If `assay = "ADT"`: `reduction_name = apca`
  * If `assay = "SCT"`: `reduction_name = sctpca`
  * If `assay = integrated`: `reduction_name = pca`
* `vars_to_regress` - Only used if running on the `ADT` assay. These are used when scaling the ADT data before running PCA.

The return of this is a seurat object with a new dimensional reduction named based on the `reduction_name` provided

**Helpful plots**

To help evaluating the PCA, a plotting function is provided that returns a handful of qc plots.

```R
qc_plots <- plot_PCA(so, HTO = FALSE, ADT = FALSE,
                     assay = "RNA", jackstraw = TRUE,
                     reduction = NULL, data_type = "RNA",
                     ffpe = TRUE)
```

This function takes your seurat object and returns helpful qc plots as a named list
* The top genes assocaited with PC1 and PC2
* PC plot colored by the sample (orig.ident)
* PC plot colored by the percent mitochondrial reads
* PC plot colored by the number of RNA features
* PC plot colored by the number of RNA counts
* PC plot colored by the number of ADT counts (if `ADT = TRUE`)
* PC plot colored by the number of ADT feature (if `ADT = TRUE`)
* The HTO classification (if `HTO = TRUE`)
* A Jackstraw plot (if `jackstraw = TRUE`) to show how many PCs to use for downstream analysis.
* An Elbow plot to show how many PCs to use for downstream analysis

This function will assume your pca reduction is named based on the default values above. If your pca reduction is not based on the reductions above, you can supply your own to `reduction`.

The `ffpe` is only used if `data_type = spatial`. It is ignored if `data_type = "RNA"` 

*Note this function requires that you have identified the percent mitochondiral reads and saved them as a column "percent.mt" in your seurat object metadata. If you have HTOs, this also requires that you have a column callsed "HTO_classification" that identifies the different hashtagged samples.*

The arguments of this function are:
* `sample_object` A seurat object
* `HTO` OPTIONAL if HTOs were included in the seurat object. Default is FALSE
* `ADT` OPTIONAL If ADTs were included in the experiment. Default is FALSE
* `assay` OPTIONAL What assay to use. This is just used to locate the PCA reduction (the reduction name provided assumes that no reduction name was supplied to PCA_dimRed, if you did provide a different reduction name, you must use "reduction" below). Default is "RNA".
* `jackstraw` OPTIONAL if a jackstraw plot should be made. This can help in determining number of PCs to use, but it can be very slow. If TRUE, a jackstraw plot will only be made if the assay is also "RNA". Default is TRUE.
* `reduction` OPTIONAL the name of the PCA reduction. Not required if you used the default reduction names in PCA_dimRed. Default is NULL
* `data_type` OPTIONAL If the data is "RNA" or "spatial". Default is "RNA"
* `ffpe` OPTIONAL If the data type is "spatial" is it frozen (FALSE) or ffpe (TRUE), This is important because ffpe samples are from a probe based approach and don't include mitochondiral genes.

#### UMAP and clustering

One function exists to perform UMAP dimensonality reduction and clustering. This is because it is best to keep the number of PCs consistent when running both of these functions.

To generate a UMAP and clusters

```R
umap_res <- group_cells(so, sample_name = NULL, save_dir = NULL,
                        nPCs = 10, resolution = 0.8, assay = "RNA",
                        HTO = FALSE, reduction = NULL)

so <- umap_res$object 
umap_plots <- umap_res$plots
```

This function both generates a UMAP and clustering and generates helpful plots to visualize the clustering and the samples. Because it does both, this function returns a list that includes the plots and your object. Be sure to follow the code above to make sure that you are able to get both results.

In this function, setting a save_dir will also save the output plots as pdfs. Leaving as NULL will not save any plots.

### Differential expression

### Naming cell types

### Gene ontology

### Plotting

#### UMAP

* ggrastr

#### Violin/jitter plot

* ggrastr

#### Heatmaps


## Full example analysis pipeline
Example scripts using this package are available at https://github.com/CUAnschutzBDC/snakemake_pipelines/tree/main/scRNA_seq/src/scripts/indiviual_analysis. Follow the scripts in order for a complete analysis.

For examples of integrating data, scripts are available https://github.com/CUAnschutzBDC/snakemake_pipelines/tree/main/scRNA_seq/src/scripts/integrated_analysis
