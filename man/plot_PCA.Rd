% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preprocessing.R
\name{plot_PCA}
\alias{plot_PCA}
\title{Make quality plots based on PCA}
\usage{
plot_PCA(
  sample_object,
  HTO = FALSE,
  ADT = FALSE,
  assay = "RNA",
  jackstraw = TRUE,
  reduction = NULL,
  data_type = "RNA",
  ffpe = TRUE
)
}
\arguments{
\item{sample_object}{A seurat object}

\item{HTO}{OPTIONAL if HTOs were included in the seurat object. Default is FALSE}

\item{ADT}{OPTIONAL If ADTs were included in the experiment. Default is FALSE}

\item{assay}{OPTIONAL What assay to use. This is just used to locate the PCA reduction
(the reduction name provided assumes that no reduction name was supplied to PCA_dimRed,
if you did provide a different reduction name, you must use "reduction" below). Default
is "RNA".}

\item{jackstraw}{OPTIONAL if a jackstraw plot should be made. This can help in
determining number of PCs to use, but it can be very slow. If TRUE, a jackstraw
plot will only be made if the assay is also "RNA". Default is TRUE.}

\item{reduction}{OPTIONAL the name of the PCA reduction. Not required if you used
the default reduction names in PCA_dimRed. Default is NULL}

\item{data_type}{OPTIONAL If the data is "RNA" or "spatial". Default is "RNA"}

\item{ffpe}{OPTIONAL If the data type is "spatial" is it frozen (FALSE) or ffpe (TRUE).
This is important because ffpe samples are from a probe based approach and don't
include mitochondiral genes.}
}
\value{
A named list of plots that differ based on the above parameters. Possibilities
include: orig.ident, percent.mt, nFeature_RNA, nCount_RNA, HTO_classification,
nFeature_ADT, nCount_ADT, nFeature_spatial, nCount_spatial, pca_loadings, and jackstraw
}
\description{
This function will create quality plots based on the PCA, including PC loadings,
the PCA colored by quality metrics, a knee plot, and a jackstraw plot.
}
\examples{
\dontrun{
plot_list <- plot_PCA(sample_object = seurat_object)
plot_list <- plot_PCA(sample_object = seurat_object,
                      HTO           = TRUE,
                      ADT           = TRUE)
plot_list <- plot_PCA(sample_object = seurat_object,
                      dta_type      = "spatial")
}
}
