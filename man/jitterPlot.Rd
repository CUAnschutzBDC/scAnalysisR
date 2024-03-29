% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{jitterPlot}
\alias{jitterPlot}
\title{Plot a jitter plot}
\usage{
jitterPlot(
  seurat_object,
  y_val,
  x_val,
  col_by = NULL,
  color = NULL,
  assay = NULL,
  size = 1,
  ggrastr = FALSE,
  raster_scale = 1,
  raster_res = 300
)
}
\arguments{
\item{seurat_object}{A seurat object}

\item{y_val}{What to use to create the y-axis}

\item{x_val}{What to use to separate along the x-axis}

\item{col_by}{OPTIONAL what to use to color the cells (or violins).}

\item{color}{OPTIONAL The color palette used to color. Default is Set1 from RColorBrewer}

\item{size}{OPTIONAL The size for the points. Default is 1}

\item{ggrastr}{OPTINAL If the plot should be rastarized. This is mostly helpful
for large datasets. All of the points will be a png while the rest will still
be editable. Default is FALSE (don't rasterize the plot)}

\item{raster_scale}{OPTIONAL The scale to use. Can be helpful if the rasterized
plot looks fuzzy. Default is 1}

\item{raster_res}{OPTIONAL. Can be helpful to change if the rasterized plot
looks fuzzy. Default is 300.}
}
\description{
This function is meant to be called by featDistPlot and not used on it's own!!!
This function will plot a jitter plot.
}
\keyword{internal}
