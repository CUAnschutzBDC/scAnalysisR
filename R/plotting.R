#' Plot a group of dimensionality reductions
#' 
#' This function will plot a dimensionality reduction that exists in your seurat
#' object (ex. UMAP or PCA). You can color the plot using any gene, adt (if they
#' exist) or meta data column. You can also highlight just one group from an
#' existing meta data column. This function will accept a vector in the col_by
#' argument. These can be any combination of genes or meta data columns. If
#' multiple options are given, a list of plots will be returned.
#' @param sample_object A seurat object
#' @param col_by What to color the plot by. May be a gene, adt, or column
#' from the metadata. It can either be discrete or continuous. This can be a vector
#' of selections (example c("cluster", "Ins1")). If a vector is given, the plots will
#' be returned as a list in the order of the supplied vector.
#' @param save_plot OPTIONAL A path to a file name if you want the plot saved. It
#' will be saved as a PDF. If nothing is provided, the plot will be returned but not
#' saved
#' @param plot_type OPTIONAL The type of plot. The possible plots can be found in
#' "reductions" in your seurat object. You can find these with names(object@reductions)
#' Default is "umap". 
#' @param dims_use OPTIONAL what dimensions to plot. Default is to plot dimensions 1 and 2
#' @param highlight_group OPTIONAL if only some cells should be highlighted on the plot.
#' when set to TRUE you must also set group and meta_data_col. When TRUE, the selected
#' group will be colored and all other cells will be grey.
#' @param meta_data_col The column in the seurat metadata (options can be found by calling
#' colnames(object[[]])) used for coloring the desired group. For example, if you want to
#' only highlight one cluster, meta_data_col = "cluster".
#' @param group The group to color. This group must be present in the selected meta data 
#' column.
#' @param reorder_cells OPTIONAL If cells should be randomly reodered to improve plotting 
#' aesthetics. Can be helpful if you don't want samples plotted in order. Default is FALSE.
#' @param ... OPTIONAL arguments supplied to plotDimRedSingle
#' @return a list of plots colored by the given parameteres
#' @import RColorBrewer
#' @import viridis
#' @export
#' @examples
#' \dontrun{
#' plotDimRed(sample_object = seurat_object,
#'            col_by        = "clusters",
#'            color         = c(1 = "#FFFFFF", 2 = "#FF0000"))
#' plotDimRed(sample_object = seurat_object,
#'            col_by        = "INS1")
#' plotDimRed(sample_object   = seurat_object,
#'            col_by          = "INS1",
#'            highlight_group = TRUE,
#'            group           = 1,
#'            meta_data_col   = "clusters")
#'}

plotDimRed <- function(sample_object, col_by, save_plot = NULL,
                       plot_type = "umap",
                       dims_use = NULL, highlight_group = FALSE,
                       group = NULL,
                       meta_data_col = "orig.ident", reorder_cells = FALSE, ...) {
  plot_list <- lapply(col_by, function(x) {
    plotDimRedSingle(seurat_object = sample_object, col_by = x, plot_type = plot_type,
                     dims_use = dims_use, highlight_group = highlight_group,
                     group = group, meta_data_col = meta_data_col,
                     reorder_cells = reorder_cells, ...)
  })
  if (!is.null(save_plot)){
    pdf(save_plot)
    print(plot_list)
    dev.off()
  }
  return(plot_list)
}

#' Plot one dimensionality reduction
#' 
#' This function is meant to be called by plotDimRed and not used on it's own!!!
#' This function will plot a dimensionality reduction that exists in your seurat
#' object (ex. UMAP or PCA). You can color the plot using any gene, adt (if they
#' exist) or meta data column. You can also highlight just one group from an
#' existing meta data column.
#' @param seurat_object A seurat object
#' @param col_by What to color the plot by. May be a gene, adt, or column
#' from the metadata. It can either be discrete or continuous.
#' @param plot_type OPTIONAL The type of plot. The possible plots can be found in
#' "reductions" in your seurat object. You can find these with names(object@reductions)
#' Default is "umap". 
#' @param dims_use OPTIONAL what dimensions to plot. Default is to plot dimensions 1 and 2
#' @param highlight_group OPTIONAL if only some cells should be highlighted on the plot.
#' when set to TRUE you must also set group and meta_data_col. When TRUE, the selected
#' group will be colored and all other cells will be grey.
#' @param meta_data_col The column in the seurat metadata (options can be found by calling
#' colnames(object[[]])) used for coloring the desired group. For example, if you want to
#' only highlight one cluster, meta_data_col = "cluster".
#' @param group OPTIONAL The group to color. This group must be present in the selected 
#' meta data column.
#' @param assay OPTIONAL The assay to use to color the plot. This is mostly useful if
#' you have the same names in your ADT and RNA slots so you can ensure you plot the
#' correct one.
#' @param reorder_cells OPTIONAL If cells should be randomly reodered to improve plotting 
#' aesthetics. Can be helpful if you don't want samples plotted in order. Default is FALSE.
#' @param ... OPTIONAL arguments supplied to groupDiscretePlots, groupContinuousPlots,
#' discretePlots, or continuousPlots
#' @keywords internal
#' @import RColorBrewer
#' @import viridis
#' @export

plotDimRedSingle <- function(seurat_object, col_by, plot_type = "umap",
                             dims_use = NULL, highlight_group = FALSE,
                             group = NULL, meta_data_col = "orig.ident",
                             assay = NULL, reorder_cells = FALSE, ...) {
  # Determine where in Seurat object to find variable to color by
  if (!is.null(assay) && col_by %in% rownames(seurat_object[[assay]])){
    DefaultAssay(seurat_object) <- assay
    col_by_data <- FetchData(object = seurat_object, vars = col_by)
  } else if (!is.null(assay)){
    stop(paste0("col_by (", col_by, ") is not present in your chosen assay"))
  } else if (col_by == "cluster" | col_by == "Cluster"){
    col_by_data <- as.data.frame(Idents(object = seurat_object))
  }else if (col_by %in% rownames(seurat_object) |
            col_by %in% colnames(seurat_object[[]])){
    col_by_data <- FetchData(object = seurat_object, vars = col_by)
  }else if ("ADT" %in% Seurat::Assays(seurat_object) &&
            col_by %in% rownames(seurat_object[["ADT"]])){
    col_by_data <- FetchData(object = seurat_object, vars = paste0("adt_", col_by))
  }else {
    stop(paste0("col_by (", col_by,
                ") must be a gene, metric from meta data or 'cluster'"))
  }
  
  # Make the name in the data frame the same regardless of what it was originally
  names(col_by_data) <- "colour_metric"
  
  if (is.null(dims_use)){
    dims_use <- c(1,2)
  }
  # Make a data frame based on the cell embeddings from the plot type of choice
  if (plot_type %in% names(seurat_object)){
    plot_coord <- Embeddings(object = seurat_object, reduction = plot_type)
    plot_names <- colnames(plot_coord)
    ndims <- length(plot_names)
    plot_cols <- lapply(dims_use, function(x){
      if (x > ndims) {
        stop("dims_use must be equal to or less than number of dimensions")
      } else {
        plot_col <- plot_names[x]
        return(plot_col)
      }
    })
    plot_cols <- unlist(plot_cols)
    plot_coord <- plot_coord[ , colnames(plot_coord) %in% plot_cols]
    axis_names <- colnames(plot_coord)
    colnames(plot_coord) <- c("dim1", "dim2")
    plot_df <- merge(plot_coord, col_by_data, by = "row.names")
    rownames(plot_df) <- plot_df$Row.names
  } else {
    stop("plot type must be a dimensional reduction in Seurat object")
  }
  # Add in group information if highlighting one group.
  if (highlight_group){
    if (is.null(group)){
      stop("if highlight_group is true, group must be a value from the meta_data
            column specified")
    }
    if (!identical(rownames(seurat_object[[]]), rownames(plot_df))) {
      plot_df <- plot_df[match(rownames(seurat_object[[]]),
                               rownames(plot_df)), , drop = FALSE]
    }
    plot_df[[meta_data_col]] <- seurat_object[[meta_data_col]][[1]]
    plot_df$all <- plot_df[[meta_data_col]]
    if (is.factor(plot_df$all)){
      plot_df$all <- factor(plot_df$all,
                            levels = c("all_samples", levels(plot_df$all)))
    }
    plot_df$all[!(plot_df$all %in% group)] <- "all_samples"
    
    # Randomly order the cells to change plotting order. Can be helpful if you don't
    # want samples plotted in order
    if(reorder_cells){
      plot_df <- plot_df[sample(nrow(plot_df)),]
    }

    # Plot as discrete
    if (!is.numeric(plot_df$colour_metric)){
      return_plot <- groupDiscretePlots(group, plot_df, axis_names = axis_names,
                                        col_by = col_by, ...)
      # Plot as continuous
    }else{
      return_plot <- groupContinuousPlots(group, plot_df, axis_names = axis_names,
                                          col_by = col_by, ...)
    }
  } else {
    # Reorder cells if necessary
    if(reorder_cells){
      plot_df <- plot_df[sample(nrow(plot_df)),]
    }
    # Plot as discrete
    if (!is.numeric(plot_df$colour_metric)){
      return_plot <- discretePlots(plot_df, axis_names = axis_names,
                                   col_by = col_by, ...)
      
      # Plot as continuous
    }else{
      return_plot <- continuousPlots(plot_df, axis_names = axis_names,
                                     col_by = col_by, ...)
    }
  }
  return(return_plot)
}

#' Plot one dimensionality reduction colored by discrete data
#' 
#' This function is meant to be called by plotDimRed and not used on it's own!!!
#' @param plot_df A dataframe created by plotDimRedSingle
#' @param col_by What to color the plot by. May be a gene, adt, or column
#' from the metadata. It can either be discrete or continuous.
#' @param axis_names OPTIONAL What will be supplied to label the x and y axis
#' @param color OPTIONAL what colors to use to color the plot. Default is Set1 from
#' RColorBrewer
#' @param show_legend OPTIONAL if a legend should be shown on the final plot. Deault is 
#' TRUE
#' @param size OPTIONAL The size of the points. Default is 0.25
#' @param ggrastr OPTINAL If the plot should be rastarized. This is mostly helpful
#' for large datasets. All of the points will be a png while the rest will still
#' be editable. Default is FALSE (don't rasterize the plot)
#' @param raster_scale OPTIONAL The scale to use. Can be helpful if the rasterized
#' plot looks fuzzy. Default is 1
#' @param raster_res OPTIONAL. Can be helpful to change if the rasterized plot
#' looks fuzzy. Default is 300.
#' @keywords internal
#' @import ggplot2
#' @import RColorBrewer
#' @export

discretePlots <- function(plot_df, col_by, axis_names = c("dim1", "dim2"),
                          color = NULL,show_legend = TRUE,
                          size = 0.25, ggrastr = FALSE,
                          raster_scale = 1, raster_res = 300){

  # Work around for "no visible binding"
  dim1 <- dim2 <- colour_metric <- NULL
  base_plot <- ggplot2::ggplot(data = plot_df,
                               ggplot2::aes(dim1, dim2)) +
    ggplot2::xlab(axis_names[1]) +
    ggplot2::ylab(axis_names[2])
  
  # Add colors based on metric chosen
  if (ggrastr){
    check_packages("ggrastr")
    base_plot <- base_plot +
      ggrastr::geom_point_rast(ggplot2::aes(colour = colour_metric),
                               show.legend = show_legend, size = size,
                               scale = raster_scale, raster.dpi = raster_res) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=3)))
  } else {
    base_plot <- base_plot +
      ggplot2::geom_point(ggplot2::aes(colour = colour_metric),
                          show.legend = show_legend, size = size) +
      ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=3)))
  }
  
  nColors <- length(levels(factor(plot_df$colour_metric)))
  
  # Color based on RColorBrewer if own palette isn't chosen
  if (is.null(color)) {
    base_plot <- base_plot + ggplot2::scale_color_manual(
      values = grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(9, "Set1"))(nColors), name = col_by)
  } else {
    base_plot <- base_plot + ggplot2::scale_color_manual(values = color, name = col_by)
  }
  
  return(base_plot)
}



#' Plot one dimensionality reduction colored by continuous data
#' 
#' This function is meant to be called by plotDimRed and not used on it's own!!!
#' @param plot_df A dataframe created by plotDimRedSingle
#' @param col_by What to color the plot by. May be a gene, adt, or column
#' from the metadata. It can either be discrete or continuous.
#' @param axis_names OPTIONAL What will be supplied to label the x and y axis
#' @param color OPTIONAL what colors to use to color the plot. Default is magma
#' from viridis
#' @param show_legend OPTIONAL if a legend should be shown on the final plot. Deault is 
#' TRUE
#' @param size OPTIONAL The size of the points. Default is 0.25
#' @param ggrastr OPTINAL If the plot should be rastarized. This is mostly helpful
#' for large datasets. All of the points will be a png while the rest will still
#' be editable. Default is FALSE (don't rasterize the plot)
#' @param raster_scale OPTIONAL The scale to use. Can be helpful if the rasterized
#' plot looks fuzzy. Default is 1
#' @param raster_res OPTIONAL. Can be helpful to change if the rasterized plot
#' looks fuzzy. Default is 300.
#' @keywords internal
#' @import ggplot2
#' @import viridis
#' @export

continuousPlots <- function(plot_df, col_by, axis_names = c("dim1", "dim2"),
                            color = NULL, show_legend = TRUE,
                            size = 0.25, ggrastr = FALSE,
                            raster_scale = 1, raster_res = 300){
  # Work around for "no visible binding"
  dim1 <- dim2 <- colour_metric <- labs <- NULL

  base_plot <- ggplot2::ggplot(data = plot_df, ggplot2::aes(dim1, dim2)) +
    ggplot2::xlab(axis_names[1]) +
    ggplot2::ylab(axis_names[2])
  
  if (ggrastr){
    check_packages("ggrastr")
    base_plot <- base_plot +
      ggrastr::geom_point_rast(ggplot2::aes(colour = colour_metric),
                               show.legend = show_legend, size = size,
                               scale = raster_scale, raster.dpi = raster_res) 
  } else {
    base_plot <- base_plot +
      ggplot2::geom_point(ggplot2::aes(colour = colour_metric),
                          show.legend = show_legend, size = size)
    
  }
  
  
  if (is.null(color)) {
    base_plot <- base_plot + ggplot2::scale_color_viridis_c(option = "magma") +
      labs(color = col_by)
  } else {
    low <- color[1]
    high <- color[2]
    base_plot <- base_plot + 
      ggplot2::scale_color_gradient(low = low, high = high, name = col_by)
  }
  
  return(base_plot)
}

#' Plot one dimensionality reduction colored by discrete data
#' 
#' This function is meant to be called by plotDimRed and not used on it's own!!!
#' @param plot_df A dataframe created by plotDimRedSingle
#' @param group The group to color.
#' @param col_by What to color the plot by. May be a gene, adt, or column
#' from the metadata. It can either be discrete or continuous.
#' @param axis_names OPTIONAL What will be supplied to label the x and y axis
#' @param color OPTIONAL what colors to use to color the plot. Default is Set1
#' from RColorBrewer
#' @param show_legend OPTIONAL if a legend should be shown on the final plot. Deault is 
#' TRUE
#' @param size OPTIONAL The size of the points. Default is 0.25
#' @param ggrastr OPTINAL If the plot should be rastarized. This is mostly helpful
#' for large datasets. All of the points will be a png while the rest will still
#' be editable. Default is FALSE (don't rasterize the plot)
#' @param raster_scale OPTIONAL The scale to use. Can be helpful if the rasterized
#' plot looks fuzzy. Default is 1
#' @param raster_res OPTIONAL. Can be helpful to change if the rasterized plot
#' looks fuzzy. Default is 300.
#' @keywords internal
#' @import ggplot2
#' @import RColorBrewer
#' @export

groupDiscretePlots <- function(group, plot_df, col_by, axis_names = c("dim1", "dim2"),
                               color = NULL, show_legend = TRUE,
                               size = 0.25, ggrastr = FALSE,
                               raster_scale = 1, raster_res = 300) {  
  # Work around for "no visible binding"
  dim1 <- dim2 <- colour_metric <- NULL
  plot1 <- plot_df[plot_df$all == "all_samples", ]
  plot2 <- plot_df[plot_df$all != "all_samples", ]

  base_plot <- ggplot2::ggplot(data = plot2, ggplot2::aes(dim1,
                                                          dim2)) +
    ggplot2::xlab(axis_names[1]) +
    ggplot2::ylab(axis_names[2])
  
  if(ggrastr){
    if (!requireNamespace("ggrastr", quietly = TRUE)){
      stop("Package \"ggrastr\" needed to make rasterized images. Please install it.",
           call. = FALSE)
    }
    
    base_plot <- base_plot + ggrastr::geom_point_rast(data = plot1, 
                                                      ggplot2::aes(dim1, dim2), 
                                                      color = "#DCDCDC",
                                                      size = size,
                                                      show.legend = FALSE,
                                                      scale = raster_scale,
                                                      raster.dpi = raster_res)
    base_plot <- base_plot + ggrastr::geom_point_rast(data = plot2,
                                                      ggplot2::aes(dim1, dim2,
                                                                   color = colour_metric),
                                                      size = size,
                                                      show.legend = show_legend,
                                                      scale = raster_scale,
                                                      raster.dpi = raster_res)  
  } else {
    base_plot <- base_plot + ggplot2::geom_point(data = plot1, 
                                                 ggplot2::aes(dim1, dim2), 
                                                 color = "#DCDCDC",
                                                 size = size,
                                                 show.legend = FALSE)
    base_plot <- base_plot + ggplot2::geom_point(data = plot2,
                                                 ggplot2::aes(dim1, dim2,
                                                              color = colour_metric),
                                                 size = size,
                                                 show.legend = show_legend)
  }
  base_plot <- base_plot + #ggplot2::theme_classic() + 
    ggplot2::ggtitle(paste(group, collapse = "_")) +
    ggplot2::xlab(axis_names[1]) +
    ggplot2::ylab(axis_names[2]) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size=3)))
  if (is.null(color)) {
    if(!is.null(levels(plot2$color_metric))){
      nColors <- length(levels(plot2$colour_metric))
    } else {
      nColors <- length(unique(plot2$colour_metric))
      
    }
    base_plot <- base_plot + ggplot2::scale_color_manual(
      values = grDevices::colorRampPalette(
        RColorBrewer::brewer.pal(9, "Set1"))(nColors), name = col_by)
  } else {
    base_plot <- base_plot +
      ggplot2::scale_color_manual(values = color, name = col_by)
  }
  
  return(base_plot)
}

#' Plot one dimensionality reduction colored by continuous data
#' 
#' This function is meant to be called by plotDimRed and not used on it's own!!!
#' @param plot_df A dataframe created by plotDimRedSingle
#' @param group The group to color.
#' @param col_by What to color the plot by. May be a gene, adt, or column
#' from the metadata. It can either be discrete or continuous.
#' @param axis_names OPTIONAL What will be supplied to label the x and y axis
#' @param color OPTIONAL what colors to use to color the plot. Default is magma from
#' viridis
#' @param show_legend OPTIONAL if a legend should be shown on the final plot. Deault is 
#' TRUE
#' @param size OPTIONAL The size of the points. Default is 0.25
#' @param ggrastr OPTINAL If the plot should be rastarized. This is mostly helpful
#' for large datasets. All of the points will be a png while the rest will still
#' be editable. Default is FALSE (don't rasterize the plot)
#' @param raster_scale OPTIONAL The scale to use. Can be helpful if the rasterized
#' plot looks fuzzy. Default is 1
#' @param raster_res OPTIONAL. Can be helpful to change if the rasterized plot
#' looks fuzzy. Default is 300.
#' @keywords internal
#' @import ggplot2
#' @import viridis
#' @export

groupContinuousPlots <- function(group, plot_df, col_by, color = NULL,
                                 limits = NULL, axis_names = c("dim1", "dim2"),
                                 save_plot = NULL, show_legend = TRUE,
                                 size = 0.25, ggrastr = FALSE,
                                 raster_scale = 1, raster_res = 300) {
  # Work around for "no visible binding"
  dim1 <- dim2 <- colour_metric <- labs <- NULL
  plot_name_comb <- paste(group, collapse = "_")
  plot1 <- plot_df[plot_df$all == "all_samples", ]
  plot2 <- plot_df[plot_df$all != "all_samples", ]
  
  base_plot <- ggplot2::ggplot(data = plot2, ggplot2::aes(dim1, dim2)) +
    ggplot2::xlab(axis_names[1]) +
    ggplot2::ylab(axis_names[2])
  
  if(ggrastr){
    if (!requireNamespace("ggrastr", quietly = TRUE)){
      stop("Package \"ggrastr\" needed to make rasterized images. Please install it.",
           call. = FALSE)
    }
    
    base_plot <- base_plot + ggrastr::geom_point_rast(data = plot1, 
                                                      ggplot2::aes(dim1, dim2), 
                                                      color = "#DCDCDC",
                                                      size = size,
                                                      show.legend = FALSE,
                                                      scale = raster_scale,
                                                      raster.dpi = raster_res)
    base_plot <- base_plot + ggrastr::geom_point_rast(data = plot2,
                                                      ggplot2::aes(dim1, dim2,
                                                                   color = colour_metric),
                                                      size = size,
                                                      show.legend = show_legend,
                                                      scale = raster_scale,
                                                      raster.dpi = raster_res)    
  } else {
    base_plot <- base_plot + ggplot2::geom_point(data = plot1, 
                                                 ggplot2::aes(dim1, dim2), 
                                                 color = "#DCDCDC",
                                                 size = size,
                                                 show.legend = FALSE)
    base_plot <- base_plot + ggplot2::geom_point(data = plot2,
                                                 ggplot2::aes(dim1, dim2,
                                                              color = colour_metric),
                                                 size = size,
                                                 show.legend = show_legend)
  }
  
  base_plot <- base_plot + #ggplot2::theme_classic() +
    ggplot2::ggtitle(paste0(plot_name_comb, " ", col_by)) +
    ggplot2::xlab(axis_names[1]) +
    ggplot2::ylab(axis_names[2])
  
  if (is.null(color)) {
    base_plot <- base_plot + ggplot2::scale_color_viridis_c(option = "magma") +
      labs(color = col_by)
  } else {
    low <- color[1]
    high <- color[2]
    if(is.null(limits)){
      base_plot <- base_plot + 
        ggplot2::scale_color_gradient(low = low, high = high, name = col_by)
    } else {
      base_plot <- base_plot + ggplot2::scale_color_gradient(low = low,
                                                             high = high, 
                                                             name = col_by,
                                                             limits = limits)
    }
  }
  
  return(base_plot)
}

#' Plot a group of violin or jitter plots
#' 
#' This function will plot a violin or jitter plot of the expression of a gene from
#' your seurat object. The plot could also be of some continuous data from your
#' meta data. The violin plot can be separated by any discrete variable in your
#' meta data. If multiple genes are supplied, a single plot conisiting of all violin 
#' or jitter plots will be returned with one plot per row.
#' @param seurat_object A seurat object
#' @param geneset What genes to use as the y axis in the plot. Can be many, but I
#' recommend not going above 3 because all plots will be returned as one object.
#' @param cell_cycle OPTIONAL Only an option if making a jitter plot. Can color each
#' cel by the cell cycle phase. I originally did this so the colors would be consistent.
#' Default is FALSE. If set to TRUE, it will override your col_by and color arguments.
#' @param plot_type OPTIONAL The type of plot. The options are "violin", "jitter", or
#' "both". "jitter" will make only a jitter plot, "violin" will make only a violin plot
#' and "both" will make a violin plot with jitter dots overlayed on tope.
#' @param col_by OPTIONAL what to use to color the cells (or violins). Generally for a violin
#' plot this is the same as sep_by, but for jitter, it can be any discrete column from meta
#' data. Default is NULL.
#' @param sep_by OPTIONAL what to use to separate the jitter or violins on the x axis. 
#' Can be any discrete value from the metadata. Deafult is "cluster"
#' @param save_plot OPTIONAL A path to a file name if you want the plot saved. It
#' will be saved as a PDF. If nothing is provided, the plot will be returned but not
#' saved
#' @param color The color palette used to color either the points for a jitter plot or
#' the violins. Default is Set1 from RColorBrewer
#' @param plot_median OPTIONAL if the median value should be included in the
#' violin plot. Default is TRUE
#' @param combine OPTIONAL if the final plots should be made into a figure or
#' returned as a list. Default (TRUE) returns a figure.
#' @param dodge OPTIONAL how to adjust the placement of the violin plot. Default is 1
#' @param width OPTIONAL how to adjust how wide the violin plots are. Default is 0.9
#' @param assay OPTIONAL what assay to pull the data from. Default is NULL.
#' @param ... arguments passed to other functions
#' @return Either a list of plots (if combine = FALSE) or a gridExtra object with all
#' plots stacked (if combine = TRUE)
#' @import ggplot2
#' @import RColorBrewer
#' @importFrom gridExtra arrangeGrob
#' @importFrom stats setNames
#' @export
#' @examples
#' \dontrun{
#' featDistPlot(seurat_object = seurat_object,
#'              geneset       = "INS1",
#'              col_by        = "cluster",
#'              sep_by        = "cluster",
#'              color         = c(1 = "#FFFFFF", 2 = "#FF0000"))
#' featDistPlot(seurat_object = seurat_object,
#'              geneset       = "INS1",
#'              col_by        = "cluster",
#'              sep_by        = "cluster",
#'              color         = c(1 = "#FFFFFF", 2 = "#FF0000"),
#'              plot_type     = "both")
#' featDistPlot(seurat_object = seurat_object,
#'              geneset       = c("INS1", "INS2", "GCG"),
#'              col_by        = "cluster",
#'              sep_by        = "cluster",
#'              color         = c(1 = "#FFFFFF", 2 = "#FF0000"),
#'              combine       = FALSE)
#'}

featDistPlot <- function(seurat_object, geneset, cell_cycle = FALSE,
                         plot_type = "violin",
                         color = NULL, sep_by = "cluster",
                         save_plot = NULL, col_by = NULL,
                         plot_median = TRUE, combine = TRUE,
                         dodge = 1, width = 0.9, assay = NULL, ...){
  # Work around for "no visible binding"
  x_value <- y_value <- NULL
  geneset <- setNames(geneset, geneset)
  if (plot_type == "jitter") {
    # Make jitter plots colored by cell cycle stage
    if(cell_cycle){
      gene_list_cycle <- lapply(geneset, function(x) jitterPlot(
        seurat_object = seurat_object, y_val = x, x_val = sep_by,
        col_by = "cycle_phase", color = c("black", "red", "purple"),
        assay = assay, ...))
      
      # Arrange all plots into one figure
      plot_list <- gridExtra::arrangeGrob(grobs = gene_list_cycle,
                                          nrow = length(geneset))
    } else {
      # Make a jitter plot based on expression of each gene given in the gene
      # set color by stage
      
      gene_list_stage <- lapply(geneset, function(x) jitterPlot(
        seurat_object = seurat_object, y_val = x, x_val = sep_by,
        color = color, col_by = col_by, assay = assay, ...))
      
      # Make a plot consisting of all plots made above
      plot_list <- gridExtra::arrangeGrob(grobs = gene_list_stage,
                                          nrow = length(geneset))
    }
    
  }
  if (plot_type == "violin" || plot_type == "both") {
    if (plot_type == "both"){
      plot_jitter <- TRUE
    } else {
      plot_jitter <- FALSE
    }
    gene_list_stage <- lapply(geneset, function(x) violinPlot(
      seurat_object = seurat_object, y_val = x, x_val = sep_by,
      color = color, plot_jitter = plot_jitter, col_by = col_by,
      plot_median = plot_median, dodge = dodge, width = width,
      assay = assay, ...))
    
    plot_list <- gridExtra::arrangeGrob(grobs = gene_list_stage,
                                        nrow = length(geneset))
    
  }
  if (!(is.null(save_plot))){
    ggplot2::ggsave(save_plot, plot = plot_list)
  }
  if(combine){
    return(plot_list)
  } else {
    return(gene_list_stage)
  }
}

#' Plot a jitter plot
#' 
#' This function is meant to be called by featDistPlot and not used on it's own!!!
#' This function will plot a jitter plot.
#' @param seurat_object A seurat object
#' @param y_val What to use to create the y-axis
#' @param x_val What to use to separate along the x-axis
#' @param col_by OPTIONAL what to use to color the cells (or violins).
#' @param color OPTIONAL The color palette used to color. Default is Set1 from RColorBrewer
#' @param size OPTIONAL The size for the points. Default is 1
#' @param ggrastr OPTINAL If the plot should be rastarized. This is mostly helpful
#' for large datasets. All of the points will be a png while the rest will still
#' be editable. Default is FALSE (don't rasterize the plot)
#' @param raster_scale OPTIONAL The scale to use. Can be helpful if the rasterized
#' plot looks fuzzy. Default is 1
#' @param raster_res OPTIONAL. Can be helpful to change if the rasterized plot
#' looks fuzzy. Default is 300.
#' @keywords internal
#' @import ggplot2
#' @import RColorBrewer
#' @export

jitterPlot <- function(seurat_object, y_val, x_val,
                       col_by = NULL, color = NULL,
                       assay = NULL, size = 1,
                       ggrastr = FALSE,
                       raster_scale = 1, 
                       raster_res = 300) {
  # Work around for "no visible binding"
  x_value <- y_value <- NULL
  plot_data <- plotDF(seurat_object, y_val, x_val,
                      col_by, assay = assay)
  # Determine the number of different colors needed.
  nColors <- length(unique(plot_data$col_by))
  
  # Plot the main plot
  plot_base <- ggplot2::ggplot(data = plot_data, ggplot2::aes(x_value, 
                                                              y_value,
                                                              color = col_by)) +
    #ggplot2::theme_classic() + 
    ggplot2::ylab(y_val) + ggplot2::xlab(x_val)
  
  if (ggrastr){
    check_packages("ggrastr")
    plot_base <- plot_base +
      ggrastr::geom_jitter_rast(shape = 16, size = size,
                                scale = raster_scale, raster.dpi = raster_res)
  } else {
    plot_base <- plot_base +
      ggplot2::geom_jitter(shape = 16, size = size) 
    
  }
  
  if(is.null(col_by)){
    plot_base <- plot_base + 
      ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                     axis.text.x=ggplot2::element_blank())
    
  } else {
    plot_base <- plot_base +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }
  
  if (is.null(color)) {
    plot_base <- plot_base +
      ggplot2::scale_color_manual(values =
                                    (grDevices::colorRampPalette(
                                      RColorBrewer::brewer.pal(9, 
                                                               "Set1")))(nColors), name = col_by)
  } else {
    plot_base <- plot_base + ggplot2::scale_color_manual(values = color, 
                                                         name = col_by)  
  }
  
  return(plot_base)
}  

#' Plot a violin plot
#' 
#' This function is meant to be called by featDistPlot and not used on it's own!!!
#' This function will plot a violin plot.
#' @param seurat_object A seurat object
#' @param y_val What to use to create the y-axis
#' @param x_val What to use to separate along the x-axis
#' @param col_by OPTIONAL what to use to color the cells (or violins).
#' @param color OPTIONAL The color palette used to color. Default is Set1 from RColorBrewer
#' @param plot_jitter OPTIONAL if a jitter plot should also be made. Default is FALSE.
#' @param plot_median OPTIONAL if the median value should be included in the
#' violin plot. Default is TRUE
#' @param dodge OPTIONAL how to adjust the placement of the violin plot. Default is 1
#' @param width OPTIONAL how to adjust how wide the violin plots are. Default is 0.9
#' @param size OPTIONAL the size of the points for the jitter plot. Default is 1.
#' @param ggrastr OPTINAL If the plot should be rastarized. This is mostly helpful
#' for large datasets. All of the points will be a png while the rest will still
#' be editable. Default is FALSE (don't rasterize the plot)
#' @param raster_scale OPTIONAL The scale to use. Can be helpful if the rasterized
#' plot looks fuzzy. Default is 1
#' @param raster_res OPTIONAL. Can be helpful to change if the rasterized plot
#' looks fuzzy. Default is 300.
#' @keywords internal
#' @import ggplot2
#' @import RColorBrewer
#' @export

violinPlot <- function(seurat_object, y_val, x_val,
                       col_by = NULL, color = NULL,
                       plot_jitter = FALSE,
                       plot_median = TRUE,
                       dodge = 1, width = 0.9,
                       assay = NULL, size = 1,
                       ggrastr = FALSE,
                       raster_scale = 1, 
                       raster_res = 300) {
  # Work around for "no visible binding"
  x_value <- y_value <- median <- NULL 
  plot_data <- plotDF(seurat_object, y_val, x_val,
                      col_by, assay = assay)
  # Determine the number of different colors needed.
  nColors <- length(unique(plot_data$col_by))
  
  # Plot the main plot
  plot_base <- ggplot2::ggplot(data = plot_data, ggplot2::aes(x_value, 
                                                              y_value,
                                                              fill = col_by)) +
    #ggplot2::theme_classic() + 
    ggplot2::ylab(y_val) + ggplot2::xlab(x_val) +
    ggplot2::geom_violin(scale = "width",
                         position = ggplot2::position_dodge(dodge),
                         width = width)
  if(is.null(col_by)){
    plot_base <- plot_base +
      ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                     axis.text.x=ggplot2::element_blank())
  } else {
    plot_base <- plot_base +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }
  if (plot_jitter) {
    if (ggrastr){
      check_packages("ggrastr")
      plot_base <- plot_base +
        ggrastr::geom_jitter_rast(shape = 16, size = size, scale = raster_scale,
                                  raster.dpi = raster_res)
    } else {
      plot_base <- plot_base +
        ggplot2::geom_jitter(shape = 16, size = size) 
      
    }
  }
  
  if (is.null(color)) {
    plot_base <- plot_base +
      ggplot2::scale_fill_manual(values =
                                   (grDevices::colorRampPalette(
                                     RColorBrewer::brewer.pal(9, 
                                                              "Set1")))(nColors), name = col_by)
  } else {
    plot_base <- plot_base + ggplot2::scale_fill_manual(values = color, 
                                                        name = col_by)  
  }
  if(plot_median){
    plot_base <- plot_base +
      ggplot2::stat_summary(fun = median, geom = "point", size = 2,
                            position = ggplot2::position_dodge(dodge))
  }
  
  return(plot_base)
}

#' Make a dataframe
#' 
#' This function is meant to be called by featDistPlot and not used on it's own!!!
#' This function will make a dataframe to pass to plotting functions
#' @param seurat_object A seurat object
#' @param y_val What to use to create the y-axis
#' @param x_val What to use to separate along the x-axis
#' @param col_by OPTIONAL what to use to color the cells (or violins).
#' @keywords internal
#' @export

plotDF <- function(seurat_object, y_val, x_val,
                   col_by = NULL, assay = NULL) {
  # Add y_value to a data frame used for plotting. This value can be a gene
  # or a value from meta data like nGene
  # Determine where in Seurat object to find variable to color by
  if (!is.null(assay) && y_val %in% rownames(seurat_object[[assay]])){
    DefaultAssay(seurat_object) <- assay
    plot_data <- FetchData(object = seurat_object, vars = y_val)
  } else if (!is.null(assay)){
    stop(paste0("gene (", y_val, ") is not present in your chosen assay"))
  } else if (y_val == "cluster" | y_val == "Cluster"){
    plot_data <- as.data.frame(Idents(object = seurat_object))
  }else if (y_val %in% rownames(seurat_object) |
            y_val %in% colnames(seurat_object[[]])){
    plot_data <- FetchData(object = seurat_object, vars = y_val)
  }else if ("ADT" %in% Seurat::Assays(seurat_object) &&
            y_val %in% rownames(seurat_object[["ADT"]])){
    plot_data <- FetchData(object = seurat_object, vars = paste0("adt_", y_val))
  }else {
    stop(paste0("gene (", y_val,
                ") must be a gene, metric from meta data or 'cluster'"))
  }
  # Name the column
  names(plot_data) <- "y_value"
  
  # Add a column contining the x_value. This should be something discrete
  # Like timepoint or cluster
  if (x_val %in% colnames(seurat_object[[]])) {
    # Should be able to fix this. Take out as and do x_val = then don't need names()
    x_plot_data <- as.data.frame(seurat_object[[]][, x_val, drop = FALSE])
    #x_plot_data <- data.frame("x_value" = seurat_object@meta.data[, x_val,
    #                                                             drop = FALSE])
    plot_data <- merge(plot_data, x_plot_data, by = "row.names")
    rownames(plot_data) <- plot_data$Row.names
    plot_data$Row.names <- NULL
  } else if (x_val == "cluster") {
    x_plot_data <- as.data.frame(Idents(seurat_object))
    #x_plot_data <- data.frame("x_value" = seurat_object@ident)
    plot_data <- merge(plot_data, x_plot_data, by = "row.names")
    rownames(plot_data) <- plot_data$Row.names
    plot_data$Row.names <- NULL
  } else {
    stop("x_val must be a metric from meta data or 'cluster'")
  }
  
  # Name the appropriate column of the plotting data
  names(plot_data)[2] <- "x_value"
  plot_data <- plot_data[match(colnames(seurat_object),
                               rownames(plot_data)), ]
  
  # Determine how to color the plot. Default is the x_value but can be any
  # discrete value.
  if (is.null(col_by)) {
    plot_data$col_by <- plot_data$x_value
  } else if (col_by %in% colnames(seurat_object[[]])) {
    plot_data$col_by <- seurat_object[[]][ , col_by]
  } else if (col_by == "cluster") {
    plot_data$col_by <- Idents(seurat_object)
  } else {
    stop("x_val must be a metric from meta data or 'cluster'")
  }
  return(plot_data)
}

#' Plot a heatmap
#' 
#' This function will plot a heatmap with annotations based on the metadata from a 
#' supplied seurat object. The heatmap is colored by the blueYellow palette from
#' the archR package.
#' @param seurat_object A seurat object
#' @param gene_list What genes to use in the heatmap. The function will take a long time
#' and the heatmap will be huge if you use all genes.
#' @param meta_col What column from the meta data should be used to annotate the x axis
#' of the heatmap.
#' @param colors OPTIONAL What colors should be used to annotate the x axis. Default is
#' Set1 from RColorBrewer
#' @param meta_df OPTIONAL If you want more complex annotations than just annotating by
#' one column of the meta data, you can include your own meta_df here. If you use this
#' option, you must include your own color list, where the names of the list match
#' each column of your meta_df and the vector of colors in each item of the list match
#' the number of items in the meta data column. Additionally, if you use this
#' option, you must make sure that the cell names from the meta data match the cell names
#' in your Seurat object. I've found that this is a particular issue when using average
#' expression of clusters.
#' @param color_list OPTIONAL only needs to be included when meta_df is set. This should be
#' a list where the names of the list match each column of your meta_df and the vector of
#' colors in each item of the list match the number of items in the meta data column.
#' @param max_val OPTIONAL What value to cap the z-score shown. This makes it so highly
#' variable genes don't completely controll the color scale. Default is 2.5
#' @param min_val OPTIONAL What value to cap the z-score shown. This makes it so highly
#' variable genes don't completely controll the color scale. Default is -2.5
#' @param cluster_rows OPTIONAL if rows should be clustered. Default is FALSE
#' @param cluster_cols OPTIONAL if columns should be clustered. Default is FALSE
#' @param average_expression OPTIONAL if the average expression of clusters should be 
#' computed before making the heatmap. This option makes much prettier plots, but
#' can be a bit buggy, I'm still working on making it work more robustly.
#' @param plot_meta_col OPTIONAL If the `meta_col` should be included in the coloring on the top
#' of the plot. I generally turn this off if I am plotting many features on the top and
#' have a custum meta_df. Default is TRUE
#' @param plot_rownames OPTIONAL if rownames should be printed. Default is TRUE. I generally
#' like to print the rownames unless there are too many to be visually pleasing.
#' @param cell_order OPTIONAL a custom order of the cells to plot. Default is NULL
#' and cells will be plotted in the order of the object or by clustering.
#' @param return_data OPTIONAL If the scaled matrix and count matrix should be returned.
#' if TRUE, a list will be returned with the heatmap and both matricies, if FALSE only
#' the heatmap will be returned. Default is FALSE
#' @param assay OPTIONAL what assay should be used for generating the plot
#' @param ... Other options passed to `pheatmap`
#' @return A pheatmap object with x axis colors based on the meta data column provided
#' or your own meta_df with color_list.
#' @import RColorBrewer
#' @import grDevices
#' @export
#' @examples
#' \dontrun{
#' plot_heatmap(seurat_object = seurat_object,
#'              gene_list     = c("INS1", "INS2", "GCG", "SST", "NKX2-2"),
#'              meta_col      = "clusters",
#'              color         = c(1 = "#FFFFFF", 2 = "#FF0000"))
#' plot_heatmap(seurat_object      = seurat_object,
#'              gene_list          = c("INS1", "INS2", "GCG", "SST", "NKX2-2"),
#'              meta_col           = "clusters",
#'              color              = c(1 = "#FFFFFF", 2 = "#FF0000"),
#'              average_expression = TRUE)
#' #' plot_heatmap(seurat_object = seurat_object,
#'              gene_list     = gene_set,
#'              meta_col      = "clusters")
#' }

plot_heatmap <- function(seurat_object, gene_list, meta_col,
                         colors = NULL, meta_df = NULL, color_list = NULL,
                         max_val = 2.5, min_val = -2.5, cluster_rows = FALSE,
                         cluster_cols = FALSE, average_expression = FALSE,
                         plot_meta_col = TRUE, plot_rownames = TRUE,
                         cell_order = NULL, return_data = FALSE,
                         assay = "RNA", ...){
  if(average_expression){
    # Find average expression of genes in clusters
    Idents(seurat_object) <- meta_col
    heatmap_df <- get_average_expression(seurat_object, seurat = FALSE,
                                         group.by = meta_col, assay = assay)
    # Test if the colnames look like integers
    character_vals <- 
      suppressWarnings(all(!is.na(as.numeric(as.character(colnames(heatmap_df))))))
    if(is.null(meta_df)){
      sample_info <- seurat_object[[meta_col]]
      # Add levels
      if(is.null(levels(sample_info[[meta_col]]))){
        sample_info[[meta_col]] <- factor(sample_info[[meta_col]])
      }
      meta_df <- data.frame(levels(sample_info[[meta_col]]))
      colnames(meta_df) <- meta_col
      if(character_vals){
        rownames(meta_df) <- paste0("X", meta_df[[meta_col]])
      } else {
        rownames(meta_df) <- meta_df[[meta_col]]
      }
      if(is.null(colors)){
        colors <- brewer.pal(length(levels(sample_info[[meta_col]])), "Set1")
        names(colors) <- levels(sample_info[[meta_col]])
      } 
      # make a list for the column labeing
      color_list <- list(colors)
      names(color_list) <- meta_col
    }
  } else {
    # Pull out data and subset to genes of interest
    heatmap_df <- get_seurat_assay(seurat_object, type = "data", assay = assay)
  }
  heatmap_df <- heatmap_df[rownames(heatmap_df) %in% gene_list, ]
  heatmap_df <- data.frame(heatmap_df)
  
  heatmap_df <- heatmap_df[order(match(rownames(heatmap_df), gene_list)), ]
  
  # remove any zero values
  heatmap_df <- heatmap_df[rowSums(heatmap_df) > 0,]
  
  if(is.null(meta_df)){
    # Make a df for the column labeling
    sample_info <- seurat_object[[meta_col]]
    # Add levels
    if(is.null(levels(sample_info[[meta_col]]))){
      sample_info[[meta_col]] <- factor(sample_info[[meta_col]])
    }
    if(is.null(colors)){
      colorcount <- length(levels(sample_info[[meta_col]]))
      colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))(colorcount)
      names(colors) <- levels(sample_info[[meta_col]])
    } 
    # make a list for the column labeing
    coloring <- list(colors)
    names(coloring) <- meta_col
  } else {
    # The sample info and color list must be provided
    sample_info <- meta_df
    coloring <- color_list
  }
  
  # Set cluster order
  
  cluster_order <- levels(sample_info[[meta_col]])
  heatmap_scale <- t(scale(t(heatmap_df), scale = TRUE))
  # Colors for heatmap (from the ArchR package)
  blueYellow <- c("#352A86", "#343DAE", "#0262E0", "#1389D2", "#2DB7A3",
                  "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D")
  
  if(!is.null(cell_order)){
    sample_info <- sample_info[order(match(rownames(sample_info),
                                           cell_order)), , drop = FALSE]
    rownames(sample_info) <- make.names(rownames(sample_info))
    if (!identical(colnames(heatmap_scale), rownames(sample_info))) {
      heatmap_scale <- heatmap_scale[, rownames(sample_info)]
    }
  } else if (!cluster_cols) {
    sample_info <- sample_info[order(match(sample_info[[meta_col]], 
                                           cluster_order)), , drop = FALSE]
    rownames(sample_info) <- make.names(rownames(sample_info))
    if (!identical(colnames(heatmap_scale), rownames(sample_info))) {
      heatmap_scale <- heatmap_scale[, rownames(sample_info)]
    }
  } else {
    rownames(sample_info) <- make.names(rownames(sample_info))
  }
  
  if(!plot_meta_col){
    sample_info[[meta_col]] <- NULL
  }
  
  # This makes the values more even
  heatmap_scale <- ifelse(heatmap_scale > max_val, max_val, heatmap_scale)
  heatmap_scale <- ifelse(heatmap_scale < min_val, min_val, heatmap_scale)
  
  heatmap <- pheatmap::pheatmap(heatmap_scale, cluster_rows = cluster_rows,
                                cluster_cols = cluster_cols,
                                show_rownames = plot_rownames,
                                show_colnames = FALSE, annotation_col = sample_info,
                                annotation_colors = coloring, color = blueYellow,
                                border_color = NA, clustering_method = "complete",
                                silent = TRUE, ...)
  if(return_data){
    return(list("heatmap" = heatmap,
                "z_score" = heatmap_scale,
                "counts" = heatmap_df))
  } else {
    return(heatmap)
  }
}

#' Create stacked barplot
#' 
#' This function will create stacked barplots based on cells in each cluster. Barplots
#' can be separated by any meta data column, for example, sample or treatment or 
#' hashtag
#' @param seurat_object A seurat object
#' @param meta_col What column to use to generate the barplot. For example, "cluster"
#' or "cell_type"
#' @param color OPTIONAL What colors should be used to color the stacked barplot.
#' Default is Set1 from RColorBrewer.
#' @param percent If the barplots should be constructed using the percent of cells in 
#' each category (TRUE) or the count (FALSE). Default is TRUE.
#' @param split_by How the barplots should be split. Can be any discrete column from the
#' meta data. For example, "sample". If nothing is provided, only one barplot will be made.
#' @param return_values If the percents and counts should be returned. If TRUE, the return
#' from this object will be a data frame with counts and percents by the specified meta data
#' columns for meta_col and split_by. If FALSE the return will be a ggplot object of the 
#' stacked barplots. Default is FALSE.
#' @return Either a stacked barplot based on the count or percent of each meta data group
#' in your sample (if return_values = FALSE) or a data frame with the count and percent
#' values of each meta data group in your sample (if return_values = TRUE).
#' @import dplyr
#' @import ggplot2
#' @import RColorBrewer
#' @export
#' @examples
#' \dontrun{
#' stacked_barplots(seurat_object = seurat_object,
#'                  meta_col      = "cluster")
#' stacked_barplots(seurat_object = seurat_object,
#'                  meta_col      = "cluster",
#'                  color         = c(1 = "#FFFFFF", 2 = "#FF0000"),
#'                  split_by      = "sample")
#' percent_df <- stacked_barplots(seurat_object = seurat_object,
#'                                meta_col      = "cluster",
#'                                return_values = TRUE)
#'}

stacked_barplots <- function(seurat_object, meta_col, color = NULL,
                             percent = TRUE, split_by = NULL,
                             return_values = FALSE){
  # Work around for "no visible binding"
  Freq <- percents <- NULL
  if(!is.null(split_by)){
    meta_data <- Seurat::FetchData(seurat_object, vars = c(meta_col, split_by)) %>%
      dplyr::rename(meta_col = 1, split_by = 2) %>%
      dplyr::group_by(split_by) %>%
      base::table() %>%
      base::data.frame()
    
  } else {
    meta_data <- Seurat::FetchData(seurat_object, vars = meta_col) %>%
      base::table() %>%
      base::data.frame() %>%
      dplyr::rename(meta_col = 1) %>%
      dplyr::mutate(split_by = "group1") 
  }
  meta_data <- meta_data %>%
    dplyr::group_by(split_by) %>%
    dplyr::mutate(percents = Freq/sum(Freq) * 100)
  # Add colors if not provided
  if(is.null(color)){
    color <- brewer.pal(length(unique(meta_data$meta_col)), "Set1")
  }
  if(percent){
    bar_plot <- ggplot2::ggplot(meta_data, ggplot2::aes(x = split_by,
                                                        y = percents,
                                                        fill = meta_col)) +
      ggplot2::ylab("Percent")
  } else {
    bar_plot <- ggplot2::ggplot(meta_data, ggplot2::aes(x = split_by,
                                                        y = Freq,
                                                        fill = meta_col)) +
      ggplot2::ylab("Count")
  }
  bar_plot <- bar_plot +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_fill_manual(values = color) +
    ggplot2::xlab(split_by) +
    ggplot2::labs(fill = meta_col)
  
  if(return_values){
    return(list(data = meta_data,
                barplot = bar_plot))
  } else {
    return(bar_plot)
  }
}
