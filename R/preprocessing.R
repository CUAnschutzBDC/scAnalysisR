#' Creates a seurat object
#' 
#' This function will create a seurat object. It assumes that the 10x output files are
#' in a directory named by the sample name followe by outs/count/filtered_feature_bc_matrix
#' This function can also generate matricies for HTO and ADTs. This function will also
#' perform CLR normlization on the ADT and HTO matricies, but will not normalize the
#' RNA data.
#' @param sample The name of the sample (also the name of the directory containing the 
#' count matricies).
#' @param count_path The full path to the 10x output files. In total this will be
#' count_path/sample/outs/count/filtered_feature_bc_matrix
#' @param ADT OPTIONAL If ADTs were included in the experiment. Default is TRUE
#' @param hashtag OPTIONAL If HTOs were included in the experiment. Default is TRUE
#' @param min_features OPTIONAL How many features per cell are necessary to keep the cell
#' will be passed to CreateSeuratObject. Default is 200
#' @param min_cells OPTIONAL the number of cells required to keep a feature, will be 
#' passed to CreateseuratObject. Default is 3.
#' @param hashtag_idents OPTIONAL the identities of the hashtags. If provided must be
#' exactly the identity that would be in the rowname. Default is to search for the 
#' word "Hashtag" in the identity.
#' @return A seurat object with assays for HTO and ADT (if HTO and ADT are true). 
#' If HTO and ADTs are included, those matricies will be normalized by CLR normalization
#' @import Seurat
#' @export
#' @examples
#' \dontrun{
#' seurat_object <- create_seurat_object(sample     = "WT_1",
#'                                       count_path = "/home/users/kwells/experiment_1/")
#' seurat_object <- create_seurat_object(sample     = "WT_1",
#'                                       count_path = "/home/users/kwells/experiment_2/",
#'                                       ADT        = FALSE,
#'                                       hashtag    = FALSE)
#' }

create_seurat_object <- function(sample, count_path, ADT = TRUE, hashtag = TRUE,
                                 min_features = 200, min_cells = 3,
                                 hashtag_idents = NULL){
  sample_path <- file.path(count_path, sample,
                           "outs", "count", "filtered_feature_bc_matrix")
  sample_data <- Read10X(data.dir = sample_path)
  if (ADT | hashtag){
    sample_object <- CreateSeuratObject(counts = sample_data[["Gene Expression"]],
                                        project = sample, min.cells = min_cells, 
                                        min.features = min_features)
    if (hashtag){
      protein_data <- sample_data[["Antibody Capture"]]
      if (is.null(hashtag_idents)){
        hashtag_data <- protein_data[grepl("Hashtag", rownames(protein_data)), ]
        ADT_data <- protein_data[!grepl("Hashtag", rownames(protein_data)), ]
      } else {
        hashtag_data <- protein_data[rownames(protein_data) %in% hashtag_idents, ]
        ADT_data <- protein_data[!(rownames(protein_data) %in% hashtag_idents), ]
      }
      sample_object[["HTO"]] <- CreateAssayObject(
        counts = hashtag_data[ ,Cells(sample_object)], min.cells = 3)
      sample_object <- NormalizeData(sample_object, assay = "HTO",
                                     normalization.method = "CLR")
    } else {
      ADT_data <- sample_data[["Antibody Capture"]]
    }
    if (ADT){
      sample_object[["ADT"]] <- CreateAssayObject(
        counts = ADT_data[ , Cells(sample_object)], min.cells = 3)
      sample_object <- NormalizeData(sample_object, assay = "ADT",
                                     normalization.method = "CLR",
                                     margin = 2)
    }
  } else {
    sample_object <- CreateSeuratObject(counts = sample_data,
                                        project = sample, min.cells = min_cells,
                                        min.features = min_features)
  }
  return(sample_object)
}

#' Performs PCA dimensionality reduction
#' 
#' This function will perform PCA dimensionality reduction and create a reduction
#' name for RNA, SCT, ADT, or integrated assays. There is currently no functionality
#' included for other assays.
#' @param sample_object A seurat object that has already been normalized
#' @param assay OPTIONAL The name of the assay to compute PCA for. Can be "RNA",
#' "ADT", "SCT", or "integrated". Default is "RNA"
#' @param reduction_name OPTIONAL The name to give the output PCA reduction. Defaults
#' are "pca" (for RNA assay), "apca" (for ADT assay), "sctpca" (for SCT assay), and
#' "integrated" (for integrated assay)
#' @param vars_to_regress OPTIONAL What variables to regress out when scaling the data.
#' Only used if assay = "ADT" (otherwise the data should already be scaled). Deefault
#' is NULL
#' @return A seurat object with an added reduction named based on the reduction_name
#' given or based on the rules described above.
#' @export
#' @examples
#' \dontrun{
#' seurat_object <- PCA_dimRed(sample_object = seurat_object)
#' seurat_object <- PCA_dimRed(sample_object = seurat_object,
#'                             assay         = "SCT")
#' }

PCA_dimRed <- function(sample_object, assay = "RNA",
                       reduction_name = NULL, vars_to_regress = NULL){
  if(assay == "RNA"){
    if(is.null(reduction_name)){
      reduction_name = "pca"
    }
    DefaultAssay(sample_object) = "RNA"
    sample_object <- RunPCA(sample_object,
                            features = VariableFeatures(object = sample_object),
                            reduction.name = reduction_name)
  } else if(assay == "ADT"){
    if(is.null(reduction_name)){
      reduction_name = "apca"
    }
    DefaultAssay(sample_object) <- "ADT"
    # Use all ADTS for dimensional reduction
    VariableFeatures(sample_object) <- rownames(sample_object[["ADT"]])
    sample_object <- Seurat::ScaleData(sample_object,
                                       vars.to.regress = vars_to_regress) %>% 
      Seurat::RunPCA(reduction.name = "apca")
  } else if(assay == "SCT"){
    if(is.null(reduction_name)){
      reduction_name = "sctpca"
    }
    DefaultAssay(sample_object) = "SCT"
    sample_object <- Seurat::RunPCA(sample_object,
                                    features =
                                      VariableFeatures(
                                        object = sample_object),
                                    reduction.name = "sctpca")
  } else if(assay == "integrated"){
    DefaultAssay(sample_object) = "integrated"
    sample_object <- RunPCA(sample_object,
                            features = VariableFeatures(object = sample_object))
  }
  return(sample_object)
}

#' Make quality plots based on PCA
#' 
#' This function will create quality plots based on the PCA, including PC loadings,
#' the PCA colored by quality metrics, a knee plot, and a jackstraw plot.
#' @param sample_object A seurat object
#' @param HTO OPTIONAL if HTOs were included in the seurat object. Default is FALSE
#' @param ADT OPTIONAL If ADTs were included in the experiment. Default is FALSE
#' @param assay OPTIONAL What assay to use. This is just used to locate the PCA reduction
#' (the reduction name provided assumes that no reduction name was supplied to PCA_dimRed,
#' if you did provide a different reduction name, you must use "reduction" below). Default
#' is "RNA".
#' @param jackstraw OPTIONAL if a jackstraw plot should be made. This can help in
#' determining number of PCs to use, but it can be very slow. If TRUE, a jackstraw
#' plot will only be made if the assay is also "RNA". Default is TRUE.
#' @param reduction OPTIONAL the name of the PCA reduction. Not required if you used
#' the default reduction names in PCA_dimRed. Default is NULL
#' @param data_type OPTIONAL If the data is "RNA" or "spatial". Default is "RNA"
#' @param ffpe OPTIONAL If the data type is "spatial" is it frozen (FALSE) or ffpe (TRUE).
#' This is important because ffpe samples are from a probe based approach and don't
#' include mitochondiral genes.
#' @return A named list of plots that differ based on the above parameters. Possibilities
#' include: orig.ident, percent.mt, nFeature_RNA, nCount_RNA, HTO_classification,
#' nFeature_ADT, nCount_ADT, nFeature_spatial, nCount_spatial, pca_loadings, and jackstraw
#' @export
#' @examples
#' \dontrun{
#' plot_list <- plot_PCA(sample_object = seurat_object)
#' plot_list <- plot_PCA(sample_object = seurat_object,
#'                       HTO           = TRUE,
#'                       ADT           = TRUE)
#' plot_list <- plot_PCA(sample_object = seurat_object,
#'                       dta_type      = "spatial")
#' }

plot_PCA <- function(sample_object, HTO = FALSE, ADT = FALSE,
                     assay = "RNA", jackstraw = TRUE,
                     reduction = NULL, data_type = "RNA",
                     ffpe = TRUE){
  if(!is.null(reduction)){
    reduction <- reduction
  } else if(assay == "RNA"){
    reduction <- "pca"
  } else if(assay == "ADT"){
    reduction <- "apca"
  } else if(assay == "SCT"){
    reduction <- "sctpca"
  } else if(assay == "integrated"){
    reduction <- "pca"
  } else {
    stop("'assay' must be 'ADT', 'SCT', 'integrated', or 'RNA' or
         a reduction must be supplied")
  }
  if(data_type == "RNA"){
    if(HTO){
      plot_list <- c("orig.ident", "percent.mt",
                     "nFeature_RNA", "nCount_RNA",
                     "HTO_classification")
    } else {
      plot_list <- c("orig.ident", "percent.mt",
                     "nFeature_RNA", "nCount_RNA")
    }
    if(ADT){
      plot_list <- c(plot_list, "nFeature_ADT", "nCount_ADT")
    }

  } else if(data_type == "spatial"){
    if(ffpe){
      plot_list <- c("orig.ident",
                     "nFeature_Spatial",
                     "nCount_Spatial")
    } else {
      plot_list <- c("orig.ident", "percent.mt",
                     "nFeature_Spatial",
                     "nCount_Spatial")   
    }
  }
  plots <- list()
  plots$pca_loadings <- VizDimLoadings(sample_object, dims = 1:2,
                                       reduction = reduction)
  feature_plots <- plotDimRed(sample_object, plot_type = reduction,
                              col_by = plot_list)
  
  names(feature_plots) <- plot_list
  
  plots <- c(plots, feature_plots)
  
  if(assay == "RNA" && jackstraw){
    sample_object <- JackStraw(sample_object, num.replicate = 100,
                               reduction = reduction)
    sample_object <- ScoreJackStraw(sample_object, dims = 1:20)
    plots$jackstraw <- JackStrawPlot(sample_object, dims = 1:20)
  }
  plots$elbow <- ElbowPlot(sample_object, reduction = reduction, ndims = 40)
  return(plots)
}

#' Run UMAP and clustering
#' 
#' This function will run both UMAP and clustering on a seurat object.
#' @param sample_object A seurat object
#' @param sample_name The name of the sample. This will be used for naming output
#' plots
#' @param save_dir Where to save plots that are generated (UMAPs of cluster and
#' original identity, and HTO identity if incldued)
#' @param nPCs OPTIONAL How many PCs to use to generate the UMAP and clustering. While
#' a default is provided, it is highly recommeded that you tune this to your own data.
#' Default is 10.
#' @param resolution OPTIONAL What resolution to use for clustering. Passed to FindClusters
#' A default is provided that tends to be a good starting point, but you are encouraged to
#' find the resolution appropriate to your dataset. Default is 0.8.
#' @param assay OPTIONAL What assay to use. This is used to locate the PCA reduction
#' (the reduction name provided assumes that no reduction name was supplied to PCA_dimRed,
#' if you did provide a different reduction name, you must use "reduction" below). 
#' This is also used to name the output plots, reductions, and clusters. Can be "RNA",
#' "ADT", "SCT", or "integrated" Default is "RNA".
#' @param HTO OPTIONAL if HTOs were included in the seurat object. Default is FALSE
#' @param reduction OPTIONAL the name of the PCA reduction. Not required if you used
#' the default reduction names in PCA_dimRed. Default is NULL
#' @param ... other arguments passed to plotDimRed
#' @return A seurat object with new cluster information in the metadata (RNA_cluster
#' if assay = "RNA", SCT_cluster if assay = "SCT", integrated_cluster if assay = 
#' "integrated", and ADT_cluster if assay = "ADT"), and a new UMAP reduction named
#' by the reduction you started with.
#' @export
#' @examples
#' \dontrun{
#' seurat_object <- group_cells(sample_object = seurat_object)
#' seurat_object <- group_cells(sample_object = seurat_object,
#'                              nPCs          = 25,
#'                              reduction     = 1.2)
#' }

group_cells <- function(sample_object, sample_name = NULL, save_dir = NULL,
                        nPCs = 10, resolution = 0.8, assay = "RNA",
                        HTO = FALSE, reduction = NULL, ...){
  if(assay == "RNA"){
    if(is.null(reduction)){
      reduction <- "pca"
      key <- "rna"
    } else {
      key <- reduction
    }
    DefaultAssay(sample_object) = "RNA"
    if(!is.null(save_dir)){
      save_plot <- file.path(save_dir, "images",
                             paste0("rnaUMAP_", sample_name, ".pdf"))
    }
    sample_object <- FindNeighbors(sample_object, dims = 1:nPCs,
                                   reduction = reduction)
    sample_object <- FindClusters(sample_object, resolution = resolution)
    sample_object <- RunUMAP(sample_object,
                             metric = "correlation", dims = 1:nPCs,
                             reduction = reduction, assay = assay,
                             reduction.key = paste0(key, "UMAP_"),
                             reduction.name = paste0(key, ".umap"))
    sample_object[["RNA_cluster"]] <- Idents(sample_object)
    col_by_list <- c("RNA_cluster", "orig.ident")
    if(HTO){
      col_by_list <- c(col_by_list, "HTO_classification")
    }
    if(is.null(save_dir)){
      save_plot <- NULL
    }
    plot_list <- plotDimRed(sample_object = sample_object,
                            save_plot = save_plot,
                            col_by = col_by_list, return_plot = TRUE,
                            plot_type = paste0(key, ".umap"), ...)
  } else if (assay == "ADT"){
    if(is.null(reduction)){
      reduction <- "apca"
      key <- "adt"
    } else {
      key <- reduction
    }
    DefaultAssay(sample_object) <- "ADT"
    save_plot <- file.path(save_dir, "images",
                           paste0("adtUMAP_", sample_name, ".pdf"))
    sample_object <- FindNeighbors(sample_object,
                                   features = rownames(sample_object),
                                   dims = 1:nPCs, reduction = reduction)
    sample_object <- FindClusters(sample_object, resolution = resolution,
                                  graph.name = "ADT_snn")
    sample_object <- RunUMAP(sample_object,
                             metric = "correlation", dims = 1:nPCs,
                             reduction = reduction, assay = assay,
                             reduction.key = paste0(key, "UMAP_"),
                             reduction.name = paste0(key, ".umap"))
    sample_object[["ADT_cluster"]] <- Idents(sample_object)
    col_by_list <- c("ADT_cluster", "orig.ident")
    if(HTO){
      col_by_list <- c(col_by_list, "HTO_classification")
    }
    if(is.null(save_dir)){
      save_plot <- NULL
    }
    plot_list <- plotDimRed(sample_object = sample_object,
                            save_plot = save_plot,
                            col_by = col_by_list, return_plot = TRUE,
                            plot_type = paste0(key, ".umap"), ...)
  } else if(assay == "SCT"){
    if(is.null(reduction)){
      reduction <- "sctpca"
      key <- "rnasct"
    } else {
      key <- reduction
    }
    DefaultAssay(sample_object) = "SCT"
    save_plot <- file.path(save_dir, "images",
                           paste0("rnasctUMAP_", sample_name, ".pdf"))
    sample_object <- FindNeighbors(sample_object, reduction = reduction,
                                   dims = 1:nPCs)
    sample_object <- FindClusters(sample_object, resolution = resolution)
    sample_object <- RunUMAP(sample_object,
                             metric = "correlation", dims = 1:nPCs,
                             reduction = reduction,  assay = assay,
                             reduction.key = paste0(key, "UMAP_"),
                             reduction.name = paste0(key, ".umap"))
    sample_object[["SCT_cluster"]] <- Idents(sample_object)
    col_by_list <- c("SCT_cluster", "orig.ident")
    if(HTO){
      col_by_list <- c(col_by_list, "HTO_classification")
    }
    if(is.null(save_dir)){
      save_plot <- NULL
    }
    plot_list <- plotDimRed(sample_object = sample_object,
                            save_plot = save_plot,
                            col_by = col_by_list, return_plot = TRUE,
                            plot_type = paste0(key, ".umap"), ...)
  } else if(assay == "integrated"){
    if(is.null(reduction)){
      reduction <- "pca"
    }
    DefaultAssay(sample_object) = "integrated"
    if(!is.null(save_dir)){
      save_plot <- paste0(save_dir, "images/integratedUMAP_", sample_name, ".pdf")
    }
    sample_object <- FindNeighbors(sample_object, dims = 1:nPCs,
                                   reduction = reduction)
    sample_object <- FindClusters(sample_object, resolution = resolution)
    sample_object <- RunUMAP(sample_object,
                             metric = "correlation", dims = 1:nPCs,
                             reduction = reduction, assay = assay,
                             reduction.key = "integratedUMAP_",
                             reduction.name = "integrated.umap")
    sample_object[["integrated_cluster"]] <- Idents(sample_object)
    col_by_list <- c("integrated_cluster", "orig.ident")
    if(HTO){
      col_by_list <- c(col_by_list, "HTO_classification")
    }
    if(is.null(save_dir)){
      save_plot <- NULL
    }
    plot_list <- plotDimRed(sample_object = sample_object,
                            save_plot = save_plot,
                            col_by = col_by_list, return_plot = TRUE,
                            plot_type = "integrated.umap", ...)
  }
  
  return(list(object = sample_object,
              plots = plot_list))
}

#' Creates a spatial seurat object
#' 
#' This function will create a seurat object. It needs the path to the directory
#' containing output from cellranger count that also includes the images.
#' @param results_dir The path to the results directory. The structure should be
#' results_dir/sample_name/outs
#' @param sample_name The name of the sample. This helps find the input data. The
#' path should be results_dir/sample_name/outs
#' @param filter_matrix OPTIONAL If the matrix should be filtered to only include
#' spots determined to be over tissue. Used by Read10X_Image.
#' @param filename OPTIONAL The name of the h5 file to load. Default is 
#' filtered_feature_bc_matrix.h5 which should generally work.
#' @return A seurat object with image data loaded.
#' @import Seurat
#' @export
#' @examples
#' \dontrun{
#' seurat_object <- create_spatial_seurat(sample     = "/home/users/kwells/experiment_1",
#'                                        count_path = "WT_1")
#' }

create_spatial_seurat <- function(results_dir, sample_name,
                                  filter_matrix = TRUE,
                                  filename = "filtered_feature_bc_matrix.h5"){
  results_path <- file.path(results_dir, sample_name, "outs")
  image_path <- file.path(results_path, "spatial")
  image_obj <- Read10X_Image(image.dir = image_path,
                             filter.matrix = filter_matrix)
  seurat_obj <- Load10X_Spatial(data.dir = results_path,
                                filename = filename,
                                filter.matrix = filter_matrix,
                                to.upper = FALSE,
                                image = image_obj)
  seurat_obj$orig.ident <- sample_name
  return(seurat_obj)
}