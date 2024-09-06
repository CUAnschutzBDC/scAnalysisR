#' Name clusters based on a reference
#' 
#' This function uses clustifyr to name clusters based on a reference dataset.
#' @param seurat_object A seurat object
#' @param ref_mat A matrix from a reference. This should be genes by cell type. If it
#' comes from a single cell dataset, you can make this reference using average_clusters
#' from clustifyr. See the clustifyr tutorial for more information.
#' @param save_dir A path to a directory to save the output matricies and plots.
#' @param save_name The name to use to save the cell type information. Example is 
#' t_cell_ref. The returned seurat object will have a column name RNA_t_cell_ref
#' in the meta data.
#' @param ADT OPTIONAL If you want to also name cell types based on clusters determined
#' based on ADTs. Only recommened if you used a large pannel of ADTs. To use this,
#' PCA_dimRed and group_cells must have been run on the ADT assay. Additionally, a combined
#' cluster must also have been generated. The steps to do this will eventually be added to
#' this package, but they have not yet.
#' @param assay OPTIONAL What assay to use to compare to the reference. "RNA" is
#' receommended. Default is "RNA"
#' @param nfeatures OPTIONAL The number of features to use to run clustifyr. Check the 
#' clustifyr tutorial for a recommendation. Normally 1000 is a good place to start. Too
#' high of a value will lead to many of the cell types having a similar correlation.
#' Default is 1000.
#' @param clusters OPTIONAL What existing clusters to use from your seurat object to identify 
#' cell types with clustifyr. Can be any discrete column from your metadata. Default is
#' "RNA_cluster"
#' @param plot_type OPTIONAL What plot type should be used to plot the results. These plots
#' will be saved in the provided directory in a subdirectory "images". Can be any plot
#' type in the seurat object. Default is "rna.umap".
#' @param cor_cutoff OPTIONAL The correlation value to use to determine if a cell type
#' is named. If no cell types have a correlation above this value for a given cluster,
#' that cluster will be named "undetermined" in the final output. Default is 0.5
#' @param features OPTIONAL A list of features to use instead of running FindVariableFeatures
#' Default is NULL. If set, will override the nfeatures argument. This is helpful if you
#' specifically want to exclude features (for example Ig genes or TCR genes).
#' @return a list containing the correlation matricies and the seurat object with a new
#' meta data column including the identified cell type.
#' @import pheatmap
#' @import dplyr
#' @importFrom viridis viridis
#' @export
#' @examples
#' \dontrun{
#' clustify_list <- name_clusters(seurat_object = seurat_object,
#'                                ref_mat       = t_cell_matrix,
#'                                save_dir      = "clustifyr",
#'                                save_name     = t_cell_ref)
#' seurat_object <- clustify_list$object
#'}

name_clusters <- function(seurat_object, ref_mat, save_dir,
                          save_name, ADT = FALSE,
                          assay = "RNA",
                          nfeatures = 1000, clusters = "RNA_cluster",
                          plot_type = "rna.umap",
                          cor_cutoff = 0.5,
                          features = NULL){
  # Work around for "no visible binding"
  r <- type <- NULL
  
  # Ask user to install clustifyr to use this function
  check_packages("clustifyr")
  requireNamespace("clustifyr")


  # Make directories if they don't exist
  ifelse(!dir.exists(file.path(save_dir, "images")),
         dir.create(file.path(save_dir, "images")), FALSE)

  ifelse(!dir.exists(file.path(save_dir, "files")),
         dir.create(file.path(save_dir, "files")), FALSE)

  # Keep only genes in seurat object
  ref_mat <- ref_mat[rownames(ref_mat) %in% rownames(seurat_object),]
  
  return_list <- list()
  
  ##########
  # Set up #
  ##########
  # get count matrix
  DefaultAssay(seurat_object) <- assay
  seurat_mat <- get_seurat_assay(object = seurat_object, type = "data", assay = assay)
  
  if(is.null(features)){
    # Find only the specified number of variable features
    seurat_var <- FindVariableFeatures(
      seurat_object,
      assay = assay,
      selection.method = "vst",
      nfeatures = nfeatures
    )
    seurat_genes <- VariableFeatures(seurat_var)    
    } else {
      seurat_genes <- features
    }

  
  # RNA
  
  seurat_metadata <- seurat_object[[clusters]]
  
  seurat_metadata[[clusters]] <- as.character(seurat_metadata[[clusters]])
  
  seurat_object[[clusters]] <- seurat_metadata[[clusters]]
  
  # Run clustify
  seurat_res <- clustifyr::clustify(
    input = seurat_mat,
    metadata = seurat_metadata,
    ref_mat = ref_mat,
    query_genes = seurat_genes,
    cluster_col = clusters
  )
  
  return_list$RNA <- seurat_res
  
  pheatmap::pheatmap(seurat_res, color = viridis::viridis(10))
  
  seurat_cluster <- clustifyr::cor_to_call(seurat_res) %>%
    dplyr::mutate(type = ifelse(r < cor_cutoff, "undetermined", type))
  new_clusters <- seurat_cluster$type
  names(new_clusters) <- seurat_cluster$cluster
  
  colname <- paste0(assay, "_", save_name)
  
  seurat_object[[colname]] <- new_clusters[seurat_object[[clusters]][[1]]]
  
  plot1 <- plotDimRed(seurat_object, col_by = colname,
                      plot_type = plot_type)
  if(ADT){
    plot2 <- plotDimRed(seurat_object, col_by = colname,
                        plot_type = "adt.umap")
    plot3 <- plotDimRed(seurat_object, col_by = colname,
                        plot_type = "wnn.umap")
  }
  
  pdf(file.path(save_dir, "images", paste0("RNA_mapping_", save_name, ".pdf")))
  print(plot1)
  if(ADT){
    print(plot2)
    print(plot3)
  }
  dev.off()
  
  if(ADT){
    # ADT
    clusters <- "ADT_cluster"
    seurat_metadata <- seurat_object[[clusters]]
    
    seurat_metadata[[clusters]] <- as.character(seurat_metadata[[clusters]])
    
    seurat_object[[clusters]] <- seurat_metadata[[clusters]]
    
    # Run clustify
    seurat_res <- clustifyr::clustify(
      input = seurat_mat,
      metadata = seurat_metadata,
      ref_mat = ref_mat,
      query_genes = seurat_genes,
      cluster_col = clusters
    )
    
    return_list$ADT <- seurat_res
    
    seurat_cluster <- clustifyr::cor_to_call(seurat_res) %>%
      mutate(type = ifelse(r < cor_cutoff, "undetermined", type))
    new_clusters <- seurat_cluster$type
    names(new_clusters) <- seurat_cluster$cluster
    
    colname <- paste0("ADT_", save_name)
    
    seurat_object[[colname]] <- new_clusters[seurat_object$ADT_cluster]
    
    plot1 <- plotDimRed(seurat_object, col_by = colname,
                        plot_type = plot_type)
    plot2 <- plotDimRed(seurat_object, col_by = colname,
                        plot_type = "adt.umap")
    plot3 <- plotDimRed(seurat_object, col_by = colname,
                        plot_type = "wnn.umap")
    
    pdf(file.path(save_dir, "images", paste0("ADT_mapping_",
                                             save_name, ".pdf")))
    print(plot1)
    print(plot2)
    print(plot3)
    dev.off()
    
    # Combined
    clusters <- "combined_cluster"
    seurat_metadata <- seurat_object[[clusters]]
    
    seurat_metadata[[clusters]] <- as.character(seurat_metadata[[clusters]])
    
    seurat_object[[clusters]] <- seurat_metadata[[clusters]]
    
    # Run clustify
    seurat_res <- clustifyr::clustify(
      input = seurat_mat,
      metadata = seurat_metadata,
      ref_mat = ref_mat,
      query_genes = seurat_genes,
      cluster_col = clusters
    )
    
    return_list$WNN <- seurat_res
    
    seurat_cluster <- clustifyr::cor_to_call(seurat_res) %>%
      mutate(type = ifelse(r < cor_cutoff, "undetermined", type))
    new_clusters <- seurat_cluster$type
    names(new_clusters) <- seurat_cluster$cluster
    
    colname <- paste0("combined_", save_name)
    
    seurat_object[[colname]] <- new_clusters[seurat_object$combined_cluster]
    
    plot1 <- plotDimRed(seurat_object, col_by = colname,
                        plot_type = plot_type)
    plot2 <- plotDimRed(seurat_object, col_by = colname,
                        plot_type = "adt.umap")
    plot3 <- plotDimRed(seurat_object, col_by = colname,
                        plot_type = "wnn.umap")
    
    pdf(file.path(save_dir, "images", paste0("combined_mapping_",
                                             save_name, ".pdf")))
    print(plot1)
    print(plot2)
    print(plot3)
    dev.off()
    
  }
  return_list$object <- seurat_object
  
  return(return_list)
}

