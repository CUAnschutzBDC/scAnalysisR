#' Run differential expression followed by pathview.
#' 
#' This function allows you to run pathview given a matrix of DE output. It can
#' also run pathview on precomputed DE output.
#' @param path_id_list A named list of KEGG path ids. The ids are the 5 number
#' ids, the name will be associated with the output file
#' @param seurat_object OPTIONAL A seurat object that has already been
#' normalized. Must provide either a seurat object or a seurat de.
#' @param seurat_de OPTIONAL the returned object from find_markers. If not 
#' provided, find markers will be run. Must provide either a seurat object or
#' a seurat de.
#' @param seurat_assay OPTIONAL which assay to pull the data from when
#' performing DE. Default is "RNA"
#' @param out_dir OPTIONAL an output directory for the pathview plots. Default
#' it "pathview"
#' @param pval OPTIONAL the pvalue cutoff for the adjusted pvalue from the DE
#' test. Default is 0.05.
#' @param ... OPTIONAL arguments supplied to find_markers
#' @return Nothing is returned, instead pathview plots are made and saved in
#' the specified directory
#' @import tidyverse
#' @export
#' @examples
#' \dontrun{
#' de_to_pathview(path_id_list  = c(NFKB_sig_path = "04064", mapk_signaling = "04010"),
#'                seurat_object = seurat_object))
#' de_to_pathview(path_id_list = c(NFKB_sig_path = "04064", mapk_signaling = "04010"),
#'                seurat_de    = de_output))
#'}

de_to_pathview <- function(path_id_list, seurat_object = NULL,
                           seurat_de = NULL,
                           out_dir = "pathview",
                           pval = 0.05, gen_id = "mmu", ...){

  # Ask user to install pathview to use this function
  if (!requireNamespace("pathview", quietly = TRUE)){
    stop("Package \"pathview\" needed for this function to work. Please install it.",
      call. = FALSE)
  }

  if(is.null(seurat_de) && is.null(seurat_object)){
    stop("must provide either a seurat object or a pre made de list")
  } else if(is.null(seurat_de)){
    marker_genes_rna <- find_markers(seurat_object, ...)
  } else {
    marker_genes_rna <- seurat_de
  }
  # Change to wide form
  marker_genes_rna_long <- marker_genes_rna %>%
    dplyr::filter(p_val_adj < pval) %>%
    tidyr::pivot_wider(names_from = cluster,
                       values_from = "avg_log2FC",
                       id_cols = "gene") %>%
    tibble::column_to_rownames(var = "gene") %>%
    as.matrix()
  
  # Run pathview on all paths
  invisible(lapply(names(path_id_list), function(x) 
    run_pathview(gene_matrix = marker_genes_rna_long,
                 path_id = path_id_list[x],
                 path_name = x,
                 out_dir = out_dir,
                 gen_id = gen_id)
  ))
}

#' Run differential expression 
#' 
#' This function runs differential expression analysis using Seurat's
#' functions. This is called by several functions but can be used on
#' its own as well. Can run a pairwise DE or a one vs all DE.
#' @param seurat_object A seurat object that has already been normalized.
#' @param test_idents OPTIONAL The identities from the metadata used to perform
#' DE. If no identities are provided, the identites already set in the object
#' will be used.
#' @param seurat_assay OPTIONAL which assay to pull the data from when
#' performing DE. Default is "RNA"
#' @param group OPTIONAL to be used when running DE between different samples.
#' The group refers to the column contining the sample info. Default is NULL.
#' @param group_ident_1 OPTIONAL The sample to use as the main sample to test
#' against. Defaults to the first value of the factor
#' @param group_ident_2 OPTIONAL which sample to use as the second group
#' if more than 2 samples are in the seurat object.
#' @return A dataframe containing all markers.
#' @import tidyverse
#' @export
#' @examples
#' \dontrun{
#' marker_df <- find_markers(seurat_object = seurat_object)
#' marker_df <- find_markers(seurat_object = seurat_object,
#'                           group         = "cluster")
#' marker_df <- find_markers(seurat_object = seurat_object,
#'                           group         = "cluster",
#'                           group_ident_1 = 1)
#'}

find_markers <- function(seurat_object,
                         test_idents = NULL, seurat_assay = "RNA",
                         group = NULL, group_ident_1 = NULL,
                         group_ident_2 = NULL){
  if(!is.null(test_idents)){
    Idents(seurat_object) <- test_idents
  }
  if(is.null(group)){
    marker_genes_rna <- FindAllMarkers(seurat_object, assay = seurat_assay)
  } else {
    # Find DE genes between samples per cluster
    if(is.null(group_ident_1)){
      if(is.null(levels(seurat_object[[group]][[1]]))){
        seurat_object[[group]][[1]] <- factor(seurat_object[[group]][[1]])
      }
      group_ident_1 <- levels(seurat_object[[group]][[1]])[1]
    }
    marker_gene_list <- lapply(unique(Idents(seurat_object)), function(x){
      genes <- FindMarkers(seurat_object, ident.1 = group_ident_1,
                           idents.2 = group_ident_2, group.by = group,
                           subset.ident = x)
      genes$cluster <- x
      genes$gene <- rownames(genes)
      return(genes)
    })
    marker_genes_rna <- do.call(rbind, marker_gene_list)
  }
  return(marker_genes_rna)
}

#' Run pathview on single cell data
#' 
#' This function is meant to be called by de_to_pathview and not used on it's own!!!
#' This function allows you to run pathview given a matrix of DE output
#' @param gene_matrix a matrix of log fold change values. Columns are samples
#' and rownames are genes. Can include NA values. Must be in a matrix format.
#' @param path_id the KEGG 5 number id for a pathway.
#' @param path_name OPTIONAL the name of the KEGG pathway. This is used to
#' name the output files. If this isn't set, output files are named based
#' on the path id
#' @param out_dir OPTIONAL the output directory to copy the files into. This
#' directory will be created if it doesn't already exist. Default is "pathways"
#' @param multi_state OPTIONAL if multiple samples are included should they
#' all be plotted on the same output image. Default is FALSE
#' @param gen_id OPTIONAL the species ID. Default is "mmu" for mouse. Set if 
#' you aren't using mouse
#' @keywords internal
#' @import tidyverse
#' @export

run_pathview <- function(gene_matrix, path_id, path_name = NULL,
                         out_dir = "pathways", multi_state = FALSE,
                         gen_id = "mmu"){
  
  # Ask user to install pathview to use this function
  if (!requireNamespace("pathview", quietly = TRUE)){
    stop("Package \"monocle\" needed for this function to work. Please install it.",
      call. = FALSE)
  }

  if(is.null(path_name)){
    path_name = path_id
  }
  
  # Create output directory
  out_dir %>%
    dir.create(showWarnings = F)
  
  pathview_out <- pathview(gene.data = gene_matrix,
                           species = gen_id,
                           pathway.id = path_id,
                           gene.idtype = "SYMBOL",
                           kegg.dir = out_dir,
                           multi.state = multi_state,
                           match.data = multi_state,
                           low = list(gene = "#225ea8", cpd = "blue"),
                           mid = list(gene = "white", cpd = "gray"),
                           high = list(gene = "#e31a1c", cpd = "yellow"),
                           na.col = "#bdbdbd")
  
  if(multi_state){
    orig_name <- str_c(gene_id, path_id, ".pathview.png")
    new_name <- str_c(out_dir, "/", gen_id, path_name,
                      ".pathview.png")
    file.rename(orig_name, new_name)
  }
  invisible(lapply(colnames(gene_matrix), function(cluster){
    orig_name <- stringr::str_c(gen_id, path_id, ".pathview.", cluster, ".png")
    new_name <- stringr::str_c(out_dir, "/", gen_id, "_", path_name,
                               ".pathview.", cluster, ".png")
    tryCatch(
      {
        file.rename(orig_name, new_name)
      },
      warning = function(cond){
        output_file <- stringr::str_c(gen_id, path_id, ".pathview.png")
        if(grepl("cannot rename file", cond) && file.exists(output_file)){
          file.remove(output_file)
        }
      }
    )
    
  }))
  
}

#' Run gost gene ontology
#' 
#' This function allows you to run gost given a matrix of DE output or seurat
#' object
#' @param seurat_de OPTIONAL A de data frame returned by find_markers.
#' Must provide either a seurat object or seurat de
#' @param seurat_object OPTIONAL A seurat object that has already been
#' normalized. Must provide either a seurat object or a seurat de.
#' @param sources OPTIONAL a list of what sources to plot. Can be any sources
#' returned by gost. Default is GO:BP, KEGG, REAC and TF
#' @param plot_colors OPTIONAL colors for the samples. Default is blue to red.
#' @param intersection_cutoff OPTIONAL how many intersecting genes there must
#' be to plot. Default is 5.
#' @return A list containing two items, gost_output and go_plots
#' @import tidyverse
#' @export
#' @examples
#' \dontrun{
#' go_list <- run_gost(seurat_object = seurat_object)
#' go_list <- run_gost(seurat_de = de_df)
#'}

run_gost <- function(seurat_de = NULL, seurat_object = NULL,
                     sources = c("GO:BP", "KEGG", "REAC", "TF"),
                     plot_colors = c("blue", "red"),
                     intersection_cutoff = 5, ...){

  # Ask user to install pathview to use this function
  if (!requireNamespace("gprofiler2", quietly = TRUE)){
    stop("Package \"gprofiler2\" needed for this function to work. Please install it.",
      call. = FALSE)
  }

  if(is.null(seurat_de) && is.null(seurat_object)){
    stop("must provide either a seurat object or a pre made de list")
  } else if(is.null(seurat_de)){
    marker_genes_rna <- find_markers(seurat_object, ...)
  } else {
    marker_genes_rna <- seurat_de
  }
  marker_genes_gost <- marker_genes_rna %>%
    dplyr::filter(p_val_adj < pval)
  
  gene_list <- lapply(unique(marker_genes_gost$cluster), function(x){
    marker_genes_short <- dplyr::filter(marker_genes_gost, cluster == x) %>%
      dplyr::arrange(p_val_adj)
    return(marker_genes_short$gene)
  })
  
  names(gene_list) <- unique(marker_genes_gost$cluster)
  
  gost_output <- gprofiler2::gost(query = gene_list,
                                  organism = "mmusculus",
                                  multi_query = FALSE,
                                  ordered_query = FALSE,
                                  user_threshold = pval,
                                  custom_bg = NULL,
                                  correction_method = "fdr")
  # Separate this out into it's own function I think, a function where you
  # only pass it the gost_output and sources.
  go_plots <- make_go_plots(gost_output = gost_output,
                            sources = sources,
                            plot_colors = plot_colors,
                            intersection_cutoff = intersection_cutoff)
  return(list(gost_output = gost_output, go_plots = go_plots))
}

#' Makes go plots
#'
#' This function is meant to be called by run_gost and not used on it's own!!!
#' This function makes plots based on different go terms after running gost. 
#' This function is primarily meant to be called by other functions.
#' @param gost_output the output from running gost.
#' @param sources OPTIONAL a list of what sources to plot. Can be any sources
#' returned by gost. Default is GO:BP, KEGG, REAC and TF
#' @param plot_colors OPTIONAL colors for the samples. Default is blue to red.
#' @param intersection_cutoff OPTIONAL how many intersecting genes there must
#' be to plot. Default is 5.
#' @keywords internal
#' @import tidyverse
#' @export

make_go_plots <- function(gost_output,
                          sources = c("GO:BP", "KEGG", "REAC", "TF"),
                          plot_colors = c("blue", "red"),
                          intersection_cutoff = 5){
  all_plots <- lapply(unique(gost_output$result$query), function(query){
    source_plots <- lapply(sources,
                           function(source) make_go_plot_single(
                             gost_output = gost_output,
                             gost_query = query,
                             gost_source = source,
                             plot_colors = plot_colors,
                             intersection_cutoff = intersection_cutoff)
    )
    names(source_plots) <- sources
    return(source_plots)
  })
  names(all_plots) <- unique(gost_output$result$query)
  return(all_plots)
}

#' Makes go plots of different terms for one DE test
#' 
#' This function is meant to be called by make_go_plots and not used on it's own!!!
#' This function makes plots based on different go terms after running gost. 
#' This function is primarily meant to be called by other functions.
#' @param gost_output the output from running gost.
#' @param gost_query the name of the cell type or sample from DE to plot
#' @param gost_source what source to plot. Can be any sources returned by gost
#' @param plot_colors OPTIONAL colors for the samples. Default is blue to red.
#' @param intersection_cutoff OPTIONAL how many intersecting genes there must
#' be to plot. Default is 5.
#' @keywords internal
#' @import tidyverse
#' @export

make_go_plot_single <- function(gost_output, gost_query, gost_source,
                                plot_colors = c("blue", "red"),
                                intersection_cutoff = 5){
  gost_output_one <- gost_output$result %>%
    dplyr::filter(query == gost_query & source == gost_source &
                    intersection_size >= intersection_cutoff) %>%
    dplyr::distinct(term_name, .keep_all = TRUE) %>%
    dplyr::top_n(-20, wt = p_value) %>%
    dplyr::arrange(precision)
  
  if(nrow(gost_output_one) > 20){
    gost_output_one <- gost_output_one[1:20,]
  }
  
  gost_output_one$term_name <- factor(gost_output_one$term_name,
                                      levels = gost_output_one$term_name)
  
  gost_plot <- ggplot2::ggplot(gost_output_one,
                               ggplot2::aes(x = precision,
                                            y = term_name,
                                            color = -log10(p_value),
                                            size = intersection_size)) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradientn(
      colors=grDevices::colorRampPalette(plot_colors)(n = 299)) +
    ggplot2::xlab("GeneRatio") +
    ggplot2::ylab(paste0(gost_source, " term"))
  return(gost_plot)
}

#' Save output from running gost
#' 
#' This function saves both plots and text output from running gost
#' @param gost_output the output from running run_gost function.
#' @param save_dir_plots the directory to save the plots to
#' @param save_dir_text OPTIONAL the directory to save the text to. If left
#' NULL, will default to the same directory as the plots
#' @param save_plots OPTIONAL if plots should be saved. Default is TRUE
#' @param save_text OPTIONAL if text should be saved to a csv file. The whole
#' gost output will be saved. Default is TRUE.
#' @param save_excel OPTIONAL if an excel file should be created. A file will
#' be created for each query and a tab will be inserted for each type of output.
#' Requres openxlsx is installed.
#' @return No return, writes either a csv or excel workbook of gene onotology results
#' and saves plots to a pdf
#' @import tidyverse
#' @export
#' @examples
#' \dontrun{
#' save_gost(gost_output    = gost_output,
#'           save_dir_plots = "plots")
#' save_gost(gost_output    = gost_output,
#'           save_dir_plots = "plots",
#'           save_dir_text  = "files",
#'           save_excel     = FALSE)
#' }

save_gost <- function(gost_output, save_dir_plots, save_dir_text = NULL,
                      save_plots = TRUE, save_text = TRUE, save_excel = TRUE){
  # Save plots
  if(save_plots){
    # Create output directory
    save_dir_plots %>%
      dir.create(showWarnings = F)
    gost_plots <- gost_output$go_plots
    invisible(lapply(names(gost_plots), function(x){
      pdf(file.path(save_dir_plots, paste0(x, "GSE.pdf")))
      gost_plots[[x]]
      dev.off()
    }))
  }
  # Save to plot directory if text directory isn't provided
  if(is.null(save_dir_text)){
    save_dir_text <- save_dir_plots
  }
  gost_text <- gost_output$gost_output$result
  gost_text$parents <- NULL
  # Save all results to csv
  if(save_text){
    # Create output directory
    save_dir_text %>%
      dir.create(showWarnings = F)
    write.csv(gost_text,
              file.path(save_dir_text, "all_GSE_results.csv"))
  }
  if(save_excel){
    if("openxlsx" %in% rownames(installed.packages()) == FALSE){
      install.packages("openxlsx")
    }
    # Create output directory
    save_dir_text %>%
      dir.create(showWarnings = F)
    invisible(lapply(unique(gost_text$query), function(gost_query){
      
      # Create excel wb
      gene_wb <- createWorkbook()
      
      # Write to excel wb
      full_list <- lapply(unique(gost_text$source), function(gost_source){
        new_df <- gost_text %>%
          dplyr::filter(source == gost_source & query == gost_query)
        gost_source_write <- sub(":", "_", gost_source)
        addWorksheet(gene_wb, gost_source_write)
        writeData(gene_wb, gost_source_write, new_df)
      })
      
      ## Save workbook to working directory
      saveWorkbook(gene_wb,
                   file = file.path(save_dir_text, paste0(gost_query,
                                 "GSE_results.xlsx")),
                   overwrite = TRUE)
    }))
  }
}

#' Run pairwise differential expression 
#' 
#' This function runs  pairwise differential expression analysis on all 
#' possible combinations in the meta data column provided using seurat's
#' FindMarkers
#' @param seurat_object A seurat object that has already been normalized.
#' @param meta_col The column in the meta data to use to identify all 
#' pairwise markers
#' @param p_val_cutoff OPTIONAL The cutoff for p value for returned genes.
#' Default is 0.05
#' @param ... OPTIONAL extra arguments passed to FindMarkers
#' @return A dataframe containing all pairwise markers.
#' @import tidyverse
#' @export
#' @examples
#' \dontrun{
#' marker_df <- pairwise_markers(seurat_object = seurat_object,
#'                               meta_col      = "clusters")
#'}

pairwise_markers <- function(seurat_object, meta_col,
                             p_val_cutoff = 0.05, ...){
  Idents(seurat_object) <- meta_col
  combinations <- combn(unique(Idents(seurat_object)), m = 2)
  all_de <- lapply(1:ncol(combinations), function(x){
    ident1 <- combinations[1, x]
    ident2 <- combinations[2, x]
    marker_genes <- FindMarkers(seurat_object, only.pos = FALSE,
                                ident.1 = ident1, ident.2 = ident2, ...)

    # Determine what cluster was idenfied as "up" and what cluster was "down"
    marker_genes_up <- marker_genes %>%
      dplyr::filter(avg_log2FC > 0) %>%
      dplyr::mutate(cluster_up = ident1, cluster_down = ident2)
    marker_genes_down <- marker_genes %>%
      dplyr::filter(avg_log2FC < 0) %>%
      dplyr::mutate(cluster_up = ident2, cluster_down = ident1) %>%
      dplyr::mutate(avg_log2FC = abs(avg_log2FC))
    marker_genes <- rbind(marker_genes_up, marker_genes_down)
    marker_genes$gene_name <- rownames(marker_genes)
    return(marker_genes)
  })
  all_de <- do.call(rbind, all_de)
  
  all_de <- all_de %>%
    dplyr::filter(p_val_adj < 0.05)
  
  return(all_de)
  
}

#' Find and write markers
#' 
#' This function finda all markers of all clusters for a specific column
#' in your seurat metadata by running FindAllMarkers (for the default) or
#' by running FindMarkers in a pairwise fashion (when pairwise is TRUE).
#' It writes these markers to both a csv file and an excel document with
#' one sheet per group from the meta data (ex 1 sheet per cluster). This
#' function will further run and write the output of hypergeometric tests
#' on supplied gene lists using the identified markers. It will perform a 
#' p-value fdr correction based on the number of lists supplied.
#' @param seurat_object A seurat object that has already been normalized.
#' @param save_dir The path to save all of the output file.
#' @param meta_col OPTIONAL The column in the meta data to use to identify  
#' markers. Default is RNA_cluster.
#' @param assay OPTIONAL the assay to run the differential expression.
#' Default is RNA.
#' @param p_val OPTIONAL The cutoff for p value for returned genes.
#' Default is 0.05
#' @param logfc OPTIONAL The cutoff for the log fold change for returned
#' genes. Default is 0.5.
#' @param gene_lists OPTIONAL gene lists to compare against differentially
#' expressed genes. Can provide many vectors of genes as a list. Provided
#' lists will be used to run a hypergeometric test against differentially
#' expressed genes identified from each cluster. If no gene lists are 
#' included, no hypergeometric tests will be run.
#' @param pairwise OPTIONAL if the differential expression should be done
#' in a pairwise fashion for all groups vs all groups (TRUE) or in a one vs
#' all fashion (FALSE). FindMarkers from Seurat defaults to the one vs all
#' approach. Default is FALSE.
#' @return A dataframe containing all markers. It will also write a csv file
#' with the DE genes, a csv file with the hypergeometric test output (if run)
#' and an excel file with both hypergeometric output and DE genes.
#' @import tidyverse
#' @import openxlsx
#' @export
#' @examples
#' \dontrun{
#' marker_df <- find_write_markers(seurat_object = seurat_object,
#'                                 save_dir      = "DE_files",
#'                                 meta_col      = "clusters",
#'                                 pairwise      = TRUE)
#'
#' # This is better with longer gene lists, but this is just an example
#' my_genes <- list(t_cell = c("Cd3d", "Cd4", "Cd8"),
#'                  beta   = c("Ins1", "Ins2", "Nkx2-2"),
#'                  alpha  = c("Gcg", "Arx"))
#' marker_df <- find_write_markers(seurat_object = seurat_object,
#'                                 save_dir      = "DE_files",
#'                                 meta_col      = "clusters",
#'                                 gene_lists    = my_genes)
#'}

find_write_markers <- function(seurat_object, save_dir,
                               meta_col = "RNA_cluster",
                               assay = "RNA", pval = 0.05,
                               logfc = 0.5, gene_lists = NULL,
                               pairwise = FALSE){
  # Create a directory for results
  ifelse(!dir.exists(file.path(save_dir, "files", "DE")),
         dir.create(file.path(save_dir, "files", "DE")), FALSE)
  
  Idents(seurat_object) <- meta_col
  
  if (pairwise){
    # Find markers for all clusters in a pairwise fashion
    marker_genes <- find_write_markers_pairwise(seurat_object = seurat_object,
                                                save_dir = save_dir,
                                                meta_col = meta_col,
                                                assay = assay,
                                                pval = pval,
                                                logfc = logfc,
                                                gene_lists = gene_lists)
  } else {
    # Find markers using a one vs all approach
    marker_genes <- find_write_markers_orig(seurat_object = seurat_object,
                                            save_dir = save_dir,
                                            meta_col = meta_col,
                                            assay = assay,
                                            pval = pval,
                                            logfc = logfc,
                                            gene_lists = gene_lists)
  }
  
  return(marker_genes)
}  


#' Find and write markers
#' 
#' This function finda all markers of all clusters for a specific column
#' in your seurat metadata by running FindAllMarkers. It writes these 
#' markers to both a csv file and an excel document with one sheet per
#' group from the meta data (ex 1 sheet per cluster). This function will
#' further run and write the output of hypergeometric tests on supplied
#' gene lists using the identified markers. It will perform a p-value fdr
#' correction based on the number of lists supplied. This function is meant
#' to be called by find_write_markers and not used directly.
#' @param seurat_object A seurat object that has already been normalized.
#' @param save_dir The path to save all of the output file.
#' @param meta_col OPTIONAL The column in the meta data to use to identify  
#' markers. Default is RNA_cluster.
#' @param assay OPTIONAL the assay to run the differential expression.
#' Default is RNA.
#' @param p_val OPTIONAL The cutoff for p value for returned genes.
#' Default is 0.05
#' @param logfc OPTIONAL The cutoff for the log fold change for returned
#' genes. Default is 0.5.
#' @param gene_lists OPTIONAL gene lists to compare against differentially
#' expressed genes. Can provide many vectors of genes as a list. Provided
#' lists will be used to run a hypergeometric test against differentially
#' expressed genes identified from each cluster. If no gene lists are 
#' included, no hypergeometric tests will be run.
#' @return A dataframe containing all markers. It will also write a csv file
#' with the DE genes, a csv file with the hypergeometric test output (if run)
#' and an excel file with both hypergeometric output and DE genes.
#' @import tidyverse
#' @import openxlsx
#' @export

find_write_markers_orig <- function(seurat_object, save_dir,
                                    meta_col = "RNA_cluster",
                                    assay = "RNA", pval = 0.05,
                                    logfc = 0.5, gene_lists = NULL) {
  
  marker_genes <- FindAllMarkers(seurat_object, assay = seurat_assay,
                                 only.pos = TRUE)
  write.csv(marker_genes, file = file.path(save_dir, "files",
                                        "DE", paste0(assay, "_markers_",
                                        meta_col, ".csv")))
  
  # Create excel wb
  gene_wb <- createWorkbook()
  
  # Write to excel wb
  full_list <- lapply(unique(marker_genes$cluster), function(x){
    x <- as.character(x)
    new_df <- marker_genes %>%
      dplyr::filter(cluster == x & p_val_adj < pval & avg_log2FC > logfc)
    addWorksheet(gene_wb, x)
    writeData(gene_wb, x, new_df)
  })
  
  # Run a hypergeomitric test
  if(!is.null(gene_lists)){
    hypergeometric <- hypergeometric_test(seurat_object = seurat_data,
                                          gene_list = gene_lists,
                                          DE_table = marker_genes,
                                          DE_p_cutoff = 0.05,
                                          DE_lfc_cutoff = 0.5,
                                          correction_method = "fdr")
    write.csv(hypergeometric, file = file.path(save_dir, "files",
                                            "DE", paste0(assay,
                                            "_hypergeometric_",
                                            meta_col, ".csv")))
    
    # Write to excel wb
    full_list <- lapply(unique(hypergeometric$cluster), function(x){
      x <- as.character(x)
      new_df <- hypergeometric %>%
        dplyr::filter(cluster == x)
      worksheet_name <-  paste0(x, "_gse")
      addWorksheet(gene_wb, worksheet_name)
      writeData(gene_wb, worksheet_name, new_df)
    })
  }
  
  ## Save workbook to working directory
  saveWorkbook(gene_wb,
               file = file.path(save_dir,
                             "files", "DE", paste0(assay, "_markers_",
                             meta_col, ".xlsx")),
               overwrite = TRUE)
  
  return(marker_genes)
}

#' Find and write markers
#' 
#' This function finda all markers of all clusters for a specific column
#' in your seurat metadata by running FindMarkers in a pairwise fashion.
#' It writes these markers to both a csv file and an excel document with
#' one sheet per group from the meta data (ex 1 sheet per cluster). This
#' function will further run and write the output of hypergeometric tests
#' on supplied gene lists using the identified markers. It will perform a 
#' p-value fdr correction based on the number of lists supplied. This
#' function is meant to be called by find_write_markers and not used directly.
#' @param seurat_object A seurat object that has already been normalized.
#' @param save_dir The path to save all of the output file.
#' @param meta_col OPTIONAL The column in the meta data to use to identify  
#' markers. Default is RNA_cluster.
#' @param assay OPTIONAL the assay to run the differential expression.
#' Default is RNA.
#' @param p_val OPTIONAL The cutoff for p value for returned genes.
#' Default is 0.05
#' @param logfc OPTIONAL The cutoff for the log fold change for returned
#' genes. Default is 0.5.
#' @param gene_lists OPTIONAL gene lists to compare against differentially
#' expressed genes. Can provide many vectors of genes as a list. Provided
#' lists will be used to run a hypergeometric test against differentially
#' expressed genes identified from each cluster. If no gene lists are 
#' included, no hypergeometric tests will be run.
#' @return A dataframe containing all markers. It will also write a csv file
#' with the DE genes, a csv file with the hypergeometric test output (if run)
#' and an excel file with both hypergeometric output and DE genes.
#' @import tidyverse
#' @import openxlsx
#' @export

find_write_markers_pairwise <- function(seurat_object, save_dir,
                                        meta_col = "RNA_cluster",
                                        assay = "RNA", pval = 0.05,
                                        logfc = 0.5, gene_lists = NULL) {
  
  marker_genes <- pairwise_markers(seurat_object, assay = seurat_assay,
                                   meta_col = meta_col)
  write.csv(marker_genes, file = file.path(save_dir, "files",
                                           "DE", paste0(assay,
                                                        "_pairwise_markers_",
                                                        meta_col, ".csv")))
  
  # Create excel wb
  gene_wb <- createWorkbook()
  
  # Write to excel wb
  values <- unique(c(marker_genes$cluster_down, marker_genes$cluster_up))
  combinations <- combn(unique(c(marker_genes$cluster_down,
                                 marker_genes$cluster_up)),
                        m = 2)
  full_list <- lapply(1:ncol(combinations), function(x){
    ident1 <- combinations[1, x]
    ident2 <- combinations[2, x]
    sheet_name <- paste0(ident1, "_vs_", ident2)
    new_df <- marker_genes %>%
      dplyr::filter((cluster_up == ident1 & cluster_down == ident2) |
                      (cluster_up == ident2 & cluster_down == ident1)) %>%
      dplyr::filter(p_val_adj < pval & avg_log2FC > logfc)
    addWorksheet(gene_wb, sheet_name)
    writeData(gene_wb, sheet_name, new_df)
  })
  
  # Alter so it will work with a hypergeometric test
  marker_genes$cluster <- paste0(marker_genes$cluster_up, "_vs_",
                                 marker_genes$cluster_down)
  
  # Run a hypergeomitric test
  if(!is.null(gene_lists)){
    hypergeometric <- hypergeometric_test(seurat_object = seurat_data,
                                          gene_list = gene_lists,
                                          DE_table = marker_genes,
                                          DE_p_cutoff = 0.05,
                                          DE_lfc_cutoff = 0.5,
                                          correction_method = "fdr")
    write.csv(hypergeometric, file = file.path(save_dir, "files",
                                               "DE", paste0(assay,
                                                            "_hypergeometric_",
                                                            meta_col, ".csv")))
    
    # Write to excel wb
    full_list <- lapply(unique(hypergeometric$cluster), function(x){
      x <- as.character(x)
      new_df <- hypergeometric %>%
        dplyr::filter(cluster == x)
      worksheet_name <-  paste0(x, "_gse")
      addWorksheet(gene_wb, worksheet_name)
      writeData(gene_wb, worksheet_name, new_df)
    })
  }
  
  ## Save workbook to working directory
  saveWorkbook(gene_wb,
               file = file.path(save_dir,
                                "files", "DE",
                                paste0(assay, "_pairwise_markers_",
                                                      meta_col, ".xlsx")),
               overwrite = TRUE)
  
  return(marker_genes)
}

#' Runs a hypergeometric test
#' 
#' Given a seurat object, results from a DE test, and a set of gene lists,
#' this runs a hypergeometric test (fishers exact test) to identify significant
#' enrichment of the provided gene lists in the markers of clusters in the 
#' seurat object. All genes present in the Seurat object are used as the background
#' gene set. This function is called by find_write_markers, but can also
#' be used independently.
#' @param seurat_object A seurat object that has already been normalized.
#' @param gene_list gene lists to compare against differentially
#' expressed genes. Can provide many vectors of genes as a list. Provided
#' lists will be used to run a hypergeometric test against differentially
#' expressed genes identified from each cluster.
#' @param DE_table A data frame with differential expression results that can
#' be made by running FindAllMarkers
#' @param DE_p_cutoff OPTIONAL The cutoff for p value for differentially expressed
#' genes to be used in the hypergeometric test. Default is 0.05
#' @param DE_lfc_cutoff OPTIONAL The cutoff for the log fold change for differentially
#' expressed genes to be used in the hypergeometric test. Default is 0.5.
#' @param correction_method OPTIONAL The method used to correct the p-values for
#' multiple correction. Default is "fdr". Can be any method used by p.adjust.
#' @return A dataframe containing hypergeometric output. Includes the expected number
#' of genes from each set seen, the enrichment of the gene set, a p-value, and a
#' corrected p-value for each cluster and each gene set.
#' @import tidyverse
#' @export
#' @examples
#' \dontrun{
#' # This is better with longer gene lists, but this is just an example
#' my_genes <- list(t_cell = c("Cd3d", "Cd4", "Cd8"),
#'                  beta   = c("Ins1", "Ins2", "Nkx2-2"),
#'                  alpha  = c("Gcg", "Arx"))
#' DE_res <- FindMarkers(seurat_object, only.pos = TRUE)
#' enrichment <- hypergeometric_test(seurat_object = seurat_object,
#'                                   gene_list     = my_list,
#'                                   DE_table      = DE_res)
#'}

hypergeometric_test <- function(seurat_object, gene_list, DE_table,
                                DE_p_cutoff = 0.05, DE_lfc_cutoff = 0.5,
                                correction_method = "fdr"){
  # Pull out gene list
  hypergeometric_list <- lapply(names(gene_list), function(list_name){
    # Pull out one gene list
    gene_list_one <- gene_list[[list_name]]
    
    DE_list <- lapply(unique(DE_table$cluster), function(cluster_name){
      # Pull out one DE test and only sig genes
      DE_one <- DE_table %>%
        dplyr::filter(cluster == cluster_name &
                        p_val_adj < DE_p_cutoff &
                        avg_log2FC > DE_lfc_cutoff)
      
      DE_genes <- unique(DE_one$gene)

      # Find number of overlaps
      x <- length(intersect(gene_list_one, DE_genes))
      
      # Length of gene list that overlaps with gene in object (total possible genes 
      # to see in comparison)
      m <- length(intersect(gene_list_one, rownames(seurat_object)))
      
      # All genes from object not in list
      n <- length(setdiff(rownames(seurat_object), gene_list_one))
      
      # Total number of genes
      total <- nrow(seurat_object)
      
      # Length of the DE list
      k <- length(DE_genes)
      
      # Calculated expected number of genes
      expected_num <- (m*k)/total
      
      # Calcluate the representation factor
      representation <- x/expected_num
      
      # Calculate the p_val
      p_val <- sum(dhyper(x, m, n, k))
      
      return_df <- data.frame(cluster = cluster_name,
                              gene_list = list_name,
                              overlaps = x,
                              gene_list_len = m,
                              DE_length = k,
                              background_len = n,
                              total_genes = total,
                              expected_overlap = expected_num,
                              overrepresentation = representation,
                              p_val = p_val)
      
      return(return_df)
    })
    full_return <- do.call(rbind, DE_list)
    return(full_return)
  })
  
  full_hypergeometric <- do.call(rbind, hypergeometric_list)
  
  # Calculate adjusted p values
  full_hypergeometric$p_adj <- p.adjust(full_hypergeometric$p_val,
                                        method = correction_method)
  
  return(full_hypergeometric)
}

#' Plots results of enrichment
#' 
#' This function takes the output of hypergeometric_test and creates
#' a heatmap based on the log of the adjusted p-values.
#' @param hypergeom_output The output from running hypergeometric_test
#' @param colors OPTIONAL The colors to use to annotate the x axis. If left
#' blank and no meta_df and color list are provided, will default to Set 1
#' from RColorBrewer
#' @param meta_df OPTIONAL A dataframe with meta data to use to color the x axis
#' of the plot. If left blank, the cluster information will be used from the
#' hypergeometric test.
#' @param color_list OPTIONAL A list of color vectors to color the x axis. The names
#' of each element of the list much match the column names of the meta_df and the
#' length of these vectors must match the levels of each colunm of the meta_df. This
#' is used in conjunction with pheatmap, this tutorial can help explain how the meta_df
#' and color_list are used
#' https://davetang.org/muse/2018/05/15/making-a-heatmap-in-r-with-the-pheatmap-package/
#' must be included if meta_df is specified.
#' @param cluster_rows OPTIONAL if the heatmap rows should be clusterd. Default is FALSE
#' @param cluster_cols OPTIONAL if the heatmap columns should be clustered. Default is
#' FALSE.
#' @param color_palette OPTIONAL The color of the heatmap. Can be "blueRed" or your
#' own palette. If no colors are provided, will use the blueYellow palette from 
#' the ArchR package.
#' @param breaks OPTIONAL How to adjust the color scale. Can use if the heatmap looks
#' washed out. Default is FALSE.
#' @param mak_val OPTIONAL where to cutoff the color scale. Default is no cutoff.
#' @return A pheatmap object
#' @import tidyverse
#' @import RColorBrewer
#' @import pheatmap
#' @export
#' @examples
#' \dontrun{
#' # This is better with longer gene lists, but this is just an example
#' my_genes <- list(t_cell = c("Cd3d", "Cd4", "Cd8"),
#'                  beta   = c("Ins1", "Ins2", "Nkx2-2"),
#'                  alpha  = c("Gcg", "Arx"))
#' DE_res <- FindMarkers(seurat_object, only.pos = TRUE)
#' enrichment <- hypergeometric_test(seurat_object = seurat_object,
#'                                   gene_list     = my_list,
#'                                   DE_table      = DE_res)
#' heatmap_out <- plot_hypergeom(enrichment)
#'}

plot_hypergeom <- function(hypergeom_output, colors = NULL, meta_df = NULL,
                           color_list = NULL,
                           cluster_rows = FALSE, cluster_cols = FALSE,
                           color_palette = NULL,
                           breaks = FALSE, max_val = NULL){
  
  hypergeom_output$log_adj_pval <- -log10(hypergeom_output$p_adj)
  
  hypergeom_output_w <- hypergeom_output %>%
    dplyr::select(c(cluster, gene_list, log_adj_pval)) %>%
    tidyr::pivot_wider(names_from = gene_list, values_from = log_adj_pval) %>%
    base::data.frame()
  
  rownames(hypergeom_output_w) <- hypergeom_output_w$cluster
  hypergeom_output_w$cluster <- NULL
  hypergeom_output_w <- t(hypergeom_output_w)
  
  if(is.null(meta_df)){
    # Make a df for the column labeling
    sample_info <- data.frame(cluster = colnames(hypergeom_output_w))
    rownames(sample_info) <- sample_info$cluster
    # Add levels
    if(is.null(levels(sample_info$cluster))){
      sample_info$cluster <- factor(sample_info$cluster)
    }
    if(is.null(colors)){
      colors <- brewer.pal(length(levels(sample_info$cluster)), "Set1")
      names(colors) <- levels(sample_info$cluster)
    } 
    # make a list for the column labeing
    coloring <- list(colors)
    names(coloring) <- "cluster"
  } else {
    # The sample info and color list must be provided
    sample_info <- meta_df
    coloring <- color_list
  }
  
  # Set cluster order
  
  cluster_order <- levels(sample_info$cluster)
  # Colors for heatmap (from the ArchR package)
  if(is.null(color_palette)){
    color_palette <- c("#352A86", "#343DAE", "#0262E0", "#1389D2", "#2DB7A3",
                       "#A5BE6A", "#F8BA43", "#F6DA23", "#F8FA0D")
  } else if (color_palette == "blueRed") {
    pal <- colorRampPalette(c("blue", "white", "red"))
    color_palette <- pal(30)
  }
  
  
  if(!cluster_cols){
    sample_info <- sample_info[order(match(sample_info$cluster,
                                           cluster_order)), , drop = FALSE]
    if(!identical(colnames(hypergeom_output_w), rownames(sample_info))){
      hypergeom_output_w <- hypergeom_output_w[ , rownames(sample_info)]
    }
  }
  
  if((breaks)){
    quantile_breaks <- function(xs, n = 30) {
      breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
      breaks[!duplicated(breaks)]
    }
    
    breaks <- quantile_breaks(hypergeom_output_w, n = 30)
  } else {
    breaks <- NULL
  }
  
  if(!is.null(max_val)){
    hypergeom_output_w <- ifelse(hypergeom_output_w > max_val,
                                 max_val, hypergeom_output_w)
  }

  heatmap <- pheatmap(hypergeom_output_w, cluster_rows = cluster_rows,
                      cluster_cols = cluster_cols,
                      show_rownames = TRUE, 
                      show_colnames = TRUE, annotation_col = sample_info,
                      annotation_colors = coloring, color = color_palette,
                      border_color = NA, clustering_method = "complete",
                      silent = TRUE, breaks = breaks)
  return(heatmap)
}

