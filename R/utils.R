#' Confusion matrix
#' 
#' This function will make a confusion matrix based on two vectors. I took this
#' function from the ArchR package.
#' @param i The first vector
#' @param j The second vector
#' @return a matrix of the confusion matrix of the two vectors
#' @importFrom Matrix sparseMatrix
#' @export
#' @examples
#' \dontrun{
#' cm <- confusionMatrix(i, j)
#'}
confusionMatrix <- function(i = NULL, j = NULL){
  ui <- unique(i)
  uj <- unique(j)
  m <- Matrix::sparseMatrix(
    i = match(i, ui),
    j = match(j, uj),
    x = rep(1, length(i)),
    dims = c(length(ui), length(uj))
  )
  rownames(m) <- ui
  colnames(m) <- uj
  m
}

#' Convert between human and mouse
#' 
#' Basic function to convert human to mouse gene name or mouse to
#' human using biomart
#' @param x The list of genes to convert
#' @param convert The type of conversion to do. Can be "human_mouse"
#' to convert human genes to mouse genes or "mouse_human" to convert
#' mouse genes to human
#' @return a dataframe containing both human and mouse genes for the 
#' input vector
#' @export
#' @examples
#' \dontrun{
#' mouse_human_list <- convertHumanGeneList(x = c("INS1", "INS2", "GCG",
#'                                                "ARX", "NKX2-2", "SST"),
#'                                          convert = "human_mouse")
#' mouse_human_list <- convertHumanGeneList(x = c("Ins1", "Ins2", "Gcg",
#'                                                "Arx", "Nkx2-2", "Sst"),
#'                                          convert = "mouse_human")
#'}

convertHumanGeneList <- function(x, convert = "human_mouse"){
  
  check_packages("biomaRt")
  requireNamespace("biomaRt")
  human = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  if(convert == "human_mouse"){
    to_symbol = "mgi_symbol"
    from_symbol = "hgnc_symbol"
    to_mart = mouse
    from_mart = human
  } else if (convert == "mouse_human"){
    to_symbol = "hgnc_symbol"
    from_symbol = "mgi_symbol"
    to_mart = human
    from_mart = mouse
  }
  
  genesV2 = biomaRt::getLDS(
    attributes = c(from_symbol),
    filters = from_symbol,
    values = x ,
    mart = from_mart,
    attributesL = c(to_symbol),
    martL = to_mart,
    uniqueRows=T
  )
  
  return(genesV2)
}

#' Check for packages
#'
#' Checks to ensure all packages required for a function are installed. If
#' any packages aren't installed, the function will quit and print an
#' error message listing packages that need to be installed.
#' @param package_list A list of packages to check
#' @param message A specialized message to print. Otherwise a generic
#' message will be printed
#' @return None
check_packages <- function(package_list,
                           message = NULL) {
  if (is.null(message)) {
    message <- "needed to run this function."
  }
  missing_packages <- c()

  # Check to see if any packages in the list aren't in the namespace
  for (package in package_list) {
    if (!requireNamespace(package, quietly = TRUE)) {
      missing_packages <- c(missing_packages, package)
    }
  }

  # Print error message for missing packages
  if (length(missing_packages) > 0) {
    stop(paste0(
      "The following packages are not installed: ",
      paste(missing_packages, collapse = ", "), 
      " They are ", message, " Please install them."
    ))
  }
}

#' Check Seurat version
#'
#' Determines if the SeuratObject package version is 5.0.0 or higher.
#' @return Logical value indicating if SeuratObject version is 5.0.0 or higher.
#' @importFrom utils packageVersion
is_seurat_v5 <- function() {
  utils::packageVersion("SeuratObject") >= "5.0.0"
}

#' Get Seurat Assay
#'
#' Determines the seurat version and grabs the assay data using
#' the appropriate arguments.
#' @param object A seurat object
#' @param type What type should be returned, data or counts
#' @param assay What assay the data should be pulled from
#' @return The resulting data matrix.
#' @import Seurat
get_seurat_assay <- function(object, 
                             type = "counts",
                              assay = "RNA") {
  if (is_seurat_v5()) {
    GetAssayData(
      object,
      layer = type, 
      assay = assay
    )
  } else {
    GetAssayData(
      object, 
      slot = type, 
      assay = assay
    )
  }
}

#' Set Seurat Assay
#'
#' Determines the seurat version and grabs the assay data using
#' the appropriate arguments.
#' @param object A seurat object
#' @param type What type should be updated, data or counts
#' @param new_data What new data should be added
#' @param assay What assay the data should be added to
#' @return The resulting data matrix.
#' @import Seurat
set_seurat_assay <- function(object, 
                             new_data, 
                             type = "counts",
                             assay = "RNA") {
  if (is_seurat_v5()) {
    SetAssayData(
      object, 
      layer = type, 
      new.data = new_data,
      assay = assay
    )
  } else {
    SetAssayData(
      object, 
      slot = type, 
      new.data = new_data,
      assay = assay
    )
  }
}

#' Get average expression
#'
#' Helper function to calculate and label the average expression
#' @param object A seurat object
#' @param cell_type_col The name of column containing cell type information.
#' Default is `cell_type`.
#' @param assay What assay the data should be pulled from
#' @param ... arguments passed to either AverageExpression (Seurat v4)
#' or PseudobulkExpression (seurat V5)
#' @return A data frame of the average expression of the gene.
#' @import Seurat
get_average_expression <- function(object,
                                   cell_type_col = "cell_type",
                                   assay = "RNA",
                                   ...) {
  # Check seurat version to determine expression function
  if (is_seurat_v5()) {
    avg_expr_df <- data.frame(
      suppressMessages(PseudobulkExpression(object,
        assay = assay,
        method = "average",
        ...
      ))
    )
  } else {
    avg_expr_df <- data.frame(
      suppressMessages(AverageExpression(object, ...)[[assay]])
    )
  }

  return(avg_expr_df)
}