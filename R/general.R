#' Confusion matrix
#' 
#' This function will make a confusion matrix based on two vectors. I took this
#' function from the ArchR package.
#' @param i The first vector
#' @param j The second vector
#' @return a matrix of the confusion matrix of the two vectors
#' @import Matrix
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
#' @import biomaRt
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
  
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
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
  
  genesV2 = getLDS(attributes = c(from_symbol), filters = from_symbol,
                   values = x , mart = from_mart, attributesL = c(to_symbol),
                   martL = to_mart, uniqueRows=T)
  
  return(genesV2)
}
