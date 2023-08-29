library(scAnalysisR)
library(Seurat)
library(tidyverse)

data_dir <- "~/Documents/Analysis/r_packages/scAnalysisR/inst/extdata"
save_dir <- "~/Documents/Analysis/r_packages/scAnalysisR/data"

# Read in data ---------------------------------------------------------------
pbmc_so <- create_seurat_object(sample = "pbmc_data",
                                count_path = data_dir,
                                ADT = TRUE,
                                hashtag = FALSE)

# Normalize --------------------------------------------------------------------
pbmc_so <- pbmc_so %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData()

pbmc_so <- SCTransform(pbmc_so, verbose = FALSE)

DefaultAssay <- "RNA"

# Dimensional reductions -------------------------------------------------------

# RNA


# SCT


# Name populations -------------------------------------------------------------


# Save data --------------------------------------------------------------------
save(pbmc_so, file = file.path(save_dir, "pbmc_so.rda"))
