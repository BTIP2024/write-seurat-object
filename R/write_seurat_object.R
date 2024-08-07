#' Load single-cell data and create a seurat object
#'
#' Supported formats are .gz and .h5 files.
#' This requires the SeuratDisk package available only via github but can be installed
#' with devtools::install_github
#' 
#' Sample file for test run: https://cellgeni.cog.sanger.ac.uk/cvidcellatlas/adata_SS2_for_download.h5ad
#' 
#' If using matrix data, make sure the input is the folder containing the matrix, features, and barcodes
#' 
#' @param input is the single-cell data file
#' @examples
#' write_seurat_object("file.h5")
#' @export
write_seurat_object <- function(input) {
  if(grepl('.gz', input)){
    obj <- Seurat::Read10X(input)
    saveRDS(obj, file = "seurat_object.rds")
  } else {
    h5_data <- Seurat::Read10X_h5(filename = input, use.names = TRUE, unique.features = TRUE)
    
    if("Gene Expression" %in% str(h5_data)) {
      cts <- h5_data$'Gene Expression'
      obj <- Seurat::CreateSeuratObject(counts = cts, project = "Seurat", min.cells = 3, min.features = 200)
    } else {
      obj <- Seurat::CreateSeuratObject(counts = h5_data, project = "Seurat", min.cells = 3, min.features = 200)
    }
  saveRDS(obj, file = "seurat_object.rds")
  }
}

