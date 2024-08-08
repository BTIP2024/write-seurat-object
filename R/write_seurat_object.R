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
write_seurat_object <- function(input){
   if(grepl('.h5', input)) {
      h5 <- hdf5r::H5File$new(input)
      obj <- Seurat::CreateSeuratObject(counts = h5, project = "Seurat", min.cells = 3, min.features = 200)
   } else if(grepl('.gz', input)){
      utils::untar(input, files = NULL, list = FALSE, exdir = "sample")
      new <- list.files(path = "./sample/", pattern = "*.gz", recursive = TRUE, full.names = T) #WORKED
      
      file.copy(from = new, to = ".", overwrite = T)
      
      obj <- Seurat::Read10X(data.dir = ".")
      saveRDS(obj, file = "seurat_object.rds")
   }
}

