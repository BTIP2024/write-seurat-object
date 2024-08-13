#' Load single-cell data and create a seurat object
#'
#' Supported formats are .gz, .h5, .h5ad, and .loom files.
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
      utils::untar(input, files = NULL, list = FALSE, exdir = "sample")
      new <- list.files(path = "./sample/", pattern = "*.gz", recursive = TRUE, full.names = T) #WORKED
      
      file.copy(from = new, to = ".", overwrite = T)
      
      obj <- Seurat::Read10X(data.dir = ".")
   } else if(grepl('.loom', input)) {
      obj <- SeuratDisk::Connect(filename = input, mode = "r")
      obj <- Seurat::as.Seurat(obj)
   } else if(grepl('.h5ad', input)){
      SeuratDisk::Convert(input, dest = "h5seurat", overwrite = TRUE)
      h5seurat <- list.files(path = ".", pattern = ".h5seurat", recursive = TRUE, full.names = T)
      obj <- SeuratDisk::LoadH5Seurat(h5seurat)
   } else {
      h5_data <- Seurat::Read10X_h5(filename = input, use.names = TRUE, unique.features = TRUE)
      
      if(class(h5_data) == "list") {
         cts <- h5_data$'Gene Expression'
         obj <- Seurat::CreateSeuratObject(counts = cts, project = "Seurat", min.cells = 3, min.features = 200)
      } else {
         obj <- Seurat::CreateSeuratObject(counts = h5_data, project = "Seurat", min.cells = 3, min.features = 200)
      } 
   }
   saveRDS(obj, file = "seurat_object.rds")
}