# write-seurat-object
This package is for the creation of a Seurat object from single-cell RNA sequences files with the following file extensions: “.h5”, “.gz”, “.h5ad”, and “.loom”. 

## Installation
The package can be installed using 
```
devtools::install_github(“BTIP2024/write-seurat-object”)
```


## Example
The output of this function would be an rds file which when loaded in R, would be a Seurat object. 
```
write_seurat_object(“scRNAdata.h5”)
```

## scRNAseq processing workflow 

The standard scRNAseq processing workflow with the R package Seurat consists of seven (7) steps. The output of this package and function should be used as input for the scRNAseq processing pipeline. 

The following are the repositories of the packages for every step of the pipeline:
1. QC and filtering: [qualitycontrolseurat package](https://github.com/BTIP2024/quality-control-seurat)
2. Normalization: [qualitycontrolseurat package](https://github.com/BTIP2024/quality-control-seurat)
3. Identification of highly variable features: [selectionscalingseurat package](https://github.com/BTIP2024/selection-scaling-seurat)
4. Scaling: [selectionscalingseurat package](https://github.com/BTIP2024/selection-scaling-seurat)
5. Linear Dimensionality Reduction (PCA): [pcaseurat package](https://github.com/BTIP2024/pca-seurat)
6. Clustering: [nonlinearreduction package](https://github.com/BTIP2024/non-linear-reduction)
7. Non-linear dimensionality reduction (t-SNE and UMAP): [nonlinearreduction package](https://github.com/BTIP2024/non-linear-reduction)

An overview of the pipeline and its outputs can be observed below:
![](https://github.com/user-attachments/assets/eff6edef-5d4b-4b16-a7f7-f5aebc5b6b17)
