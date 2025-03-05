# ScTypeWorkflow

## Description
The `ScTypeWorkflow` function automates cell type annotation for single-cell RNA sequencing (scRNA-seq) data using the ScType workflow (https://sctype.app/). The function takes in a Seurat object and returns cell type annotations for each cluster in the dataset. It also provides an option to export the resulting cell type scores to a CSV file. 

The gene sets for various tissues can be found here: https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx

## Arguments
- `seurat`: The Seurat object containing your single-cell RNA sequencing data.
- `tissue`: The tissue type for which to generate gene sets. Available tissue types include Immune system, Pancreas, Liver, Eye, Kidney, Brain, Lung, Adrenal, Heart, Intestine, Muscle, Placenta, Spleen, Stomach, and Thymus (default is "Immune system").
- `resolution`: The cluster resolution to use for assigning annotations (default is "seurat_clusters").
- `export`: If TRUE, exports the resulting cell type scores to a CSV file (default is FALSE).
- `export_path`: Specifies the path for saving the exported CSV file (if export = TRUE).

## Citations
- Ianevski A, Giri AK, Aittokallio T (2022). "Fully-automated and ultra-fast cell-type identification using specific marker combinations from single-cell transcriptomic data." _Nature Communications_, **13**, 1246. https://doi.org/10.1038/s41467-022-28803-w.
