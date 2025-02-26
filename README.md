# ScTypeWorkflow

## Description
The ScTypeWorkflow() function is designed to automate cell type annotation on single-cell RNA-sequencing data using the ScType framework (https://sctype.app/). The function takes in a Seurat object and returns cell type annotations for each cluster in the dataset. Optionally, it can also export the resulting scores to a CSV file.

The list of tissues and corresponding gene sets can be found here: https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx

## Arguments
- **seurat**: The Seurat object containing your single-cell RNA sequencing data.
- **tissue**: The tissue type for which to generate gene sets. Available tissue types include Immune system, Pancreas, Liver, Eye, Kidney, Brain, Lung, Adrenal, Heart, Intestine, Muscle, Placenta, Spleen, Stomach, and Thymus (default is "Immune system").
- **resolution**: The cluster resolution to use for assigning annotations (default is "seurat_clusters").
- **export**: If TRUE, exports the resulting cell type scores to a CSV file.
- **export_path**: Specifies the path for saving the exported CSV file (if export = TRUE).
