#' Run_ScType
#'
#' The `Run_ScType()` function automates cell type annotation for single-cell RNA sequencing (scRNA-seq) data using the ScType workflow (https://sctype.app/).
#' The function takes in a Seurat object and returns cell type annotations for each cluster in the dataset.
#' It also provides an option to export the resulting cell type scores to a CSV file. 
#' The gene sets for various tissues can be found here: https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx
#'
#' @param seurat A Seurat object containing the single-cell RNA sequencing data.
#' @param tissue The tissue type for which to generate gene sets. Available tissue types include Immune system, Pancreas, Liver, Eye, Kidney, Brain, Lung, Adrenal, Heart, Intestine, Muscle, Placenta, Spleen, Stomach, and Thymus (default is "Immune system").
#' @param resolution The cluster resolution to use for assigning annotations (default is "seurat_clusters").
#' @param export If TRUE, exports the resulting cell type scores to a CSV file (default is FALSE)
#' @param export_path Specifies the path for saving the exported CSV file (if export = TRUE).
#' 
#' @return A data frame containing the top 10 annotated cell types for each cluster, along with their corresponding scores and the number of cells per cluster.
#'
#' @references
#' Ianevski A, Giri AK, Aittokallio T (2022). "Fully-automated and ultra-fast cell-type identification using specific marker combinations from single-cell transcriptomic data." *Nature Communications*, 13, 1246. https://doi.org/10.1038/s41467-022-28803-w.
#' 
#' @author Vicki Do
#' @lastUpdated 2025-11-7
#'
#' @examples
#' sctype_scores <- Run_ScType(seurat = seurat_obj, tissue = "Kidney", resolution = "RNA_snn_res.0.5", export = TRUE, export_path = "tables/Kidney_scType_scores.csv")

Run_ScType() <- function(seurat,
                         tissue = "Immune system",
                         resolution = "seurat_clusters",
                         export = FALSE,
                         export_path = "sctype_scores.csv") {
  
  # Load gene set preparation and cell type annotation functions
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
  source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
  
  # DB file
  db_ <- "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx"
  
  # Prepare gene sets
  gs_list <- gene_sets_prepare(db_, tissue)
  
  # Join layers
  seurat_joined <- JoinLayers(seurat)
  
  # Check Seurat object version (scRNA-seq matrix extracted differently in Seurat v4/v5)
  seurat_package_v5 <- isFALSE('counts' %in% names(attributes(seurat_joined[["RNA"]])))
  print(sprintf("Seurat object %s is used", ifelse(seurat_package_v5, "v5", "v4")))
  
  # Extract scaled counts
  scaled_counts <- if (seurat_package_v5) as.matrix(seurat_joined[["RNA"]]$scale.data) else as.matrix(seurat_joined[["RNA"]]@scale.data)
  
  # Run ScType
  es.max <- sctype_score(scRNAseqData = scaled_counts, scaled = TRUE, gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)
  
  # Merge by cluster
  sctype_scores <- do.call("rbind", lapply(unique(seurat_joined@meta.data[[resolution]]), function(cl) {
    es.max.cl = sort(rowSums(es.max[, rownames(seurat_joined@meta.data[seurat_joined@meta.data[[resolution]] == cl, ])]), decreasing = TRUE)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_joined@meta.data[[resolution]] == cl)), 10)
  }))
  
  # Order by cluster
  sctype_scores$cluster <- factor(sctype_scores$cluster, levels = sort(as.numeric(levels(sctype_scores$cluster))))
  sctype_scores <- sctype_scores[order(sctype_scores$cluster), ]
  
  # Set low-confident (low ScType score) clusters to "unknown"
  sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells / 4] <- "Unknown"
  
  # If export is TRUE, save the sctype_scores as a CSV file
  if (export) {
    # Create directory if it doesn't exist
    dir.create(dirname(export_path), recursive = TRUE, showWarnings = FALSE)
    
    write.csv(sctype_scores, file = export_path, row.names = FALSE)
    print(paste("sctype_scores table has been exported to", export_path))
  }
  
  return(sctype_scores)
}
