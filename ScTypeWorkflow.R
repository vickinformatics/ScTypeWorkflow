### Description: The ScTypeWorkflow() function is designed to automate cell type annotation on single-cell RNA-sequencing data using the ScType framework (https://sctype.app/). The function takes in a Seurat object and returns cell type annotations for each cluster in the dataset. Optionally, it can also export the resulting scores to a CSV file.
### The list of tissues and corresponding gene sets can be found here: https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx

### Original code: https://github.com/IanevskiAleksandr/sc-type

### Last updated: 2025-2-26

### Author: Vicki Do

ScTypeWorkflow <- function(seurat,
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
    write.csv(sctype_scores, file = export_path, row.names = FALSE)
    print(paste("sctype_scores table has been exported to", export_path))
  }
  
  return(sctype_scores)
}