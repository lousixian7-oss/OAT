###############################################################################
# Monocle2 + CytoTRACE Pseudotime Analysis
# Purpose: Use CytoTRACE to assist Monocle2 trajectory inference and 
#          determine differentiation root state
# Data: obj.fib (Seurat fibroblast subset)
###############################################################################
library(Seurat)
library(monocle)
library(Matrix)

##----Step 1. Convert Seurat object to Monocle CDS----
obj.fib <- readRDS("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/obj.fib.rds")
sce_fib_new_monocle_fun <- seurat_to_monocle(obj.fib, assay = "RNA", slot = "counts")
# Normalize and estimate dispersion
sce_fib_new_monocle <- estimateSizeFactors(sce_fib_new_monocle_fun)
sce_fib_new_monocle <- estimateDispersions(sce_fib_new_monocle)
# QC: keep genes expressed in ≥10 cells (min_expr = 0.1)
sce_fib_new_monocle <- detectGenes(sce_fib_new_monocle, min_expr = 0.1)
expressed_genes <- row.names(subset(fData(sce_fib_new_monocle), num_cells_expressed >= 10))
# Add UMI counts
pData(sce_fib_new_monocle)$Total_mRNAs <- Matrix::colSums(exprs(sce_fib_new_monocle))
sce_fib_new_monocle <- sce_fib_new_monocle[, pData(sce_fib_new_monocle)$Total_mRNAs < 1e6]

##----Step 2. Select ordering genes & reduce dimension----
cds_DGT <- sce_fib_new_monocle
diff_test_res <- differentialGeneTest(cds_DGT, fullModelFormulaStr = "~annotation3")
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
cds_DGT <- setOrderingFilter(cds_DGT, ordering_genes)
plot_ordering_genes(cds_DGT)

# Dimension reduction & cell ordering
cds_DGT <- reduceDimension(cds_DGT, max_components = 2, reduction_method = 'DDRTree')
cds_DGT <- orderCells(cds_DGT, reverse = TRUE)

saveRDS(cds_DGT, "D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H10. Final/Figure2/cds_DGT_monocle.rds")

##----Step 3. Visualization of trajectories----
#初步可视化
plot_cell_trajectory(cds_DGT, color_by = "annotation3")
plot_cell_trajectory(cds_DGT, color_by = "Pseudotime")
plot_cell_trajectory(cds_DGT, color_by = "diseaseStatus")
plot_cell_trajectory(cds_DGT, markers = "OAT", use_color_gradient = TRUE)

# Save PDF plots
pdf(".../Pserdotime_umap.pdf", width = 6, height = 6)
plot_cell_trajectory(cds_DGT, color_by = "Pseudotime")
dev.off()

pdf(".../Pserdotime_OAT_time.pdf", width = 5, height = 4)
plot_genes_in_pseudotime(cds_DGT["OAT", ], color_by = "annotation3")
dev.off()


##----Step 4. CytoTRACE scoring----
library(CytoTRACE)
mat <- as.matrix(obj.fib@assays$RNA@counts)
results <- CytoTRACE(mat = mat)   # Run CytoTRACE
obj.fib$CytoTRACE <- results$CytoTRACE

# Map CytoTRACE score into Monocle trajectory
pData(cds_DGT)$cytotrace <- results$CytoTRACE[rownames(pData(cds_DGT))]

# Determine root state = state with highest mean CytoTRACE
state_scores <- aggregate(pData(cds_DGT)$cytotrace, by = list(State = pData(cds_DGT)$State), mean)
root_state <- state_scores$State[which.max(state_scores$x)]
cds_DGT <- orderCells(cds_DGT, root_state = root_state)

# Visualize CytoTRACE-colored trajectory
library(RColorBrewer)
cols <- rev(brewer.pal(11, "Spectral"))
plot_cell_trajectory(cds_DGT, color_by = "cytotrace", size = 1, show_backbone = TRUE) +
  scale_colour_gradientn(colours = cols)

pdf(".../cytotrace_Spectral.pdf", width = 6, height = 6)
plot_cell_trajectory(cds_DGT, color_by = "cytotrace", size = 1, show_backbone = TRUE) +
  scale_colour_gradientn(colours = cols)
dev.off()


##----Step 5. UMAP visualization in Seurat----
library(SCP)
NM <- SCP::FeatureDimPlot(obj.fib, features = "CytoTRACE", reduction = "umap")
pdf(".../cytotrace_umap.pdf", width = 5, height = 5)
SCP::FeatureDimPlot(obj.fib, features = "CytoTRACE", reduction = "umap")
dev.off()

## Function
seurat_to_monocle <- function(otherCDS, assay, slot, lowerDetectionLimit = 0, import_all = FALSE) {
  if(class(otherCDS)[1] == 'Seurat') {
    requireNamespace("Seurat")
    data <- GetAssayData(otherCDS, assay = assay, slot = slot)
    data <- data[rowSums(as.matrix(data)) != 0,]
    pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)
    if(all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- uninormal()
    } else {
      expressionFamily <- tobit()
    }
    valid_data <- data[, row.names(pd)]
    monocle_cds <- newCellDataSet(data,
                                  phenoData = pd, 
                                  featureData = fd,
                                  lowerDetectionLimit=lowerDetectionLimit,
                                  expressionFamily=expressionFamily)
    if(import_all) {
      if("Monocle" %in% names(otherCDS@misc)) {
        otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
        otherCDS@misc$Monocle@auxClusteringData$scran <- NULL
        monocle_cds <- otherCDS@misc$Monocle
        mist_list <- otherCDS
      } else {
        mist_list <- otherCDS
      }
    } else {
      mist_list <- list()
    } 
  } 
  return(monocle_cds)
}
