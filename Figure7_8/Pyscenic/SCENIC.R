##---------------------------------------------------
## SCENIC analysis in R (transcription factor activity inference)
##---------------------------------------------------

library(Matrix)
library(SCENIC)

#--------------------------------------
# Function: Convert Seurat object to SCENIC-compatible input files
#--------------------------------------
seurat_to_adata <- function(object,                    # Input Seurat object
                            Dimension=c('UMAP','TSNE'),# Dimension reduction method
                            path){                     # Output directory
  seurat_obj <- object
  seurat_obj$barcode <- colnames(seurat_obj)
  
  # Extract UMAP or tSNE coordinates
  if(Dimension=='UMAP'){
    cell.embeddings <- seurat_obj@reductions$umap@cell.embeddings
    seurat_obj$UMAP_1 <- cell.embeddings[,1]
    seurat_obj$UMAP_2 <- cell.embeddings[,2]
  } else {
    cell.embeddings <- seurat_obj@reductions$tsne@cell.embeddings
    seurat_obj$TSNE_1 <- cell.embeddings[,1]
    seurat_obj$TSNE_2 <- cell.embeddings[,2]
  }
  
  # Save metadata
  write.csv(seurat_obj@meta.data, file=paste0(path,'metadata.csv'), quote=F, row.names=F)
  # Save count matrix (MTX format)
  counts_matrix <- GetAssayData(seurat_obj, assay='RNA', slot='counts')
  writeMM(counts_matrix, file=paste0(path, 'counts.mtx'))
  # Save PCA embeddings
  write.csv(seurat_obj@reductions$pca@cell.embeddings, file=paste0(path,'pca.csv'), quote=F,row.names=F)
  # Save gene names
  write.table(data.frame('gene'=rownames(counts_matrix)),
              file=paste0(path,'gene_names.csv'),
              quote=F,row.names=F,col.names=F)
}

# Example: convert fibroblast Seurat object
adata <- seurat_to_adata(obj.fib, Dimension='UMAP',
                         path = 'D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H9. SCENIC/')

# Also save full counts as CSV
write.csv(t(as.matrix(obj.fib@assays$RNA@counts)), file = "sce_exp.csv")

#--------------------------------------
# Step 1: Load SCENIC .loom results
#--------------------------------------
sce_SCENIC <- open_loom("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H9.1. SCENIC/pythonProject/sce_SCENIC.loom")

# Extract expression matrix and regulons
exprMat <- get_dgem(sce_SCENIC)
exprMat_log <- log2(exprMat+1)
regulons_incidMat <- get_regulons(sce_SCENIC, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)

# Extract regulon AUC activity scores
regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(sce_SCENIC)

#--------------------------------------
# Step 2: Prepare metadata for downstream visualization
#--------------------------------------
obj.fib <- readRDS("d:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H10. Final/Figure2/obj.fib.scp1006.rds")
human_data <- obj.fib
cellinfo <- human_data@meta.data[,c('annotation2','diseaseStatus',"nFeature_RNA","nCount_RNA","annotation3")]
colnames(cellinfo) <- c('celltype', 'group','nGene' ,'nUMI','annotation3')

#--------------------------------------
# Step 3: Calculate TF specificity (OAT+ vs OAT− fibroblasts)
#--------------------------------------
cellTypes <- as.data.frame(subset(cellinfo,select = 'annotation2'))
selectedResolution <- "annotation2"
sub_regulonAUC <- regulonAUC

rss <- calcRSS(
  AUC=getAUC(sub_regulonAUC),
  cellAnnotation=cellTypes[colnames(sub_regulonAUC), selectedResolution]
)
rss <- na.omit(rss)

# Subset example TFs for plotting
rss1 <- rss[c("IKZF2(+)","TFE3(+)","STAT1(+)","KLF6(+)","EGR1(+)","FOS(+)","PURA(+)","NFYA(+)","ZBED1(+)","STAT6(+)"),]

rssPlot <- plotRSS(
  rss1,
  zThreshold = 0.1,
  cluster_columns = FALSE,
  order_rows = TRUE,
  thr=0.1,
  varName = "cellType",
  col.low = '#330066',
  col.mid = '#66CC66',
  col.high = '#FFCC33'
)

#--------------------------------------
# Step 4: Add TF activity scores to Seurat metadata
#--------------------------------------
next_regulonAUC <- regulonAUC[,match(colnames(human_data),colnames(regulonAUC))]
library(SummarizedExperiment)
regulon_AUC <- regulonAUC@NAMES
human_data@meta.data <- cbind(human_data@meta.data ,t(assay(next_regulonAUC[regulon_AUC,])))

# Visualization of selected TFs
TF_plot <- c("KLF6(+)","EGR1(+)","FOS(+)")
DotPlot(human_data, features = TF_plot)+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(hjust =1,vjust=1, angle = 45))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))

FeaturePlot(human_data, features="KLF6(+)")
FeaturePlot(human_data, features="EGR1(+)")
FeaturePlot(human_data, features="FOS(+)")
FeaturePlot(human_data, features="FOSB(+)")

#--------------------------------------
# Step 5: Parse TF-target gene interactions (example: FOS)
#--------------------------------------
sce_regulons <- read.csv("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H9.1. SCENIC/pythonProject/sce.regulons.csv")
sce_regulons <- sce_regulons[-2,]
colnames(sce_regulons) <- sce_regulons[1,]
sce_regulons <- sce_regulons[-1,]
colnames(sce_regulons) <- c("TF","ID","AUC","NES",
                            "MotifsimilarityQvalue",
                            "OrthologousIdentity",
                            "Annotation",
                            "Context",
                            "TargetGenes","RankAtMax")

# Extract FOS targets with AUC > 0.05
FOS <- subset(sce_regulons, TF == "FOS")
FOS <- FOS[FOS$AUC > 0.05, c("TF","TargetGenes")]
FOS$TargetGenes <- gsub("\\[","",FOS$TargetGenes)

# Parse comma-separated target genes into data frames
library(stringr)
split_FOS <- str_split(FOS$TargetGenes,",")
FOS1 <- as.data.frame(split_FOS[[1]])
...
# (Repeat for FOS2–FOS5, combine, clean duplicates, etc.)

#--------------------------------------
# Step 6: Build TF regulatory network
#--------------------------------------
TF_target <- rbind(FOS_gene,FOSB_gene)
TF_target <- na.omit(TF_target)
TF_target$score <- as.numeric(TF_target$score)

# Build node list (example for FOS, KLF6, EGR1)
path <- c("FOS","KLF6","EGR1")
nodelist <- list()
for (i in 1:length(path)){
  node <- subset(TF_target, tf == path[i])
  nodes <- data.frame(name = unique(union(node$tf,node$gene)))
  nodes$value <- c(sum(node$score)/10,node$score)
  nodelist[[i]] <- nodes
}

# Build graph and visualize with ggraph
library(ggraph)
library(tidygraph)

edges <- TF_target[c("tf","gene","score")]
edges$class <- edges$tf
edges <- na.omit(edges)

layout_cir <- tbl_graph(nodes = nodes, edges = edges)
ggraph(layout_cir, layout = "linear", circular = TRUE)+
  geom_node_point(aes(size=value,colour=cluster)) +
  geom_node_text(aes(label=name, color=cluster), repel=TRUE) + 
  geom_edge_arc(aes(colour=class)) +
  theme_void()
