##--------------------------------------------------------------------
## GSEA+UCell
##--------------------------------------------------------------------
library(Seurat)
library(GSVA) 
library(GSEABase)
library(msigdbr)
library(clusterProfiler)
# Load Seurat object
obj <- readRDS("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H10. Final/obj18966.rds")
# Check cluster annotations
table(obj$annotation)
Idents(obj) <- obj$anntation

##----Step 1. Average expression per cluster----
av <- AverageExpression(obj,assays="RNA")
av=av[[1]]
cg=names(tail(sort(apply(av,1,sd)),1000)) # Select top 1000 most variable genes across clusters
pheatmap::pheatmap(cor(av[cg,])) # Heatmap of correlations across clusters (based on variable genes)


##----Step 2. Prepare gene sets (MSigDB Reactome pathways)----
gs <- msigdbr_collections()
all_gene_sets = msigdbr(species = "Homo sapiens",
                        category = "C2", 
                        subcategory = "CP:REACTOME")


##----Step 3. Run GSEA for each cluster----
geneList = av[,1] 
gl = apply(av, 2, function(geneList){
  geneList=sort(geneList,decreasing = T)
  print(head(geneList))
  print(tail(geneList))
  # Filter: keep genes with average expression > 0.1
  geneList=geneList[geneList>0.1]
  egmt <- GSEA(geneList, TERM2GENE= all_gene_sets[,c('gs_name','gene_symbol')] , 
               minGSSize = 20, 
               pvalueCutoff = 1,
               verbose=FALSE)
  head(egmt)
  egmt@result 
  gsea_results_df <- egmt@result 
  return(gsea_results_df)  # return GSEA result table
})

##————Step 4. Collect enrichment scores into matrix----
path = unique(all_gene_sets$gs_name)

es.max <- do.call(cbind,
                  lapply(gl, function(x){
                    x[path,'enrichmentScore']
                  }))
rownames(es.max) = path
es.max = na.omit(es.max)

# Heatmap of enrichment scores
pheatmap::pheatmap(es.max, show_colnames = TRUE, show_rownames = FALSE)


##----Step 5. Identify top 5 cluster-specific pathways----
library(dplyr)
df = do.call(rbind,
             lapply(1:ncol(es.max), function(i){
               data.frame(
                 path  = rownames(es.max),
                 cluster = colnames(es.max)[i],
                 sd.1 = es.max[,i],
                 sd.2 = apply(es.max[,-i], 1, median)  # median of other clusters
               )
             }))
# Define fold change as enrichment difference
df$fc = df$sd.1 - df$sd.2
# Select top 5 enriched pathways per cluster
top5 <- df %>% group_by(cluster) %>% top_n(5, fc)
n=es.max[unique(top5$path),]
rownames(n)
rownames(n)[grepl('FGFR',rownames(n))]
rownames(n)=gsub('REACTOME_','',rownames(n))
rownames(n)=substring(rownames(n),1,30)
pheatmap::pheatmap(n,show_rownames = T)


##----Step 6. Validation with UCell scoring----
library(UCell)
library(msigdbr)
obj <- readRDS("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H10. Final/obj18966.rds")
# Prepare pathway gene sets
geneSets <- msigdbr(species="Homo sapiens", category="C2", subcategory="CP:REACTOME") %>%
  split(x=.$gene_symbol, f=.$gs_name)

pathway1 <- geneSets[["REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION"]]
pathway2 <- geneSets[["REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX"]]
pathway3 <- geneSets[["REACTOME_SIGNALING_BY_MET"]]

# Run UCell to score pathways at single-cell level
DefaultAssay(obj) <- 'RNA'
obj <- AddModuleScore_UCell(obj, features = pathway1, name="_score")
obj <- AddModuleScore_UCell(obj, features = pathway2, name="_score")
obj <- AddModuleScore_UCell(obj, features = pathway3, name="_score")

# Visualize pathway activity on UMAP
P1.7.1 <- FeaturePlot(obj, features = "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION_score",
                      order = TRUE, cols = viridis::viridis(256))
P1.7.2 <- FeaturePlot(obj, features = "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX_score",
                      order = TRUE, cols = viridis::viridis(256))
P1.7.3 <- FeaturePlot(obj, features = "REACTOME_SIGNALING_BY_MET_score",
                      order = TRUE, cols = viridis::viridis(256))
