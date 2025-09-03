setwd("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H1. harmony/GSE171213/")
library(tidyverse)
library(data.table)
library(Seurat)

##----Step 1. Extract fibroblasts and classify into OAT+ / OAT− ----
Idents(obj) <- obj$anntation
obj.fib <- subset(obj, idents = "Fibroblasts")   # subset fibroblast cluster
# Extract OAT gene expression
oat_expression <- FetchData(obj.fib, vars = "OAT")
# Add OAT expression as metadata
obj.fib <- AddMetaData(obj.fib, metadata = oat_expression, col.name = "OAT_Expression")
# Define fibroblasts as OAT+ if OAT expression > 0, otherwise OAT−
obj.fib$annotation2 <- ifelse(obj.fib$OAT_Expression > 0, "OAT+ Fib", "OAT- Fib")
# Save intermediate fibroblast object
saveRDS(obj.fib, "obj.fib1006.rds")
obj.fib <- readRDS("obj.fib1006.rds")


##----Step 2. Reclustering fibroblasts (harmony integration + UMAP)----
library(harmony)
library(clustree)
# Standard preprocessing
obj.fib <- NormalizeData(obj.fib, normalization.method = "LogNormalize", scale.factor = 10000)
obj.fib <- FindVariableFeatures(obj.fib, selection.method = "vst", nfeatures = 1500)
obj.fib <- ScaleData(obj.fib, vars.to.regress = c("percent.mt", "percent.rb"))
obj.fib <- RunPCA(obj.fib, npcs = 30, verbose = TRUE)
Seurat::ElbowPlot(obj.fib, ndims = 30)   # elbow plot to determine PCs
# Batch correction with Harmony
obj.fib <- RunHarmony(obj.fib, group.by.vars = "orig.ident")
# Dimensionality reduction and clustering
obj.fib <- RunUMAP(obj.fib, reduction = "harmony", dims = 1:15)
obj.fib <- FindNeighbors(obj.fib, reduction = "harmony", dims = 1:15)
obj.fib <- FindClusters(obj.fib, resolution = seq(from = 0.1, to = 1, by = 0.1))
# Visualize cluster tree across resolutions
clustree(obj.fib, prefix = "RNA_snn_res.") + coord_flip()
# Select clustering resolution 0.9
Idents(obj.fib) <- obj.fib$RNA_snn_res.0.9
table(obj.fib$diseaseStatus)
# Visualization: UMAP split by disease status
DimPlot(obj.fib, pt.size = 1, label = TRUE, split.by = "diseaseStatus")
# OAT expression visualization
FeaturePlot(obj.fib, features = "OAT")  # fibroblasts only
FeaturePlot(obj, features = "OAT")      # all cell types

##----Step 3. Identify marker genes and assist annotation----
marker_cosg <- cosg(
  obj.fib,
  groups = 'all',
  assay = 'RNA',
  slot = 'data',
  mu = 1,
  remove_lowly_expressed = TRUE,
  expressed_pct = 0.1,
  n_genes_user = 50
)
# Alternative: FindAllMarkers
obj.fib <- readRDS(".../obj.fib.scp1006.rds")
DimPlot(obj.fib, group.by = "cluster")
Idents(obj.fib) <- obj.fib$cluster
FAM <- FindAllMarkers(obj.fib,
                      min.pct = 0.25,
                      logfc.threshold = 0.25,
                      only.pos = TRUE)
# Save results
data.table::fwrite(FAM, "obj.fib.FAMresults.csv")
top100 <- FAM %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
top50  <- FAM %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
data.table::fwrite(top50, "obj.fib.FAM.top50.csv")
data.table::fwrite(top100, "obj.fib.FAM.top100.csv")

# Manually assign clusters into subtypes Fib1–Fib6
DotPlot(obj.fib,features = "OAT",group.by = "annotation3")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
obj.fib$annotation3 <- NA
obj.fib$annotation3 <- ifelse(obj.fib$RNA_snn_res.0.9 %in% c(0,9),"Fib1",obj.fib$annotation3)
obj.fib$annotation3 <- ifelse(obj.fib$RNA_snn_res.0.9 %in% c(3,6,8),"Fib2",obj.fib$annotation3)
obj.fib$annotation3 <- ifelse(obj.fib$RNA_snn_res.0.9 %in% c(1,2,7),"Fib3",obj.fib$annotation3)
obj.fib$annotation3 <- ifelse(obj.fib$RNA_snn_res.0.9 %in% c(5),"Fib4",obj.fib$annotation3)
obj.fib$annotation3 <- ifelse(obj.fib$RNA_snn_res.0.9 %in% c(4),"Fib5",obj.fib$annotation3)
obj.fib$annotation3 <- ifelse(obj.fib$RNA_snn_res.0.9 %in% c(10),"Fib6",obj.fib$annotation3)
DimPlot(obj.fib,pt.size = 1,label = T,group.by = "annotation3")
saveRDS(obj.fib,"D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H10. Final/Figure2/obj.fib.rds")
Idents(obj.fib)=obj.fib$annotation3

DotPlot(obj.fib,features = markers2,group.by = "annotation3")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


##----Step 4. Visulization with SCP and scplotter R packages----
Idents(obj.fib) <- obj.fib$diseaseStatus
obj.fib.HC <- subset(obj.fib,ident = "HC")
obj.fib.PD <- subset(obj.fib,ident = "PD")
obj.fib.HC$annotation3 <- factor(obj.fib.HC$annotation3,levels = c("Fib1","Fib2","Fib3","Fib4","Fib5","Fib6"))
obj.fib.PD$annotation3 <- factor(obj.fib.PD$annotation3,levels = c("Fib1","Fib2","Fib3","Fib4","Fib5","Fib6"))


# UMAP by fibroblast subtypes (Fib1–Fib6)
p2.1.1 <- scplotter::CellDimPlot(obj.fib.HC, group_by = "annotation3", theme = "theme_blank")
p2.1.2 <- scplotter::CellDimPlot(obj.fib.PD, group_by = "annotation3", theme = "theme_blank")

# Direct UMAP plots
scplotter::CellDimPlot(obj.fib.HC, group_by = "annotation3", reduction = "umap")
scplotter::CellDimPlot(obj.fib.PD, group_by = "annotation3", reduction = "umap")

# Barplot: fibroblast subtypes across disease status
obj.fib$cluster <- obj.fib$annotation3
p2.2 <- scplotter::CellStatPlot(obj.fib, ident = "annotation3",
                                group_by = "diseaseStatus", 
                                frac = "group", position = "stack")

# Clustering tree
p2.3 <- scplotter::ClustreePlot(obj.fib, prefix = "RNA_snn_res.")


# Add OAT expression to metadata
OAT_expression <- FetchData(obj.fib, vars = "OAT")
obj.fib$OAT <- OAT_expression

# Expression statistics
p2.4 <- scplotter::FeatureStatPlot(obj.fib, features = "OAT",
                                   ident = "diseaseStatus",
                                   facet_scales = "free_y", add_point = TRUE)

# SCP version with statistical comparison
p2.4.1 <- SCP::FeatureStatPlot(srt = obj.fib, slot = "data",
                               group.by = "diseaseStatus", stat.by = c("OAT"),
                               add_point = TRUE, comparisons = list(c("HC","PD")))


# OAT expression across fibroblast subtypes
p2.4.2 <- DotPlot(obj.fib, features = "OAT", group.by = "annotation3") +
  RotatedAxis() + ggtitle("Expression of OAT")

# Multiple marker genes
features <- c("CXCL6","CXCL1","CXCL5","CXCL13",
              "CD34","APOD","GPX3",
              "COL11A1","CILP2","OGN","THBS1",
              "MAFB","TMEM26",
              "LAMP5","CST1","CST2",
              "MKI67","TOP2A")

DotPlot(obj.fib, features = features, group.by = "annotation3") +
  RotatedAxis() + ggtitle("Annotation")

# Differential expression test
obj.fib.scp <- SCP::RunDEtest(srt = obj.fib, group_by = "annotation3", fc.threshold = 1)

# Volcano plot
SCP::VolcanoPlot(srt = obj.fib.scp, group_by = "annotation3")

# Extract DEGs
DEGs <- obj.fib.scp@tools$DEtest_annotation3$AllMarkers_wilcox
DEGs <- DEGs[with(DEGs, avg_log2FC > 1 & p_val_adj < 0.05), ]

# Annotate features (TFs + surface proteins)
obj.fib.scp <- SCP::AnnotateFeatures(obj.fib.scp, species = "Homo_sapiens", db = c("TF", "CSPA"))

# Feature heatmap with functional annotation
ht <- SCP::FeatureHeatmap(srt = obj.fib.scp,
                          group.by = "annotation3",
                          features = DEGs$gene,
                          db = c("GO_BP", "KEGG", "WikiPathway"),
                          anno_terms = TRUE)
print(ht$plot)
