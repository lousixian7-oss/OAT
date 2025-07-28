setwd("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H1. harmony/GSE171213/")
library(tidyverse)
library(data.table)
library(Seurat)

obj <- readRDS("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H1. harmony/GSE171213/annotation_finished.rds")
DimPlot(obj,group.by = "anntation",label = T,split.by = "diseaseStatus")
dim(obj)
# [1] 24540 18966

Idents(obj) <- obj$anntation
obj.fib <- subset(obj,idents = "Fibroblasts") 
oat_expression <- FetchData(obj.fib, vars = "OAT")
obj.fib <- AddMetaData(object = obj.fib, metadata = oat_expression, col.name = "OAT_Expression")
obj.fib$annotation2 <- ifelse(obj.fib$OAT_Expression > 0,"OAT+ Fib","OAT- Fib")
saveRDS(object = obj.fib,"obj.fib1006.rds")
obj.fib <- readRDS("obj.fib1006.rds")

##-----------------------
## 成纤维细胞继续分群探究
##----------------------——
library(harmony)
library(clustree)
obj.fib <- NormalizeData(obj.fib, normalization.method = "LogNormalize", scale.factor = 10000)
obj.fib <- FindVariableFeatures(obj.fib, selection.method = "vst",nfeatures = 1500) 
obj.fib <- ScaleData(obj.fib, vars.to.regress = c("percent.mt","percent.rb")) 
obj.fib <- RunPCA(obj.fib, npcs = 30,verbose=T) 
Seurat::ElbowPlot(obj.fib, ndims = 30) 
obj.fib <- RunHarmony(obj.fib, group.by.vars = "orig.ident")
obj.fib <- RunUMAP(obj.fib, reduction = "harmony", dims = 1:15)
obj.fib <- FindNeighbors(obj.fib, reduction = "harmony", dims = 1:15)
obj.fib <- FindClusters(obj.fib, resolution = seq(from = 0.1, to = 1, by = 0.1)) 
clustree(obj.fib, prefix = "RNA_snn_res.")+coord_flip() 

Idents(obj.fib)=obj.fib$RNA_snn_res.0.9 
table(obj.fib$diseaseStatus)
DimPlot(obj.fib,pt.size = 1,label = T,split.by = "diseaseStatus") 
FeaturePlot(obj.fib,features = "OAT")
Seurat::FeaturePlot(obj,features = "OAT")

library(COSG)
table(Idents(obj.fib)) 
marker_cosg <- cosg(
  obj.fib,
  groups='all', 
  assay='RNA',
  slot='data',
  mu=1,         
  remove_lowly_expressed=TRUE,   
  expressed_pct=0.1,            
  n_genes_user=50      
)

## Findallmarkers找到差异gene 取top100 进行富集
## GO/KEGG富集之后，去定义每个细胞群的功能
obj.fib <- readRDS("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H10. Final/Figure2/obj.fib.scp1006.rds")
DimPlot(obj.fib,group.by = "cluster")
Idents(obj.fib) <- obj.fib$cluster
FAM <- FindAllMarkers(obj.fib,
                      min.pct = 0.25,
                      logfc.threshold = 0.25,
                      only.pos = T)
getwd()
data.table::fwrite(FAM,"obj.fib.FAMresults.csv")
data.table::fwrite(top50,"obj.fib.FAM.top50.csv")
data.table::fwrite(top100,"obj.fib.FAM.top100.csv")
top100 <-  FAM %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
top50 <-  FAM  %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)


markers.fib <- c(
  "DCN",
  "ACTA2","ASPN","COL1A1",
  "LAMP5",
  "AHRR","PDE10A","MAFB","TGFB2",
  "APOE","LPL","CCL19",
  "CXCL6","CXCL1", 
  "MKI67","TOP2A",
  "MMP11","MMP13","IGSF10","COL11A1"
)

markers2 <- c(
  "CXCL6","CXCL1","CXCL5","CXCL13",
  "CD34","APOD","GPX3",
  "COL11A1","CILP2","OGN","THBS1",
  "MAFB","TMEM26",
  "LAMP5","CST1","CST2",
  "MKI67","TOP2A",
  "OAT"
)

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

DimPlot(obj.fib,pt.size = 1,label = T,group.by = "annotation2")
dittoBarPlot(obj.fib, "annotation2", group.by = "annotation3",main = '')
dittoBarPlot(obj.fib, "annotation3", group.by = "annotation2",main = '')
dittoBarPlot(obj.fib, "annotation3", group.by = "diseaseStatus",main = '')


# 绘制OAT基因表达量比较图片
library(ggplot2)
library(ggpubr)
p <- ggplot(obj.fib@meta.data, aes(x = diseaseStatus, y = OAT_Expression, fill = diseaseStatus)) +
  geom_violin(trim = FALSE, color = "black") + 
  geom_boxplot(width = 0.2, color = "white", fill = NA) +  
  geom_jitter(width = 0.15, size = 1, color = "black") +  
  theme_classic() +  # 使用简洁主题
  theme(
    axis.text = element_text(size = 12, color = 'black'),
    axis.line = element_line(size = 1),
    axis.title.y = element_text(size = 15, color = 'black'),
    plot.title = element_text(size = 15, hjust = 0.5, color = 'black')
  ) +
  labs(
    title = "OAT Expression by Disease Status",
    x = "Disease Status",
    y = "Expression"
  ) +
  scale_fill_manual(values = c("#F8766D", "#00BA38")) +
  ylim(0,4)
# 添加显著性比较
p <- p + stat_compare_means(
  aes(group = diseaseStatus), 
  method = "wilcox.test",  
  label = "p.signif"  
)

print(p)

