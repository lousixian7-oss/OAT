##--------------------------------------------------------------------
## 完整的GSEA二次分析，针对obj所有细胞群体 2024-10-28
## GSEA是否适合单细胞分析存在争议，推荐的哪怕GSVA都不是首选，大概都是AUCell或者UCell打分
## 解决办法：一方面通过GSEA找到通路 然后用UCell验证；一方面直接从GSVA下手
##--------------------------------------------------------------------
library(Seurat)
library(GSVA) 
library(GSEABase)
library(msigdbr)
library(clusterProfiler)
obj <- readRDS("obj18966.rds")
table(obj$annotation_new)
Idents(obj) <- obj$anntation

av <- AverageExpression(obj,
                        assays="RNA")
av=av[[1]]
cg=names(tail(sort(apply(av,1,sd)),1000))
pheatmap::pheatmap(cor(av[cg,]))

gs <- msigdbr_collections()
all_gene_sets = msigdbr(species = "Homo sapiens",
                        category='H')   ## 采用Hallmarker的数据
all_gene_sets = msigdbr(species = "Homo sapiens",
                        category = "C2", subcategory =   "CP:REACTOME"  ) 
geneList = av[,1] 
gl = apply(av, 2, function(geneList){
  geneList=sort(geneList,decreasing = T)
  print(head(geneList))
  print(tail(geneList))
  geneList=geneList[geneList>0.1]
  egmt <- GSEA(geneList, TERM2GENE= all_gene_sets[,c('gs_name','gene_symbol')] , 
               minGSSize = 20, 
               pvalueCutoff = 1,
               verbose=FALSE)
  head(egmt)
  egmt@result 
  gsea_results_df <- egmt@result 
  return(gsea_results_df) 
  
})

path = unique(all_gene_sets$gs_name)

es.max <- do.call(cbind,
                  lapply(gl, function(x){
                    x[path,'enrichmentScore']
                  }))

rownames(es.max) = path
head(es.max)  
es.max=na.omit(es.max)
pheatmap::pheatmap(es.max,show_colnames =T,show_rownames = F) 
#每个单细胞亚群的特异性top5基因集的 富集分析结果
library(dplyr)
df = do.call(rbind,
             lapply(1:ncol(es.max), function(i){
               dat= data.frame(
                 path  = rownames(es.max),
                 cluster =   colnames(es.max)[i],
                 sd.1 = es.max[,i], #每一列原值
                 sd.2 = apply(es.max[,-i], 1, median)  #除当列以外每行的中位值
               )
             })) 
df$fc = df$sd.1 - df$sd.2#两值相减，变化越大说明越有意义（从中挑出top5）
top5 <- df %>% group_by(cluster) %>% top_n(5, fc)#找出每个细胞类型的前五个通路
n=es.max[unique(top5$path),]
rownames(n)
rownames(n)[grepl('FGFR',rownames(n))]
rownames(n)=gsub('REACTOME_','',rownames(n))
rownames(n)=substring(rownames(n),1,30)
pheatmap::pheatmap(n,show_rownames = T)


av.fib = av[,c(7,8)]
gl.fib = apply(av.fib, 2, function(geneList){
  geneList=sort(geneList,decreasing = T)
  print(head(geneList))
  print(tail(geneList))
  geneList=geneList[geneList>0.1]
  egmt <- GSEA(geneList, TERM2GENE= all_gene_sets[,c('gs_name','gene_symbol')] , 
               minGSSize = 20, 
               pvalueCutoff = 1,
               verbose=FALSE)
  head(egmt)
  return(egmt) 
  
})
gsea.fib <- gl.fib[["Fibroblasts"]]

library(GseaVis)
gseaNb(object = gsea.fib,
       geneSetID = 'REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION',
       subPlot = 2)

gseaNb(object = gsea.fib,
       geneSetID = 'REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX',
       subPlot = 2)
gseaNb(object = gsea.fib,
       geneSetID = 'REACTOME_SIGNALING_BY_MET',
       subPlot = 2)

gseaNb(object = gsea.fib,
       geneSetID = c('REACTOME_SIGNALING_BY_MET','REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX','REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION'),
       subPlot = 2)




library(msigdbr)
all_gene_sets = msigdbr(species = "Homo sapiens",
                        category = "C2", 
                        subcategory =  "CP:REACTOME")
get_geneSets.msigdbr <- function(order){
  ## Program: 获得注释数据集
  i <- order
  msigdbr_collections <- msigdbr_collections()
  m_df<- msigdbr(species = "Homo sapiens",   # 注意选择种类msigdbr::msigdbr_show_species()
                 category = msigdbr_collections$gs_cat[i], 
                 subcategory = msigdbr_collections$gs_subcat[i])
  geneSets<- m_df %>% 
    split(x = .$gene_symbol, f = .$gs_name) 
  return(geneSets)
}
geneSets <- get_geneSets.msigdbr(order = 7)
pathway1 <- subset(geneSets,names(geneSets) == "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION")
pathway2 <- subset(geneSets,names(geneSets) == "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX")
pathway3 <- subset(geneSets,names(geneSets) == "REACTOME_SIGNALING_BY_MET")
library(UCell)
DefaultAssay(obj) <- 'RNA'
obj <- AddModuleScore_UCell(obj,
                            features= pathway1,
                            name="_score")
obj <- AddModuleScore_UCell(obj,
                            features= pathway2,
                            name="_score")
obj <- AddModuleScore_UCell(obj,
                            features= pathway3,
                            name="_score")
P1.7.1 <- FeaturePlot(obj,features = "REACTOME_EXTRACELLULAR_MATRIX_ORGANIZATION_score",
                      order = T,cols = viridis::viridis(256))
P1.7.2 <- FeaturePlot(obj,features = "REACTOME_DEGRADATION_OF_THE_EXTRACELLULAR_MATRIX_score",
                      order = T,cols = viridis::viridis(256))
P1.7.3 <- FeaturePlot(obj,features = "REACTOME_SIGNALING_BY_MET_score",
                      order = T,cols = viridis::viridis(256))
