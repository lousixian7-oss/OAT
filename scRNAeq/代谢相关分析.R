##-----------------------------
## 针对irGSEA做代谢相关内容
## 使用代谢相关数据集
##-----------------------------
library(UCell)
library(irGSEA)
library(msigdbr)
library(clusterProfiler)
gmt_kegg <- read.gmt("D://R/R-4.3.2/library/scMetabolism/data/KEGG_metabolism_nc.gmt") %>% as.data.frame()
gmt_kegg$gene <- str_to_title(gmt_kegg$gene)
geneSets_kegg <- split(gmt_kegg$gene, gmt_kegg$term)

gmt_reactome <- read.gmt("D://R/R-4.3.2/library/scMetabolism/data/REACTOME_metabolism.gmt") %>% as.data.frame()
gmt_reactome$gene <- str_to_title(gmt_reactome$gene)
geneSets_reactome <- split(gmt_reactome$gene, gmt_reactome$term)

obj.fib.irGSEA <- irGSEA.score(object = obj.fib, 
                               assay = "RNA", 
                               slot = "data", 
                               seeds = 123, 
                               ncores = 1,
                               min.cells = 3,
                               min.feature = 0,
                               custom = T, 
                               geneset = geneSets_kegg, 
                               #msigdb = F, 
                               #species = "Homo sapiens", category = "H",  
                               #subcategory = NULL, 
                               #geneid = "symbol",
                               method = c("AUCell", "UCell", "singscore", "ssgsea"), # ssGSEA很慢
                               aucell.MaxRank = NULL, 
                               ucell.MaxRank = NULL, 
                               kcdf = 'Gaussian')

obj.fib.irGSEA <- irGSEA.score(object = obj.fib, 
                               assay = "RNA", 
                               slot = "data", 
                               seeds = 123, 
                               ncores = 1,
                               min.cells = 3,
                               min.feature = 0,
                               custom = T, 
                               geneset = geneSets_reactome, 
                               #msigdb = F, 
                               #species = "Homo sapiens", category = "H",  
                               #subcategory = NULL, 
                               #geneid = "symbol",
                               method = c("AUCell", "UCell", "singscore", "ssgsea"), # ssGSEA很慢
                               aucell.MaxRank = NULL, 
                               ucell.MaxRank = NULL, 
                               kcdf = 'Gaussian')
result.dge <- irGSEA.integrate(object = obj.fib.irGSEA, 
                               group.by = "annotation2",
                               metadata = NULL, col.name = NULL,
                               method = c("AUCell", "UCell", "singscore", "ssgsea"))
saveRDS(result.dge,file = "irGSEA.result.deg.kegg.rds")
getwd()
irGSEA.heatmap.plot <- irGSEA.heatmap(object = result.dge, 
                                      method = "RRA",
                                      top = 50, 
                                      show.geneset = NULL)
irGSEA.bubble.plot <- irGSEA.bubble(object = result.dge, 
                                    method = "RRA", 
                                    top = 50)
irGSEA.density.scatterplot <- irGSEA.density.scatterplot(object = obj.fib.irGSEA,
                                                         method = "UCell",
                                                         show.geneset = "Glycolysis / Gluconeogenesis",
                                                         reduction = "umap")

png("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H10. Final/Figure2/irGSEA.kegg.pathway3.png", 
    width = 5, height = 5, units = "in", res = 600)
print(irGSEA.density.scatterplot) 
dev.off()
pdf("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H10. Final/Figure2/irGSEA.kegg.pathway3.pdf", 
    width = 5, height = 5)
print(irGSEA.density.scatterplot) 
dev.off()


# 绘制相关性分析
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggExtra)
library(Seurat)
Idents(obj.fib.irGSEA) <- obj.fib.irGSEA$annotation2
obj.fib.OATpos <- subset(obj.fib.irGSEA,idents = "OAT+ Fib")
genes_of_interest <- c("OAT","DPP8")
pathway_of_interest <- c("Arginine and proline metabolism",
                         "Glycolysis / Gluconeogenesis",
                         "N-Glycan biosynthesis")
pathway_of_interest <- c("Glutamate and glutamine metabolism",
                         "Metabolism of amino acids and derivatives",
                         "Metabolism of nitric oxide NOS3 activation and regulation",
                         "Glycogen metabolism")
expression_matrix <- GetAssayData(object = obj.fib.OATpos, slot = "data")
extracted_matrix <- expression_matrix[genes_of_interest, ]
extracted_df <- as.data.frame(t(extracted_matrix))
pw_score_matrix <- obj.fib.OATpos@assays[["UCell"]]@data
pw_score_matrix <- pw_score_matrix[pathway_of_interest, ]
extracted_pw <- as.data.frame(t(pw_score_matrix))
extracted_df$RowNames <- rownames(extracted_df)
extracted_pw$RowNames <- rownames(extracted_pw)
merge_pw_df <- merge(extracted_df, 
                     extracted_pw, 
                     by = "RowNames", 
                     all = TRUE)


p <- 
  ggplot(merge_pw_df, aes_string(x = 'OAT',y = '`N-Glycan biosynthesis`')) +
  geom_point(size = 2,color = '#C1E0DB',alpha = 0.5) +
  theme_bw() +
  # 主题细节调整
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25,'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()
  ) +
  # 添加回归线
  geom_smooth(method = 'lm',se = T,color = '#F9B208',size = 1.5,fill = '#FEA82F') +
  # 添加相关性系数及p值
  stat_cor(method = "pearson",digits = 3,size=6)
# 添加边际柱形密度图
p2 <- ggMarginal(p,type = "densigram",
                 xparams = list(binwidth = 0.1, fill = "#B3E283",size = .7),
                 yparams = list(binwidth = 0.1, fill = "#8AB6D6",size = .7))

png("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H10. Final/Figure2/irGSEA.reactome.correlation4.png", 
    width = 6, height = 5, units = "in", res = 600)
print(p) 
dev.off()
pdf("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H10. Final/Figure2/irGSEA.reactome.correlation1.pdf", 
    width = 8, height = 8)
print() 
dev.off()
