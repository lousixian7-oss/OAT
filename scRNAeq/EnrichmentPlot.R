#------------------------
## SCP 差异分析 GO分析
##------------------------
obj.fib.scp <- SCP::RunDEtest(srt = obj.fib, 
                              group_by = "annotation3", 
                              fc.threshold = 1, 
                              only.pos = FALSE)
SCP::VolcanoPlot(srt = obj.fib.scp,group_by = "annotation3")

DEGs <- obj.fib.scp@tools$DEtest_annotation3$AllMarkers_wilcox
DEGs <- DEGs[with(DEGs, avg_log2FC > 1 & p_val_adj < 0.05), ]
# Annotate features with transcription factors and surface proteins
obj.fib.scp <- SCP::AnnotateFeatures(obj.fib.scp, 
                                     species = "Homo_sapiens", 
                                     db = c("TF", "CSPA"))
ht <- SCP::FeatureHeatmap(
  srt = obj.fib.scp, 
  group.by = "annotation3", 
  features = DEGs$gene, 
  feature_split = DEGs$group1,
  species = "Homo_sapiens", 
  db = c("GO_BP", "KEGG", "WikiPathway"), 
  anno_terms = TRUE,
  feature_annotation = c("TF", "CSPA"), 
  feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  height = 5, width = 4
)
print(ht$plot)


obj.fib.scp <- SCP::RunEnrichment(
  srt = obj.fib.scp, 
  group_by = "annotation3", 
  db = "GO_BP", 
  species = "Homo_sapiens",
  DE_threshold = "avg_log2FC > log2(1.5) & p_val_adj < 0.05"
)

SCP::EnrichmentPlot(
  srt = obj.fib.scp, 
  group_by = "annotation3", 
  group_use = c("Fib5","Fib6"),
  plot_type = "bar"
)

obj.fib.scp <- SCP::RunGSEA(
  srt = obj.fib.scp, 
  group_by = "annotation3", 
  db = "GO_BP", 
  species = "Homo_sapiens",
  DE_threshold = "p_val_adj < 0.05"
)
SCP::GSEAPlot(srt = obj.fib.scp, 
              group_by = "annotation3", 
              group_use = "Fib1", 
              id_use = "GO:0007186")

SCP::GSEAPlot(
  srt = obj.fib.scp, 
  group_by = "annotation3", 
  group_use = "Fib1", 
  plot_type = "bar",
  direction = "both", 
  topTerm = 20
)

saveRDS(obj.fib.scp,"D://scRNAseq/06.Proj_processing/Sen_yazhouyan/H10. Final/Figure2/obj.fib.scp1006.rds")




##-----------------
## 补充分析
## 针对obj.fib 的 OAT+ 和 OAT- 的成纤维细胞的DEG 
##---------------------------------------------
table(obj.fib$annotation2)
# OAT- Fib OAT+ Fib 
# 294      712 

Idents(obj.fib) <- obj.fib$annotation2
degs <- FindAllMarkers(obj.fib,
                       #min.pct = 0.25,
                       #logfc.threshold = 0.25,
                       test.use = "wilcox",
                       only.pos = F)
library(CellChat)
library(EnhancedVolcano)
library(ggrepel)
degs <- degs[-1,]

EnhancedVolcano(degs,                
                lab = row.names(degs),                
                x = 'avg_log2FC',                
                y = 'p_val_adj')

?EnhancedVolcano()
EnhancedVolcano(degs,    
                lab = row.names(degs),    
                x = 'avg_log2FC',    
                y = 'p_val_adj',    
                title = 'OAT+ versus OAT-',
                pCutoff = 0.1,###修改p截断值    
                FCcutoff = 1,# 修改log2FC截断值    
                pointSize = 3.0,##点大小    
                labSize = 6.0###标签大小
)

all.markers = markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

top_gene = degs %>% 
  group_by(cluster) %>% 
  top_n(n = 5, wt = avg_log2FC)


#logfc.threshold = 0.25
#test.use = "wilcox"
#min.pct = 0.1

markers_2 <- FindMarkers(obj.fib,ident.1="OAT+ Fib",ident.2="OAT- Fib",assay = 'RNA',slot = 'data')
markers_2 <- markers_2[-1,]
EnhancedVolcano(markers_2,    
                lab = row.names(markers_2),    
                x = 'avg_log2FC',    
                y = 'p_val_adj',    
                title = 'OAT+ versus OAT-',
                pCutoff = 0.05,###修改p截断值    
                FCcutoff = 0.3,# 修改log2FC截断值    
                pointSize = 3.0,##点大小    
                labSize = 3.0,###标签大小
                xlim = c(-1, 1),
                ylim = c(0,8)
)

# 筛选 avg_log2FC > 0.25 且 p_val_adj < 0.05 的基因
key_genes <- markers_2 %>% 
  filter(avg_log2FC > 0.25, p_val_adj < 0.05)

# 查看筛选后的结果
print(key_genes)
key_genes$Genesymbol <- row.names(key_genes)
gene.df <- bitr(key_genes$Genesymbol, ## 需要转化的IDD，你也可以单独提取gene名作为向量
                fromType = "SYMBOL", #fromType是指你的数据ID类型是属于哪一类的
                toType = c("ENTREZID", "SYMBOL"), #toType是指你要转换成哪种ID类型，可以写多种，也可以只写一种
                OrgDb = org.Hs.eg.db) #Orgdb是指对应的注释包是哪个
colnames(gene.df)[1] <- "Genesymbol"
DEG_data1 <- left_join(gene.df,key_genes) ## 完善关于DEG的信息


GO_all <- enrichGO(gene = DEG_data1$ENTREZID,  #基因列表(转换的ID)
                   keyType = "ENTREZID",  #指定的基因ID类型，默认为ENTREZID
                   OrgDb=org.Hs.eg.db,  #物种对应的org包
                   ont = "ALL",   #CC细胞组件，MF分子功能，BF生物学过程，ALL以上三个
                   pvalueCutoff = 0.05,  #p值阈值
                   pAdjustMethod = "fdr",  #多重假设检验校正方式，所以针对的是padjust！！！
                   minGSSize = 10,   #注释的最小基因集，默认为10
                   maxGSSize = 500,  #注释的最大基因集，默认为500
                   qvalueCutoff = 0.05,  #q值阈值
                   readable = TRUE)  #基因ID转换为基因名
GO_result <- data.frame(GO_all)  


kegg <- enrichKEGG(DEG_data1$ENTREZID,                    
                   organism = 'hsa',                    
                   pvalueCutoff = 0.05,                    
                   pAdjustMethod = 'BH',                    
                   minGSSize = 10,                    
                   maxGSSize = 500,                    
                   qvalueCutoff = 0.05,                    
                   use_internal_data = FALSE)
kegg_result <- data.frame(kegg)
write.csv(kegg@result, file = './RawData/Enrichanalysis/KEGG.result.csv')#把基因ID转为基因名
KEGG <- setReadable(kegg, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")


p5 <- barplot(GO_all,showCategory = 20) 
go_enrichment_pathway <- GO_result %>% group_by(ONTOLOGY) %>% top_n(n = 12, wt = -p.adjust)
p7 <- ggplot(go_enrichment_pathway, aes(x=reorder(Description, -Count), y=Count, fill=ONTOLOGY)) +
  geom_bar(stat="identity", position=position_dodge(width=0.8), width=0.6) +
  theme_minimal() +
  labs(x="GO Term", y="Gene_Number", title="Top 10 Enriched GO Terms") +
  facet_grid(ONTOLOGY~., scale = 'free_y', space = 'free_y')+
  coord_flip() +  #让柱状图变为纵向
  scale_fill_manual(values=c("CC"="skyblue","BP"="pink","MF"="lightgreen"))+
  theme_bw()
p7

p1 <- barplot(kegg, drop = TRUE, showCategory = 15,color = "p.adjust",title = "KEGG Pathway")
p1

