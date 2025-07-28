##---------------------------------------------------
## SCENIC R部分
##--------------------------------------------------
library(Matrix)
library(SCENIC)
seurat_to_adata <- function(object,#seurat对象
                            Dimension=c('UMAP','TSNE'),#降维方式
                            path){#文件保存路径
  seurat_obj <- object
  seurat_obj$barcode <- colnames(seurat_obj)
  if(Dimension=='UMAP'){
    cell.embeddings<- seurat_obj@reductions$umap@cell.embeddings
    seurat_obj$UMAP_1 <- cell.embeddings[,1]
    seurat_obj$UMAP_2 <- cell.embeddings[,2]
  }else{
    
    cell.embeddings<- seurat_obj@reductions$tsne@cell.embeddings
    seurat_obj$TSNE_1 <- cell.embeddings[,1]
    seurat_obj$TSNE_2 <- cell.embeddings[,2]
  }
  #保存metadat
  write.csv(seurat_obj@meta.data, file=paste0(path,'metadata.csv'), quote=F, row.names=F)
  #保存matrix
  counts_matrix <- GetAssayData(seurat_obj, assay='RNA', slot='counts')
  writeMM(counts_matrix, file=paste0(path, 'counts.mtx'))
  #PCA降维信息
  write.csv(seurat_obj@reductions$pca@cell.embeddings, file=paste0(path,'pca.csv'), quote=F,row.names=F)
  
  #保存gene name
  write.table(data.frame('gene'=rownames(counts_matrix)),file=paste0(path,'gene_names.csv'),
              quote=F,row.names=F,col.names=F)
}

adata <- seurat_to_adata(obj.fib,Dimension='UMAP',path = 'D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H9. SCENIC/')


## 保存矩阵
write.csv(t(as.matrix(obj.fib@assays$RNA@counts)),file = "sce_exp.csv")

## 得到sce_SCENIC.loom文件之后开始进行下游分析
sce_SCENIC <- open_loom(
  "D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H9.1. SCENIC/pythonProject/sce_SCENIC.loom"
)
exprMat <- get_dgem(sce_SCENIC)#从sce_SCENIC文件提取表达矩阵
exprMat_log <- log2(exprMat+1) # log处理
regulons_incidMat <- get_regulons(sce_SCENIC, column.attr.name="Regulons")
regulons <- regulonsToGeneLists(regulons_incidMat)
class(regulons)

regulonAUC <- get_regulons_AUC(sce_SCENIC, column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(sce_SCENIC)

## 可视化操作
obj.fib <- readRDS("d:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H1. harmony/GSE171213/obj.fib1006.rds")
obj.fib <- readRDS("d:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H10. Final/Figure2/obj.fib.scp1006.rds")
human_data <- obj.fib ## 原先的seurat对象
cellinfo <- human_data@meta.data[,c('annotation2','diseaseStatus',"nFeature_RNA","nCount_RNA","annotation3")]#细胞meta信息
colnames(cellinfo)=c('celltype', 'group','nGene' ,'nUMI','annotation3')


######计算细胞特异性TF(这里实际计算的是OAT+ 和 OAT- 成纤维细胞的区别)
cellTypes <-  as.data.frame(subset(cellinfo,select = 'annotation2'))
selectedResolution <- "annotation2"
sub_regulonAUC <- regulonAUC

rss <- calcRSS(AUC=getAUC(sub_regulonAUC),
               cellAnnotation=cellTypes[colnames(sub_regulonAUC),
                                        selectedResolution])
rss=na.omit(rss)
rss1 <- rss[c("IKZF2(+)","TFE3(+)","STAT1(+)","KLF6(+)","EGR1(+)","FOS(+)","PURA(+)","NFYA(+)","ZBED1(+)","STAT6(+)"),]
rssPlot <- 
  plotRSS(
    rss1, ## rss全局，主要参数Z值和 RSS值
    zThreshold = 0.1,
    cluster_columns = FALSE,
    order_rows = TRUE,
    thr=0.1,
    varName = "cellType",
    col.low = '#330066',
    col.mid = '#66CC66',
    col.high = '#FFCC33')
rssPlot
dim(rss)
summary(rss)

## 提取转录因子打分到seurat对象
next_regulonAUC <- regulonAUC[,match(colnames(human_data),colnames(regulonAUC))]
dim(next_regulonAUC)
library(SummarizedExperiment)
regulon_AUC <- regulonAUC@NAMES
human_data@meta.data = cbind(human_data@meta.data ,t(assay(next_regulonAUC[regulon_AUC,])))

#选定重要的转录因子
TF_plot <- c("KLF6(+)","EGR1(+)","FOS(+)")

DotPlot(human_data, features = TF_plot)+
  theme_bw()+
  theme(panel.grid = element_blank(), 
        axis.text.x=element_text(hjust =1,vjust=1, angle = 45))+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=3))

FeaturePlot(human_data, features ="KLF6(+)")
FeaturePlot(human_data, features ="EGR1(+)")
FeaturePlot(human_data, features ="FOS(+)")
FeaturePlot(human_data, features ="FOSB(+)")
DimPlot(human_data,group.by = "anntation")


sce_regulons <- read.csv("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H9.1. SCENIC//pythonProject/sce.regulons.csv")
sce_regulons <- sce_regulons[-2,]
colnames(sce_regulons) <- sce_regulons[1,]
sce_regulons <- sce_regulons[-1,]
colnames(sce_regulons) <- c("TF","ID","AUC","NES",
                            "MotifsimilarityQvalue",
                            "OrthologousIdentity",
                            "Annotation",
                            "Context",
                            "TargetGenes","RankAtMax")



FOS <- subset(sce_regulons,TF == "FOS")
FOS <- FOS[which(FOS$AUC > 0.05),]
FOS <- FOS[,c("TF","TargetGenes")]
FOS$TargetGenes <- gsub("\\[","",FOS$TargetGenes)

library(stringr)
split_FOS <- str_split(FOS$TargetGenes,",")
FOS1 <- as.data.frame(split_FOS[[1]])
FOS2 <- as.data.frame(split_FOS[[2]])
FOS3 <- as.data.frame(split_FOS[[3]])
FOS4 <- as.data.frame(split_FOS[[4]])
FOS5 <- as.data.frame(split_FOS[[5]])

names(FOS1) <- "TF"
names(FOS2) <- "TF"
names(FOS3) <- "TF"
names(FOS4) <- "TF"
names(FOS5) <- "TF"


FOS.new <- rbind(FOS1,FOS2,FOS3,FOS4,FOS5)
FOS_target <- FOS.new[seq(1,nrow(FOS.new),2),]
FOS_score <- FOS.new[seq(0,nrow(FOS.new),2),]
FOS_gene <- data.frame(FOS_target,FOS_score)
FOS_gene <- FOS_gene[!duplicated(FOS_gene$FOS_target),]
FOS_gene$gene <- "FOS"
colnames(FOS_gene) <- c("target","score","tf")

# 假设你的数据框名为 df
FOS_gene <- FOS_gene %>%
  mutate(
    # 去除第一列的 ( 和 '，并重命名为 "gene"
    gene = str_replace_all(target, "[\\(')]", ""),
    # 去除第二列的最后一个 )
    score = str_replace(score, "\\)$", "")
  ) %>%
  # 保留所需列并移除原来的 target 列
  select(gene, score, tf)

TF_target <- rbind(FOS_gene,FOSB_gene)
TF_target <- na.omit(TF_target)
TF_target$score <- as.numeric(TF_target$score)

# 绘制网络图
path <- c("FOS","KLF6","EGR1")
nodelist <- list()
for (i in 1:length(path)){
  node <- subset(TF_target,tf == path[i]) # 提取数据
  nodes <- data.frame(name = unique(union(node$tf,node$gene))) # 整理为data
  nodes$value <- c(sum(node$score)/10,node$score) # 加上values
  nodelist[[i]] <- nodes
}

nodes <- rbind(nodelist[[1]],nodelist[[2]]) # 合并节点文件
nodes$cluster <- c(rep("FOS",1),rep("FOS_gene",462),
                   rep("FOSB",1),rep("FOSB_gene",279))

edges <- TF_target[c("tf","gene","score")]
edges$class <- edges$tf
sum(is.na(edges))      
edges <- na.omit(edges)


library(ggraph)
library(tidygraph)
colnames(edges)
colnames(nodes) <- c("gene","score","cluster")
layout_cir <- tbl_graph(nodes = nodes,edges = edges)
nodes$cluster <- gsub("_gene$", "",nodes$cluster)
ggraph(layout_cir,layout = "linear",cirular = T)+
  geom_node_point(aes(size=value,colour = cluster)) +
  geom_node_text(aes(x = 1.03 * x,
                     y = 1.03 * y,
                     label = name,
                     color =cluster,
                     angle = -((-node_angle(x,y)+90) %% 180)+90),
                 hjust = 'outward') + 
  geom_edge_arc(aes(colour = class)) +
  theme_void()+
  theme(legend.position = "none")+
  scale_colour_manual(values = c('#407972',
                                 '#961E28',
                                 '#D46724',
                                 '#0f8096'))+
  scale_edge_colour_manual(values = c('#961E28',
                                      '#D46724',
                                      '#0f8096'))+
  scale_size_continuous(range = c(2,8))+
  coord_cartesian(xlim = c(-1.5,-1.5),ylim = c(-1.5,-1.5))

## 创建自己的节点和边数据
edges1 <- edges %>% filter(.,score>1.74)
edges1 <- edges1 %>%
  mutate(gene = str_trim(gene)) %>%  # 去除前后空格
  distinct(gene, .keep_all = TRUE)    # 删除重复的基因行

nodes1 <- edges1 %>%
  select(tf, gene) %>%
  pivot_longer(cols = c(tf, gene), names_to = "type", values_to = "name") %>%
  distinct(name) %>%
  mutate(cluster = if_else(name %in% edges$tf, "tf", "gene"))

# 创建 tbl_graph 对象
graph <- tbl_graph(nodes = nodes1, edges = edges1)
# 使用 ggraph 绘制图形
ggraph(graph, layout = "linear", circular = TRUE) +
  geom_edge_arc(aes(width = score, color = class), show.legend = FALSE) +
  geom_node_point(aes(color = cluster, size = if_else(cluster == "TF", 5, 3))) +
  geom_node_text(aes(label = name, color = cluster), repel = TRUE) +
  scale_color_manual(values = c("TF" = "#00A08A", "Gene" = "#F98400")) +
  theme_void() +
  theme(legend.position = "none")

# 使用 graph_from_data_frame 函数创建图对象
graph <- graph_from_data_frame(TF_final, directed = TRUE)
# 设置节点类型 (转录因子或靶基因)
V(graph)$type <- ifelse(V(graph)$name %in% TF_final$tf, "TF", "Target")
# 自动生成颜色：为每个转录因子分配一个唯一颜色
tf_colors <- setNames(colorRampPalette(c("#00A08A", "#F98400", "#5BBCD6", "#E79F00", "#F2AD00"))(length(unique(TF_final$tf))), unique(TF_final$tf))
# 绘制图形
ggraph(graph, layout = "linear", circular = TRUE) +
  geom_edge_arc(aes(width = score, color = tf), alpha = 0.6) +  # 边的宽度和颜色映射到 tf
  geom_node_point(aes(color = type), size = 5) +  # 使用不同颜色区分 TF 和 Target
  geom_node_text(aes(label = name), repel = TRUE, size = 3) +  # 在节点旁标注基因名称
  scale_color_manual(values = c("TF" = "#00A08A", "Target" = "#F98400")) +  # 设置节点颜色
  scale_edge_color_manual(values = tf_colors) +  # 自动为转录因子分配颜色
  theme_void() +
  theme(legend.position = "none") +
  labs(title = "Transcription Factor Regulatory Network", subtitle = "Based on SCENIC Analysis")