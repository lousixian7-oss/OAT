##--------------------------------
##Cellchat HC vs PD
##--------------------------------
library(CellChat)
library(Seurat)
Idents(obj2) <- obj2$diseaseStatus #obj2为成纤维细胞亚群
obj.HC = subset(obj2,idents = "HC")
obj.PD = subset(obj2,idents = "PD")
Idents(obj.HC) <- obj.HC$anntation
Idents(obj.PD) <- obj.PD$anntation
cellchat.HC <- createCellChat(object = obj.HC@assays[["RNA"]]@data, meta = obj.HC@meta.data, group.by ="anntation")
cellchat.PD <- createCellChat(object = obj.PD@assays[["RNA"]]@data, meta = obj.PD@meta.data, group.by ="anntation")

## 针对HC的处理
cellchat=cellchat.HC
CellChatDB <- CellChatDB.human
cellchat@DB  <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat,do.fast = FALSE) ##识别过表达基因
cellchat <- identifyOverExpressedInteractions(cellchat) ##识别过表达受体配体对
cellchat <- projectData(cellchat, PPI.human) ##受体配体表达值投射到PPI
cellchat@idents <- droplevels(cellchat@idents)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE,population.size =T)##使用表达值推测细胞互作的概率
cellchat <- filterCommunication(cellchat, min.cells = 3)##过滤
cellchat <- computeCommunProbPathway(cellchat)##计算每个信号通路相关的所有配体-受体相互作用的通信结果
cellchat <- aggregateNet(cellchat)##汇总
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchat.HC = cellchat ## 短暂保留结果

## 针对PD的处理
cellchat=cellchat.PD
CellChatDB <- CellChatDB.human
cellchat@DB  <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat,do.fast = FALSE) ##识别过表达基因
cellchat <- identifyOverExpressedInteractions(cellchat) ##识别过表达受体配体对
cellchat <- projectData(cellchat, PPI.human) ##受体配体表达值投射到PPI
cellchat@idents <- droplevels(cellchat@idents)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE,population.size =T)##使用表达值推测细胞互作的概率
cellchat <- filterCommunication(cellchat, min.cells = 3)##过滤
cellchat <- computeCommunProbPathway(cellchat)##计算每个信号通路相关的所有配体-受体相互作用的通信结果
cellchat <- aggregateNet(cellchat)##汇总
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchat.PD = cellchat ## 短暂保留结果

saveRDS(cellchat.HC,"D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/cellchat.HC.rds")
saveRDS(cellchat.PD,"D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/cellchat.PD.rds")

cellchat.HC <- readRDS("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/cellchat.HC.rds")
cellchat.PD <- readRDS("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/cellchat.PD.rds")
cco.list <- list(HC = cellchat.HC, 
                 PD = cellchat.PD)
cellchat_merge <- mergeCellChat(cco.list, add.names = names(cco.list), cell.prefix = TRUE)


##--------------------------------
##可视化
##--------------------------------
gg1 <- compareInteractions(cellchat_merge, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat_merge, show.legend = F, group = c(1,2), measure = "weight")
p <- gg1 + gg2
p

par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat_merge, weight.scale = T)
netVisual_diffInteraction(cellchat_merge, weight.scale = T, measure = "weight")


par(mfrow = c(1,1))
h1 <- netVisual_heatmap(cellchat_merge)
h2 <- netVisual_heatmap(cellchat_merge, measure = "weight")
h1+h2

## 通路信号强度对比分析
gg1 <- rankNet(cellchat_merge, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat_merge, mode = "comparison", stacked = F, do.stat = TRUE)
p2 <- gg1 + gg2
p2

pathway.union <- union(cco.list[[1]]@netP$pathways, cco.list[[2]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(cco.list[[1]], pattern = "all", signaling = pathway.union, 
                                        title = names(cco.list)[1], width = 8, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(cco.list[[2]], pattern = "all", signaling = pathway.union,
                                        title = names(cco.list)[2], width = 8, height = 10)
ht1+ht2
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))

#cicrle图，主要关心PD中独有的通路
pathways.show <- c("WNT") 
par(mfrow = c(1,1), xpd=TRUE)
netVisual_aggregate(cellchat.PD, signaling = "WNT")
#和弦图
netVisual_aggregate(cellchat.PD, signaling ="WNT", layout = "chord")
netVisual_aggregate(cellchat.PD, signaling ="FASLG", layout = "chord")
netVisual_aggregate(cellchat.PD, signaling ="NRG", layout = "chord")
netVisual_aggregate(cellchat.PD, signaling ="CHEMERIN", layout = "chord") # ANNEXIN、CXCL和TNF
netVisual_aggregate(cellchat.PD, signaling ="ANNEXIN", layout = "chord")
netVisual_aggregate(cellchat.PD, signaling ="CXCL", layout = "chord")
netVisual_aggregate(cellchat.PD, signaling ="TNF", layout = "chord")
netVisual_aggregate(cellchat.HC, signaling ="SEMA3", layout = "chord")
## 信号通路点图展示
netAnalysis_signalingRole_scatter(cellchat.HC)
netAnalysis_signalingRole_scatter(cellchat.PD)

#####################
## 计算受体和配体得分
#####################
cellchat.PD.new <- netAnalysis_computeCentrality(cellchat.PD, slot.name = "netP")
par(mfrow = c(2,2))
a1 <- netAnalysis_signalingRole_network(cellchat.PD.new, signaling = "WNT", width = 8, height = 2.5, font.size = 10)
a2 <- netAnalysis_signalingRole_network(cellchat.PD.new, signaling = "FASLG", width = 8, height = 2.5, font.size = 10)
a3 <- netAnalysis_signalingRole_network(cellchat.PD.new, signaling = "NRG", width = 8, height = 2.5, font.size = 10)
a4 <- netAnalysis_signalingRole_network(cellchat.PD.new, signaling = "CHEMERIN", width = 8, height = 2.5, font.size = 10)
P <- a1+a2+a3+a4
##从所有信号通路对聚合的细胞-细胞通信网络的信号作用分析
gg1 <- netAnalysis_signalingRole_scatter(cellchat.PD.new)
###从所有信号通路对聚合的细胞-细胞通信网络的信号作用分析
gg2 <- netAnalysis_signalingRole_scatter(cellchat.PD.new, signaling = c("WNT", "FASLG","NRG","CHEMERIN"))
gg1 + gg2


##识别对某些细胞群的传出或传入信号贡献最大的信号，从所有信号通路对聚合的细胞-细胞通信网络的信号作用分析。
ht1 <- netAnalysis_signalingRole_heatmap(cellchat.PD.new, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat.PD.new, pattern = "incoming")
ht1 + ht2

##气泡图(主要展示具体的受体配体在细胞之间的作用)
levels(cellchat.PD.new@idents)
netVisual_bubble(cellchat.PD.new, sources.use = 3, targets.use = c(1:15), remove.isolate = FALSE)
##sources.use = 2 是值第二个细胞亚群
netVisual_bubble(cellchat.PD.new, sources.use =c(3,7), targets.use = c(1:5), remove.isolate = FALSE)
##大小表示P，互作强度看颜色
##指定信号通路
cellchat.PD.new@netP$pathways 
netVisual_bubble(cellchat.PD.new, sources.use =c(1:15), targets.use =c(3),signaling =  c("WNT", "FASLG","NRG","CHEMERIN"), remove.isolate = FALSE)

##用小提琴图绘制信号基因的表达分布 参与某条信号通路（如TGFb）的所有基因在细胞群中的表达情况展示
plotGeneExpression(cellchat.PD.new, signaling = "WNT")
plotGeneExpression(cellchat.PD.new, signaling = "FASLG")
plotGeneExpression(cellchat.PD.new, signaling = "NRG")
plotGeneExpression(cellchat.PD.new, signaling = "CHEMERIN")


##-----------------------------------
## OAT和受体配体相关性分析
##-----------------------------------
library(tidyverse)
expr <- obj2@assays$RNA
gene_name <- c("OAT","FZD4","LRP5","FZD6")
gene_expression <- expr %>% 
  .[gene_name,] %>% 
  t() %>% 
  as.data.frame()
colnames(gene_expression) <- paste0(gene_name)
identical(colnames(obj2),row.names(gene_expression))
obj2$OAT <- gene_expression[,paste0("OAT")]
obj2$FZD4 <- gene_expression[,paste0("FZD4")]
obj2$FZD6 <- gene_expression[,paste0("FZD6")]
obj2$LRP5 <- gene_expression[,paste0("LRP5")]

identical(obj2@meta.data[,paste0(gene_name)],gene_expression[,paste0(gene_name)])
meta <- obj2@meta.data

meta <- meta[,c("orig.ident","diseaseStatus","anntation","OAT","FZD4","FZD6","LRP5")]
meta1 <- meta[,c(2,3,4,5)]

cor.data <- meta1[which(meta1$anntation == "Fibroblasts"),]
cor.data$OAT = as.numeric(cor.data$OAT)
cor.data$FZD4 = as.numeric(cor.data$FZD4)
PD.cor<-cor.test(cor.data[which(cor.data$diseaseStatus=="PD"),"OAT"],cor.data[which(cor.data$diseaseStatus == "PD"), "FZD4"])
HC.cor<-cor.test(cor.data[which(cor.data$diseaseStatus=="HC"),"OAT"],cor.data[which(cor.data$diseaseStatus == "HC"), "FZD4"])


cor.data_filtered <- cor.data[cor.data$OAT > 0 & cor.data$FZD4 > 0, ]
cor.data_filtered$OAT = as.numeric(cor.data_filtered$OAT)
cor.data_filtered$FZD4 = as.numeric(cor.data_filtered$FZD4)
PD.cor<-cor.test(cor.data_filtered[which(cor.data_filtered$diseaseStatus=="PD"),"OAT"],cor.data_filtered[which(cor.data_filtered$diseaseStatus == "PD"), "FZD4"])
HC.cor<-cor.test(cor.data_filtered[which(cor.data_filtered$diseaseStatus=="HC"),"OAT"],cor.data_filtered[which(cor.data_filtered$diseaseStatus == "HC"), "FZD4"])

cor.data_filtered <- cor.data[cor.data$OAT > 0 & cor.data$FZD4 > 0, ]


p1 <- ggplot(cor.data_filtered, aes(OAT, FZD4)) +  geom_point(color = "#00AFBB")

p2 <- p1 +  
  geom_smooth(method="lm", se=T) +  
  geom_hline(yintercept = 1.5, linetype = "dashed", color = "blue") +  
  geom_vline(xintercept = 1.5, linetype = "dashed", color = "red") +  
  annotate("text", x=3, y=2.5, parse=TRUE,label="r^2 == 0.0138 * ' p-value = 0.1529' ")                    

ggMarginal(p2, type = "histogram", fill = "#00AFBB")

library(ggstatsplot)           
ggscatterstats(   
  data = cor.data_filtered,                                           
  x = Sepal.Length,                                                   
  y = Sepal.Width,  
  xlab = "Sepal Length",  
  ylab = "Sepal Width",  
  marginal = TRUE,      
  marginal.type = "densigram",  
  margins = "both",  
  xfill = "blue", # 分别设置颜色  
  yfill = "#009E73",  
  title = "Relationship between Sepal Length and Sepal Width",  
  messages = FALSE
)#这里我们还需要调用一下R包ggstatsplot，否则就会报错

