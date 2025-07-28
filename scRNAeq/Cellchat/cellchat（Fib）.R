##---------------------------------------
## 更新cellchat，OAT+ OAT- Fib加入分析
## 先在全局的metadata里annotation里写清楚 OAT+ OAT-
##--------------------------------------
table(obj.fib@meta.data$annotation2)
# OAT- Fib OAT+ Fib 
# 294      712
sample <- row.names(obj.fib@meta.data)
OAT1 <- row.names(obj.fib@meta.data)[obj.fib$annotation2 == "OAT- Fib"]
OAT2 <- row.names(obj.fib@meta.data)[obj.fib$annotation2 == "OAT+ Fib"]

obj$annotation_new <- obj$anntation
obj$annotation_new <- as.character(obj$annotation_new)
obj$names <- row.names(obj@meta.data)
obj$annotation_new[obj$names %in% OAT1] <- "OAT- Fib"
obj$annotation_new[obj$names %in% OAT2] <- "OAT+ Fib"
obj$annotation_new <- as.factor(obj$annotation_new)
DimPlot(obj,group="annotation_new")

library(CellChat)
library(Seurat)
Idents(obj) <- obj$diseaseStatus
obj.HC = subset(obj,idents = "HC")
obj.PD = subset(obj,idents = "PD")
Idents(obj.HC) <- obj.HC$anntation
Idents(obj.PD) <- obj.PD$anntation
cellchat.HC <- createCellChat(object = obj.HC@assays[["RNA"]]@data, meta = obj.HC@meta.data, group.by ="annotation_new")
cellchat.PD <- createCellChat(object = obj.PD@assays[["RNA"]]@data, meta = obj.PD@meta.data, group.by ="annotation_new")

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
cellchat.HC = cellchat ## 短暂保留结

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
cellchat.PD = cellchat ## 短暂保留结

saveRDS(cellchat.HC,"cellchat.HC.rds")
saveRDS(cellchat.PD,"cellchat.PD.rds")

cellchat.HC <- readRDS("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/cellchat.HC.rds")
cellchat.PD <- readRDS("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/cellchat.PD.rds")
cco.list <- list(HC = cellchat.HC, 
                 PD = cellchat.PD)
cellchat_merge <- mergeCellChat(cco.list, add.names = names(cco.list), cell.prefix = TRUE)


## 【第二步】可视化
gg1 <- compareInteractions(cellchat_merge, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat_merge, show.legend = F, group = c(1,2), measure = "weight")
p <- gg1 + gg2
p

par(mfrow = c(1,1))
h1 <- netVisual_heatmap(cellchat_merge)
h2 <- netVisual_heatmap(cellchat_merge, measure = "weight")
h1+h2

## 通路信号强度对比分析
gg1 <- rankNet(cellchat_merge, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat_merge, mode = "comparison", stacked = F, do.stat = TRUE)
p2 <- gg1 + gg2
p2

##气泡图(主要展示具体的受体配体在细胞之间的作用)
levels(cellchat.PD.new@idents)
netVisual_bubble(cellchat.PD.new, sources.use = 11, targets.use = c(4,1,16,7,8), remove.isolate = FALSE)
##sources.use = 2 是值第二个细胞亚群
netVisual_bubble(cellchat.PD.new, sources.use =c(3,7), targets.use = c(1:5), remove.isolate = FALSE)
##大小表示P，互作强度看颜色
##指定信号通路
cellchat.PD.new@netP$pathways 
netVisual_bubble(cellchat.PD.new, sources.use =c(1:16), targets.use =c(11),signaling =  c("WNT", "FASLG","NRG","CHEMERIN"), remove.isolate = FALSE)

##用小提琴图绘制信号基因的表达分布 参与某条信号通路（如TGFb）的所有基因在细胞群中的表达情况展示
# plotGeneExpression(cellchat.PD.new, signaling = "WNT")
plotGeneExpression(cellchat.PD.new, signaling = "FASLG")
plotGeneExpression(cellchat.PD.new, signaling = "NRG")
plotGeneExpression(cellchat.PD.new, signaling = "CHEMERIN")