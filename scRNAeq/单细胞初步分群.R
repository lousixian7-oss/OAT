scRNAlist = readRDS("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H1. harmony/GSE171213/scRNAlist.rds")

wd = getwd()
objNames <- sapply(scRNAlist, function(x) x@project.name)
objs <- lapply(seq_along(scRNAlist), function(i){
  objName <- objNames[i]
  obj <- scRNAlist[[i]]
  scRNAdataPreProcessing(
    obj, 
    objName, 
    plotDir = wd, 
    minFeatures = 200, 
    maxFeatures = 5000, 
    minCounts = 1000, 
    maxCounts = 25000, 
    maxPctMito = 20, 
    nfeatures = 2000, 
    dims = 1:20, 
    #res = 0.4, 
    runDoubletFinder = T, 
    runDecontX = T, 
    estDubRate = 0.02, 
    ncores=1,
    use_logfile = F
  )
})

names(objs)=objNames
# Merge individual Seurat objects
obj <- merge(x=objs[[1]], 
             y=objs[2:length(objs)], 
             add.cell.ids=objNames, 
             project="scalp")
# Remove doublets
nPreDub <- dim(obj)[2]
obj <- subset(obj, subset = (DF.classify == "Singlet"))
nPostDub <- dim(obj)[2]
message(sprintf("Removed %s suspected doublets from %s total cells...", nPreDub - nPostDub, nPreDub))
useDecontX <- T
# If contaminated RNA removed, some cells may have very little RNA now
if(useDecontX){
  nPreDecon <- dim(obj)[2]
  obj <- subset(obj, subset = (nFeature_RNA > 200 & nCount_RNA > 500))
  nPostDecon <- dim(obj)[2]
  message(sprintf("Removed %s cells with too little RNA after decontamination. %s cells remaining.", 
                  nPreDecon - nPostDecon, nPostDecon))
}

# Save intermediate Seurat object:
message("Pre-processing complete. Saving intermediate Seurat object...")
saveRDS(obj, file = paste0(wd, "/preprocessed.rds"))


library(harmony)
library(clustree)
obj = readRDS("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H1. harmony/GSE171213/preprocessed.rds")
obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
obj <- FindVariableFeatures(obj, selection.method = "vst",nfeatures = 2000) 
obj <- ScaleData(obj, vars.to.regress = c("percent.mt","percent.rb")) 
obj <- RunPCA(obj, npcs = 50,verbose=T) 
Seurat::ElbowPlot(obj, ndims = 50) 
obj <- RunHarmony(obj, group.by.vars = "orig.ident")
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:20)
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:20)
obj <- FindClusters(obj, resolution = seq(from = 0.1, to = 1.5, by = 0.1)) 
saveRDS(obj, file = paste0(wd, "/preprocessed.rds"))
clustree(obj, prefix = "RNA_snn_res.")+coord_flip() 
Idents(obj)=obj$RNA_snn_res.0.7 
obj$diseaseStatus <- "HC"
obj$diseaseStatus <- ifelse(grepl("PD", obj$orig.ident), "PD", obj$diseaseStatus)
table(obj$diseaseStatus)

library(dittoSeq)
dittoPlot(obj, # seurat对象
          "OAT", # 你需要对比的gene
          group.by = "diseaseStatus", # 分组依据
          plots=c("vlnplot","boxplot"), # 绘制的图
          boxplot.fill=F,
          boxplot.color='white',
          color.panel = dittoColors(),
          colors = c(1,3),
          theme = theme(axis.text = element_text(size = 12, color = 'black'),
                        axis.line = element_line(size = 1),
                        axis.title.y = element_text(size = 15, color = 'black'),
                        plot.title = element_text(size=15,hjust=0.5, color = 'black')),
          ylab = 'Expression',
          y.breaks = seq(0,4,1),
          xlab = '',
          x.labels = c("HC","PD"),
          x.labels.rotate =F,
          max=4,
          min=0,
          main = "OAT",
          legend.show = F)

DimPlot(obj,group.by = "RNA_snn_res.0.7",pt.size = 0.6,label = T,split.by = "diseaseStatus")
DimPlot(obj,group.by = "orig.ident",pt.size = 0.6,label = T)
FeaturePlot(obj,features = "OAT",pt.size = 0.6,label = T)
FeaturePlot(obj,features = "LUM",pt.size = 0.6,label = T)
FeaturePlot(obj,features = "VWF",pt.size = 0.6,label = T)

library(scCustomize)
# Set color palette
pal <- viridis::viridis(n = 10, option = "D")
FeaturePlot_scCustom(seurat_object = obj, 
                     features = "OAT", order = F,colors_use = pal,
                     pt.size = 0.6,split.by = "diseaseStatus")

Idents(obj) <- obj$RNA_snn_res.0.7
Idents(obj) <- obj$RNA_snn_res.1.2
## 重新注释
library(COSG)
table(Idents(obj)) #聚类分群结果
marker_cosg <- cosg(
  obj,
  groups='all', #考虑全部分组
  assay='RNA',
  slot='data',
  mu=1,         #惩罚项参数，值越大
  remove_lowly_expressed=TRUE,   #是否过滤低表达基因
  expressed_pct=0.1,             #设置低表达的阈值
  n_genes_user=100      #每个cluster定义Top-N个marker gene
)

FeaturePlot(obj,features = "DPP8",pt.size = 0.6)
table(obj$RNA_snn_res.0.7)

obj$diseaseStatue <- ifelse(grepl(""))


obj2 = readRDS("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/sce_harmony_unannotated.rds")
library(Seurat)
library(Nebulosa)
library(ggnetwork)
library(dplyr)
plot_density(obj2, 
             features = c("OAT"),
             pal = 'magma', 
             raster = T,  # 多余参数
             size = 0.8,
             reduction = "umap") &
  theme_blank()&
  theme(legend.frame = element_rect(colour = "black"),
        legend.ticks = element_line(colour = "black", linewidth  = 0),
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.8, "cm"),
        legend.title = element_text(color = 'black', face = "bold", size=8))

plot_density(obj2, 
             features = c("DPP8"),
             pal = 'magma', 
             raster = T,  
             size = 0.8,
             reduction = "umap") &
  theme_blank()&
  theme(legend.frame = element_rect(colour = "black"),
        legend.ticks = element_line(colour = "black", linewidth  = 0),
        legend.key.width = unit(0.3, "cm"),
        legend.key.height = unit(0.8, "cm"),
        legend.title = element_text(color = 'black', face = "bold", size=8))

##------------------------------------------------------------------------------------
## 开始进行细胞注释
## 采用CSGO 注释 呈现采用气泡图 第一次全部大类注释 结合原文
## T B Plasama Endothelial  Neutropil Monocytic Fibroblasts Mast cell Epithelial MDSCs
## 细胞比例图也要采用
##------------------------------------------------------------------------------------
markers <- c(
  "CD3D","TRAC", # T 0 1 2 
  "MS4A1","BANK1","CD19", # B 10群
  "PECAM1","VWF","CDH5", # 4群 endo
  "FCGR3B","CSF3R", # Neu 3群+ 17群
  "FCGR3A","GNLY", # NK 6群
  "LYZ","MS4A6A", # 9+12 monocytes 22+11 mono/mac
  "LUM","DCN","COL1A1","THBS2", # 7群Fibroblasts
  "TPSAB1","CPA3","TPSB2", # 11群mastcell
  "XCR1","CLEC9A","IDO1", # cDC1s # 21
  "CD1C","CD1E","FCER1A", # cDC2s # 19
  "IGHG1","XBP1","JCHAIN", # 5 + 8 群plasma cells
  "LILRA4","CLEC4C", # 15群 pDC
  "MKI67","TOP2A", # 13群增殖 proliferative cells
  "KRT6A","KRT19", # 14 Epithelial
  "RGS5","ACTA2" # 16 SMC
)

DotPlot(obj,
        features = markers,
        group.by = "anntation")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
obj$anntation = factor(obj$anntation,levels = c("T cells",
                                                "B cells",
                                                "Endothelials",
                                                "Neutrophils",
                                                "NKs",
                                                "mono/mac",
                                                "Fibroblasts",
                                                "Mast cells",
                                                "cDC1s",
                                                "cDC2s",
                                                "Plasma cells",
                                                "pDCs",
                                                "Proliferative cells",
                                                "Epithelials",
                                                "SMCs",
                                                "unsure"))

obj$anntation = "mono/mac"
obj$anntation <- ifelse(obj$RNA_snn_res.0.7 %in% c("0","1","2"),"T cells",obj$anntation)
obj$anntation <- ifelse(obj$RNA_snn_res.0.7 %in% c("10"),"B cells",obj$anntation)
obj$anntation <- ifelse(obj$RNA_snn_res.0.7 %in% c("8","5"),"Plasma cells",obj$anntation)
obj$anntation <- ifelse(obj$RNA_snn_res.0.7 %in% c("11"),"Mast cells",obj$anntation)
obj$anntation <- ifelse(obj$RNA_snn_res.0.7 %in% c("7"),"Fibroblasts",obj$anntation)
obj$anntation <- ifelse(obj$RNA_snn_res.0.7 %in% c("3","17"),"Neutrophils",obj$anntation)
obj$anntation <- ifelse(obj$RNA_snn_res.0.7 %in% c("15"),"pDCs",obj$anntation)
obj$anntation <- ifelse(obj$RNA_snn_res.0.7 %in% c("13"),"Proliferative cells",obj$anntation)
obj$anntation <- ifelse(obj$RNA_snn_res.0.7 %in% c("16"),"SMCs",obj$anntation)
obj$anntation <- ifelse(obj$RNA_snn_res.0.7 %in% c("6"),"NKs",obj$anntation)
obj$anntation <- ifelse(obj$RNA_snn_res.0.7 %in% c("4"),"Endothelials",obj$anntation)
obj$anntation <- ifelse(obj$RNA_snn_res.0.7 %in% c("14"),"Epithelials",obj$anntation)
obj$anntation <- ifelse(obj$RNA_snn_res.0.7 %in% c("18","19"),"unsure",obj$anntation)
obj$anntation <- ifelse(obj$RNA_snn_res.1.2 %in% c("11","22"),"mono/mac",obj$anntation)
obj$anntation <- ifelse(obj$RNA_snn_res.1.2 %in% c("21"),"cDC1s",obj$anntation)
obj$anntation <- ifelse(obj$RNA_snn_res.1.2 %in% c("19"),"cDC2s",obj$anntation)
DimPlot(obj,group.by = "annotation",label = T,pt.size = 0.5)
Idents(obj)=obj$anntation
obj2 = subset(obj,idents = c("T cells",
                             "B cells",
                             "Endothelials",
                             "Neutrophils",
                             "NKs",
                             "mono/mac",
                             "Fibroblasts",
                             "Mast cells",
                             "cDC1s",
                             "cDC2s",
                             "Plasma cells",
                             "pDCs",
                             "Proliferative cells",
                             "Epithelials",
                             "SMCs"))
DimPlot(obj2,group.by = "anntation",label = T,split.by = "diseaseStatus")
saveRDS(obj2,"D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H1. harmony/GSE171213/annotation_finished.rds")
obj2 <- readRDS("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H1. harmony/GSE171213/annotation_finished.rds")
scplotter::CellStatPlot(obj2, ident = "anntation",group_by = "diseaseStatus", frac = "group",swap = TRUE, position = "stack")
scplotter::CellStatPlot(obj2,ident = "anntation",plot_type = "ring", group_by = "diseaseStatus",palette = "Spectral")
##---------------------------------------------------
## 细胞比例图(甜甜圈图) + UMAP美化图 + Featureplot美化图
##---------------------------------------------------
library(tidyverse)
library(data.table)
ratio <- table(obj2@meta.data$anntation,obj2@meta.data$diseaseStatus) %>% melt()
colnames(ratio) <- c("Cluster","Sample","Number")
for (i in 1:15){
  ratio$percentage[i] = ratio$Number[i]/7701
}
for (i in 16:31){
  ratio$percentage[i] = ratio$Number[i]/11265
}
ratio$Percentage <- round(as.numeric(ratio$percentage) * 100, 2)
ratio$Percentage <- sprintf("%.2f%%", ratio$percentage * 100)
ratio = ratio[-c(16,32),]

library(ggplot2)
library(cols4all)

P1 = ggplot(ratio, aes(x = 3, 
                       y = Percentage, 
                       fill = Cluster)) + 
  geom_col(width = 1.5, 
           color = 'white') + 
  facet_grid(.~Sample) 

P2 = P1 + coord_polar(theta = "y") + xlim(c(0.2, 3.8))

#挑选自定义配色：
c4a_gui()
mycol <- c4a('pastel',15)
mycol2 <- c4a('20',15)
mycol3 <- c4a('miller_stone',15)
#颜色、主题更改：
P3 <- P2 +
  scale_fill_manual(values = mycol2) +
  theme_void()+ #空白主题
  theme(
    strip.text.x = element_text(size = 14), 
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 14) 
  )
P3

P4 <- P3 +
  geom_text(aes(label = paste0(Percentage,'%')),
            position = position_stack(vjust = 0.5),
            size = 4)
P4 

##---------------------------------
## 主要在成纤维细胞和SMC细胞表达。
## 那么我可以研究，PD 和 HC 对表达量作比较
##----------------------------------
Idents(obj2) <- obj2$anntation
obj3 = subset(obj2,idents = c("Endothelials"))
library(dittoSeq)
dittoPlot(obj2, 
          "OAT", 
          group.by = "diseaseStatus", 
          plots=c("vlnplot"),
          boxplot.fill=F,
          boxplot.color='white',
          color.panel = dittoColors(),
          colors = c(1,3),
          theme = theme(axis.text = element_text(size = 12, color = 'black'),
                        axis.line = element_line(size = 1),
                        axis.title.y = element_text(size = 15, color = 'black'),
                        plot.title = element_text(size=15,hjust=0.5, color = 'black')),
          ylab = 'Expression',
          y.breaks = seq(0,4,1),
          xlab = '',
          x.labels = c("HC","PD"),
          x.labels.rotate =F,
          max=4,
          min=0,
          main = "DPP8",
          legend.show = F) +stat_compare_means(method = "wilcox.test", label = "p.signif")

FeaturePlot(obj2, 
            features = "OAT", order = T,
            pt.size = 0.6,split.by = "diseaseStatus")


library(scCustomize)
sample_colors <- c("dodgerblue", "firebrick1")

# Create Plots
Stacked_VlnPlot(seurat_object = obj2, 
                features = c("DPP8","OAT"), 
                x_lab_rotate = TRUE,
                colors_use = sample_colors, 
                split.by = "diseaseStatus")

Stacked_VlnPlot(seurat_object = obj2, 
                features = c("DPP8","OAT"), 
                x_lab_rotate = TRUE,
                colors_use = sample_colors, 
                split.by = "diseaseStatus")
Stacked_VlnPlot(seurat_object = obj2, 
                features = c("DPP8","OAT"), 
                x_lab_rotate = TRUE,
                plot_spacing = 0.3)

VlnPlot_scCustom(seurat_object = obj2, features = "DPP8", plot_median = TRUE) & NoLegend()
VlnPlot_scCustom(seurat_object = obj2, features = "OAT", plot_boxplot = TRUE) & NoLegend()

##函数区域
scRNAdataPreProcessing <- function(
    obj, objPrefix, plotDir, 
    minFeatures=200, maxFeatures=Inf, minCounts=1000, maxCounts=Inf, maxPctMito=10, 
    nfeatures=2500, dims=1:15, res=0.5, 
    runDoubletFinder=TRUE, 
    estDubRate=0.075, 
    runDecontX=TRUE, 
    assays="RNA",
    ncores=1, use_logfile=TRUE
){
  
  # Perform some basic filtering on a seurat object.
  # Seurat object should be created from a single sample (e.g. one 10x run)
  # Optionally run other pre-processing tools:
  #
  # - DoubletFinder to estimate likely doublets 
  #   https://github.com/chris-mcginnis-ucsf/DoubletFinder
  #   https://www-cell-com.stanford.idm.oclc.org/cell-systems/fulltext/S2405-4712(19)30073-0
  #
  # - DecontX to reduce ambient RNA contamination
  #   https://github.com/campbio/celda
  #   https://genomebiology.biomedcentral.com/articles/10.1186/s13059-020-1950-6
  # 
  ################################
  # obj = seurat object
  # objPrefix = prefix for raq_qc plots
  # plotDir = directory for plotting
  # minFeatures = minimum number of features (genes) per cell
  # maxFeatures = maximum number of features (genes) per cell
  # minCounts = minimum number of UMIs per cell
  # maxCounts = maximum number of UMIs per cell
  # maxPctMito = maximum mitochondrial read percentage
  # nfeatures = number of variable features to be used in DoubletFinder
  # dims = which PCA dimensions will be used in DoubletFinder
  # estDubRate = estimated percent doublets (DF.classify will identify exactly this percent as doublets!)
  
  if(use_logfile){
    logfile <- paste0(plotDir, sprintf("/%s_preprocess_log_%s.txt", objPrefix, format(Sys.time(), "%Y%m%d-%H%M%S")))
    con <- file(logfile, open = "wt")
    sink(con, type="output")
    sink(con, type="message")
  }
  
  # Add percent.mt percent.rb percent.hb cellcycle
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj[["percent.rb"]] <- PercentageFeatureSet(obj, pattern = "^RP[SL]")
  obj[["percent.hb"]] <- PercentageFeatureSet(obj, pattern = "^HB[^(P)]") 
  s.genes = Seurat::cc.genes.updated.2019$s.genes
  g2m.genes = Seurat::cc.genes.updated.2019$g2m.genes
  obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features =  g2m.genes, set.ident = TRUE)
  cellsBeforeFiltering <- dim(obj)[2]
  
  # Save some quick plots of raw qc
  histBreaks <- 100
  
  pdf(paste0(plotDir, sprintf("/%s_nCountHist.pdf", objPrefix)))
  df <- data.frame(cells=Cells(obj), log10nCount=log10(obj$nCount_RNA))
  p <- qcHistFilter(df, cmap = "blue", bins=histBreaks, border_color="black", lower_lim=log10(minCounts), upper_lim=log10(maxCounts))
  print(p)
  dev.off()
  
  pdf(paste0(plotDir, sprintf("/%s_nFeatureHist.pdf", objPrefix)))
  df <- data.frame(cells=Cells(obj), log10nFeatures=log10(obj$nFeature_RNA))
  p <- qcHistFilter(df, cmap = "blue", bins=histBreaks, border_color="black", lower_lim=log10(minFeatures), upper_lim=log10(maxFeatures))
  print(p)
  dev.off()
  
  pdf(paste0(plotDir, sprintf("/%s_pctMitoHist.pdf", objPrefix)))
  df <- data.frame(cells=Cells(obj), PctMito=obj$percent.mt)
  p <- qcHistFilter(df, cmap = "blue", bins=histBreaks, border_color="black", upper_lim=maxPctMito)
  print(p)
  dev.off()
  
  # Perform basic hard threshold filters
  message("Will filter based on:")
  message(sprintf("%s < unique genes < %s", minFeatures, maxFeatures))
  message(sprintf("%s < UMIs (counts) < %s", minCounts, maxCounts))
  message(sprintf("Percent mitochondrial < %s", maxPctMito))
  
  obj <- subset(obj, 
                subset = (
                  nFeature_RNA > minFeatures & 
                    nFeature_RNA < maxFeatures & 
                    nCount_RNA > minCounts & 
                    nCount_RNA < maxCounts & 
                    percent.mt < maxPctMito
                )
  )
  cellsAfterFiltering <- dim(obj)[2]
  message(sprintf("%s filtered down to %s (%s%% remaining)", 
                  cellsBeforeFiltering, cellsAfterFiltering, 
                  round(100*(cellsAfterFiltering/cellsBeforeFiltering), 2)))
  
  # Perform standard Seurat pre-processing:
  obj <- seuratPreProcess(obj, selectionMethod="vst", nFeatures = nfeatures, dims = dims)
  
  # Run DoubletFinder, if indicated
  if(runDoubletFinder){
    message("Running DoubletFinder...")
    obj <- runDoubletFinder(obj, dims, estDubRate=estDubRate, ncores=ncores)
  }
  
  # Run DecontX, if indicated
  if(runDecontX){
    message("Running DecontX...")
    obj <- runDecontX(obj)
    assays <- c(assays, "origCounts")
  }
  
  # Return filtered Seurat object
  obj <- DietSeurat(obj, counts=TRUE, data=TRUE, scale.data=FALSE, assays=assays)
  
  # Close connections
  message("Finished preprocessing...")
  if(use_logfile){
    on.exit({ sink(type = "message"); sink(type = "output"); close(con) })
  }
  
  # Return processed Seurat object
  return(obj)
}
runDoubletFinder <- function(obj, dims, estDubRate=0.075, ncores=1){
  # Run DoubletFinder on a provided (preprocessed) Seurat object
  # Return the seurat object with the selected pANN parameter and the 
  # DoubletFinder doublet classifications
  
  ### pK Identification (parameter-sweep) ###
  # "pK ~ This defines the PC neighborhood size used to compute pANN (proportion of artificial nearest neighbors), 
  # expressed as a proportion of the merged real-artificial data. 
  # No default is set, as pK should be adjusted for each scRNA-seq dataset"
  
  sweep.res.list <- paramSweep(obj, PCs=dims, sct=FALSE, num.cores=ncores)
  sweep.stats <- summarizeSweep(sweep.res.list, GT=FALSE)
  bcmvn <- find.pK(sweep.stats)
  pK <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
  message(sprintf("Using pK = %s...", pK))
  
  # Get expected doublets (DF.classify will identify exactly this percent as doublets!)
  nExp_poi <- round(estDubRate * length(Cells(obj)))
  
  # DoubletFinder:
  obj <- doubletFinder(obj, PCs = dims, pN = 0.25, pK = pK, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
  
  # Rename results into more useful annotations
  pann <- grep(pattern="^pANN", x=names(obj@meta.data), value=TRUE)
  message(sprintf("Using pANN = %s...", pann))
  classify <- grep(pattern="^DF.classifications", x=names(obj@meta.data), value=TRUE)
  obj$pANN <- obj[[pann]]
  obj$DF.classify <- obj[[classify]]
  obj[[pann]] <- NULL
  obj[[classify]] <- NULL
  
  return(obj)
}
runDecontX <- function(obj, seed=1){
  # Run DecontX on a provided Seurat object
  # From the DecontX vignette: 
  # "**Only the expression profile of *"real"* cells after cell calling are required to run DecontX. 
  # Empty cell droplet information (low expression cell barcodes before cell calling) are not needed.**"
  
  # DecontX can take either `SingleCellExperiment` object... or a single counts matrix as input. 
  # `decontX` will attempt to convert any input matrix to class `dgCMatrix` before beginning any analyses.
  counts <- GetAssayData(object = obj, slot = "counts")
  clusters <- Idents(obj) %>% as.numeric()
  
  # Run on only expressed genes
  x <- counts[rowSums(counts)>0,]
  message(sprintf("Running decontX on %s cells with %s non-zero genes...", dim(x)[2], dim(x)[1]))
  decon <- decontX(x, z=clusters, verbose=TRUE, seed=seed)
  
  # Save desired information back to Seurat Object
  # We will place the estimated 'decontaminated counts' in place of the original counts ('RNA')
  # and keep the original counts as a separate assay called 'origCounts'
  obj[["origCounts"]] <- CreateAssayObject(counts = counts)
  newCounts <- decon$decontXcounts
  # Add back unexpressed genes and sort according to original counts
  counts_rows <- rownames(counts) ## rownames1
  newCounts_rows <- rownames(newCounts) ## rownames2
  left <- setdiff(counts_rows,newCounts_rows) ## DEG
  left_matrix <- counts[left,, drop = FALSE]
  row.names(left_matrix) <- left
  newCounts <- rbind2(newCounts, left_matrix)[counts_rows,]
  identical(rownames(newCounts), counts_rows)
  
  obj[["RNA"]]@counts <- as(round(newCounts), "sparseMatrix")
  obj$estConp <- decon$contamination # Estimated 'contamination proportion, 0 to 1'
  
  return(obj)
}
