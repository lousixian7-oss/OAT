##----Step 1.Data Preprocess----
library(data.table)
library(Seurat)

# Define input directory and sample IDs
base_dir <- "D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H1. harmony/GSE171213/"
sample_ids <- c("HC1","HC2","HC3","HC4","PD1","PD2","PD3","PD4","PD5")

# Initialize an empty list to store Seurat objects
scRNAlist <- list()

# Loop through each sample and construct Seurat objects
for (sample in sample_ids) {
  message("Processing ", sample, " ...")
  
  # Define file paths for count matrix and cell name list
  count_file <- file.path(base_dir, sample, list.files(file.path(base_dir, sample), pattern="counts.tsv.gz"))
  name_file  <- file.path(base_dir, sample, list.files(file.path(base_dir, sample), pattern="cellname.list.txt.gz"))
  
  # Load raw count data (gene × cell index matrix) and cell name mapping
  df <- fread(count_file)   # raw counts, columns = gene + cell indices (C1, C2, …)
  name <- fread(name_file)  # mapping table: CellName ↔ CellIndex
  
  # Replace cell indices (C1, C2, …) with actual cell barcodes
  cell_index <- colnames(df)[-1]                 # exclude "gene" column
  map <- setNames(name$CellName, name$CellIndex) # create mapping
  new_colnames <- c("gene", map[cell_index])     # assign new colnames
  colnames(df) <- new_colnames
  
  # Save processed count matrix (optional)
  write.csv(df, file = file.path(base_dir, sample, paste0(sample, "_count.csv")), row.names = FALSE)
  
  # Convert to numeric matrix (genes as rownames, cells as colnames)
  mat <- as.matrix(df[,-1])    # exclude "gene" column
  rownames(mat) <- df$gene
  
  # Create Seurat object with basic QC filters
  obj <- CreateSeuratObject(
    counts = mat, 
    project = sample,
    min.cells = 3, #  retain genes expressed in ≥3 cells
    min.features = 200 # retain cells with ≥200 detected genes
  )
  
  # Store object into list
  scRNAlist[[sample]] <- obj
}

# Save final list of Seurat objects for downstream analysis
saveRDS(scRNAlist, file.path(base_dir, "scRNAlist.rds"))



##----Step 2. Load Seurat object list and preprocessing----
scRNAlist = readRDS("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H1. harmony/GSE171213/scRNAlist.rds")
wd = getwd()  # working directory
objNames <- sapply(scRNAlist, function(x) x@project.name)   # Extract project names
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
    runDoubletFinder = T, # Optional processing steps
    runDecontX = T, # Optional processing steps
    estDubRate = 0.02, 
    ncores = 1,
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


##----Step 3. Harmony integration and clustering----
library(harmony)
library(clustree)
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

DimPlot(obj,group.by = "RNA_snn_res.0.7",pt.size = 0.6,label = T,split.by = "diseaseStatus")
DimPlot(obj,group.by = "orig.ident",pt.size = 0.6,label = T)
FeaturePlot(obj,features = "OAT",pt.size = 0.6,label = T)

##----Step 4. Cell Type Annotation-------
Idents(obj) <- obj$RNA_snn_res.0.7
table(obj$RNA_snn_res.0.7)
## 重新注释
library(COSG)
table(Idents(obj)) 
marker_cosg <- cosg(
  obj,
  groups='all', 
  assay='RNA',
  slot='data',
  mu=1,         
  remove_lowly_expressed=TRUE,   
  expressed_pct=0.1,             
  n_genes_user=100      
)

FeaturePlot(obj,features = "DPP8",pt.size = 0.6)
FeaturePlot(obj,features = "OAT",pt.size = 0.6)


markers <- c(
  "CD3D","TRAC", # T
  "MS4A1","BANK1","CD19", # B
  "PECAM1","VWF","CDH5", # endo
  "FCGR3B","CSF3R", # Neu 
  "FCGR3A","GNLY", # NK 
  "LYZ","MS4A6A", # mono/mac
  "LUM","DCN","COL1A1","THBS2", # Fibroblasts
  "TPSAB1","CPA3","TPSB2", # mastcell
  "XCR1","CLEC9A","IDO1", # cDC1s 
  "CD1C","CD1E","FCER1A", # cDC2s 
  "IGHG1","XBP1","JCHAIN", # plasma cells
  "LILRA4","CLEC4C", #  pDC
  "MKI67","TOP2A", #  proliferative cells
  "KRT6A","KRT19", #  Epithelial
  "RGS5","ACTA2" #  SMC
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
##----Step 5. Visualization with SCP & scplotter----
## 1.1 UMAP 美化
obj <- readRDS("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/GSE171213/annotation_finished.rds")
Idents(obj) <- obj$anntation

# 1. UMAP plots
P1 <- scplotter::CellDimPlot(obj,
                             group_by = "anntation", 
                             #highlight = 'anntation == "Fibroblasts"',
                             #highlight = TRUE, 
                             theme = "theme_blank",
                             #label = TRUE,
                             #label_size = 3,
                             #label_fg = "orange",
                             #label_bg = "red",
                             theme_args = list(base_size = 16))
P1 

?CellDimPlot
P2 <- scplotter::CellDimPlot(obj, group_by = "anntation", reduction = "umap")
P2

pal7 <- SCP::palette_scp(1:15, palette = "Paired", type = "discrete")
SCP::show_palettes(pal7)
cluster <- pal7
cluster_col <- c("#A6CEE3","#3B8ABE","#72B29C","#84C868","#4FA435", "#EEBC6A","#FE911F","#FD8C4C",
                 "#F47575","#E12429","#CD9CBB","#8C66AF","#A99099","#EEDB80","#B15928")



# 2. Proportion plots
?CellStatPlot
Idents(obj) <- obj$anntation
P1.2.1 <- scplotter::CellStatPlot(obj, plot_type = "ring",ident = "anntation",group_by = "diseaseStatus",palette = "Spectral")
P1.2.2 <- scplotter::CellStatPlot(obj, plot_type = "pie",ident = "anntation", split_by = "diseaseStatus")


# 3. Feature expression plots
P1.3.1 <- SCP::FeatureDimPlot(obj, 
                              features = "OAT", 
                              reduction = "umap", 
                              pt.size = 0.5,
                              theme_use = ggplot2::theme_classic, 
                              theme_args = list(base_size = 16),
                              #lower_cutoff =  1,
                              #upper_cutoff = 2
)
P1.3.2 <- SCP::FeatureDimPlot(obj, 
                              features = "DPP8", 
                              reduction = "umap", 
                              pt.size = 0.5,
                              theme_use = ggplot2::theme_classic, 
                              theme_args = list(base_size = 16),
                              #lower_cutoff =  1,
                              #upper_cutoff = 2
)

P1.3.3 <- scplotter::FeatureStatPlot(obj, 
                                     ident = "anntation",
                                     features = "OAT", 
                                     reduction = "umap", 
                                     pt.size = 0.5,
                                     theme_use = ggplot2::theme_classic, 
                                     theme_args = list(base_size = 16),
                                     highlight = 'anntation == "Fibroblasts"')

P1.3.4 <- scplotter::FeatureStatPlot(obj, 
                                     ident = "anntation",
                                     features = "DPP8", 
                                     reduction = "umap", 
                                     pt.size = 0.5,
                                     theme_use = ggplot2::theme_classic, 
                                     theme_args = list(base_size = 16),
                                     highlight = 'anntation == "Endothelials"'
)

# 4. DotPlot of cell type markers
markers <- c(
  "CD3D","TRAC", 
  "MS4A1","BANK1","CD19",
  "PECAM1","VWF","CDH5", 
  "FCGR3B","CSF3R", 
  "FCGR3A","GNLY", 
  "LYZ","MS4A6A", 
  "LUM","DCN","COL1A1","THBS2", 
  "TPSAB1","CPA3","TPSB2", 
  "XCR1","CLEC9A","IDO1", 
  "CD1C","CD1E","FCER1A", 
  "IGHG1","XBP1","JCHAIN",
  "LILRA4","CLEC4C", 
  "MKI67","TOP2A",
  "KRT6A","KRT19", 
  "RGS5","ACTA2" 
)

library(ggh4x)
library(ggplot2)
library(dplyr)
df <- obj@reductions$umap@cell.embeddings%>% 
  as.data.frame() %>%
  cbind(clusters = obj@meta.data$RNA_snn_res.0.7)%>%
  cbind(celltype = obj@meta.data$anntation)
Idents(obj)=obj$anntation
p = DotPlot(obj, 
            features = markers,
            col.min = 0,
            group.by = "anntation")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dat <- p$data

anno <- distinct(df, obj$RNA_snn_res.0.7,obj$anntation)
colnames(anno) <- c("id","celltype")

df <- left_join(dat, anno, by = "id")

P1.4.2 <- ggplot(dat, aes(features.plot, id,size=pct.exp, fill=avg.exp.scaled)) + 
  geom_point(shape = 21, colour="black", stroke=0.5) +#气泡
  guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill=NA))) + #legend
  theme(
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA),
    panel.grid.major.x = element_line(color = "grey80"),
    panel.grid.major.y = element_line(color = "grey80"),
    axis.title = element_blank(),
    axis.text.y = element_text(color='black',size=12),
    axis.text.x = element_text(color='black',size=12, angle = 90, hjust = 1, vjust = 0.5))+
  scale_fill_gradientn(colours = c('#5749a0', '#0f7ab0', '#00bbb1',
                                   '#bef0b0', '#fdf4af', '#f9b64b',
                                   '#ec840e', '#ca443d', '#a51a49'))
P1.4.2


















































































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
