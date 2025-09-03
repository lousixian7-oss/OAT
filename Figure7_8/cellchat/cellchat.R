##----Defining OAT+/OAT− Fibroblast Subgroups----
##---------------------------------------
## Step 1. Check fibroblast subgroup annotations
##---------------------------------------
table(obj.fib@meta.data$annotation2)
# OAT- Fib   OAT+ Fib 
#   294        712

##---------------------------------------
## Step 2. Extract cell IDs for each subgroup
##---------------------------------------
sample <- row.names(obj.fib@meta.data)    # Extract all cell IDs
OAT1 <- row.names(obj.fib@meta.data)[obj.fib$annotation2 == "OAT- Fib"]  # IDs of OAT- fibroblasts
OAT2 <- row.names(obj.fib@meta.data)[obj.fib$annotation2 == "OAT+ Fib"]  # IDs of OAT+ fibroblasts

##---------------------------------------
## Step 3. Create a new annotation column in the Seurat object
##---------------------------------------
obj$annotation_new <- obj$anntation         # Copy the existing cluster annotation
obj$annotation_new <- as.character(obj$annotation_new)  # Convert factor → character
obj$names <- row.names(obj@meta.data)       # Store cell IDs in a new column "names"

##---------------------------------------
## Step 4. Relabel OAT+ and OAT− fibroblast cells
##---------------------------------------
obj$annotation_new[obj$names %in% OAT1] <- "OAT- Fib"   # Mark OAT- fibroblasts
obj$annotation_new[obj$names %in% OAT2] <- "OAT+ Fib"   # Mark OAT+ fibroblasts
obj$annotation_new <- as.factor(obj$annotation_new)     # Convert back to factor







############################################################
##----Creating and Saving CellChat Objects for HC and PD----
##---------------------------------------
## Step 1. Load required packages
##---------------------------------------
library(CellChat)
library(Seurat)
##---------------------------------------
## Step 2. Split dataset into HC and PD groups
##---------------------------------------
Idents(obj) <- obj$diseaseStatus              # Set cell identities to disease status (HC / PD)
obj.HC = subset(obj, idents = "HC")           # Subset fibroblast cells from healthy controls
obj.PD = subset(obj, idents = "PD")           # Subset fibroblast cells from patients

Idents(obj.HC) <- obj.HC$anntation            # Assign cluster identities for HC group
Idents(obj.PD) <- obj.PD$anntation            # Assign cluster identities for PD group
##---------------------------------------
## Step 3. Initialize CellChat objects
##---------------------------------------
cellchat.HC <- createCellChat(object = obj.HC@assays[["RNA"]]@data, 
                              meta = obj.HC@meta.data, 
                              group.by = "annotation_new")

cellchat.PD <- createCellChat(object = obj.PD@assays[["RNA"]]@data, 
                              meta = obj.PD@meta.data, 
                              group.by = "annotation_new")
##---------------------------------------
## Step 4. Run CellChat workflow for HC
##---------------------------------------
cellchat = cellchat.HC
CellChatDB <- CellChatDB.human                              # Load human ligand–receptor database
cellchat@DB  <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# Focus on secreted signaling pathways
cellchat <- subsetData(cellchat)                           # Extract relevant signaling genes
cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE) 
# Identify overexpressed genes
cellchat <- identifyOverExpressedInteractions(cellchat)    # Identify overexpressed ligand–receptor pairs
cellchat <- projectData(cellchat, PPI.human)               # Map signaling interactions onto PPI network
cellchat@idents <- droplevels(cellchat@idents)             # Remove unused cell identity levels
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE) 
# Compute communication probability
cellchat <- filterCommunication(cellchat, min.cells = 3)   # Remove low-confidence interactions
cellchat <- computeCommunProbPathway(cellchat)             # Infer pathway-level communications
cellchat <- aggregateNet(cellchat)                         # Aggregate communication networks
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
# Compute network centrality
cellchat.HC = cellchat                                     # Save HC result
##---------------------------------------
## Step 5. Run CellChat workflow for PD
##---------------------------------------
cellchat = cellchat.PD
CellChatDB <- CellChatDB.human
cellchat@DB  <- subsetDB(CellChatDB, search = "Secreted Signaling")
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat, do.fast = FALSE)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat@idents <- droplevels(cellchat@idents)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
cellchat.PD = cellchat                                     # Save PD result
##---------------------------------------
## Step 6. Save and reload results
##---------------------------------------
saveRDS(cellchat.HC, "cellchat.HC.rds")                    # Save HC CellChat object
saveRDS(cellchat.PD, "cellchat.PD.rds")                    # Save PD CellChat object

cellchat.HC <- readRDS("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H10. Final/cellchat.HC.rds")
cellchat.PD <- readRDS("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H10. Final/cellchat.PD.rds")
##---------------------------------------
## Step 7. Merge HC and PD CellChat objects
##---------------------------------------
cco.list <- list(HC = cellchat.HC, 
                 PD = cellchat.PD)

cellchat_merge <- mergeCellChat(cco.list, 
                                add.names = names(cco.list), 
                                cell.prefix = TRUE)





############################################################
##----Visualization and Correlation Analysis of HC vs PD CellChat Results----

##---------------------------------------
## Step 1. Interaction network visualization
##---------------------------------------

# Compare number (count) and strength (weight) of cell-cell interactions
a1 <- compareInteractions(cellchat_merge, show.legend = F, group = c(1,2), measure = "count")   
a2 <- compareInteractions(cellchat_merge, show.legend = F, group = c(1,2), measure = "weight")  
cellchat.p1 <- a1 + a2
cellchat.p1   # Barplots comparing HC vs PD

# Heatmaps of interaction probability (counts vs weights)
par(mfrow = c(1,1))
b1 <- netVisual_heatmap(cellchat_merge)                          
b2 <- netVisual_heatmap(cellchat_merge, measure = "weight")      
cellchat.p2 <- b1 + b2
cellchat.p2

# Differential interaction network between HC and PD
par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat_merge, weight.scale = T, comparison = c(1,2))
netVisual_diffInteraction(cellchat_merge, weight.scale = T, measure = "weight")
c1 <- netVisual_diffInteraction(cellchat_merge, weight.scale = T, comparison = c(1,2))
c2 <- netVisual_diffInteraction(cellchat_merge, weight.scale = T, measure = "weight")
cellchat.p3 <- c1 + c2
cellchat.p3

##---------------------------------------
## Step 2. Pathway activity comparison
##---------------------------------------

# Compare signaling pathway activities between HC and PD
d1 <- rankNet(cellchat_merge, mode = "comparison", stacked = T, do.stat = TRUE)  
d2 <- rankNet(cellchat_merge, mode = "comparison", stacked = F, do.stat = TRUE)  
cellchat.p4 <- d1 + d2
cellchat.p4

# Scatter plots showing signaling role analysis
cellchat.p5.1 <- netAnalysis_signalingRole_scatter(cellchat.HC)
cellchat.p5.2 <- netAnalysis_signalingRole_scatter(cellchat.PD)

##---------------------------------------
## Step 3. Bubble plots of ligand–receptor pairs
##---------------------------------------

levels(cellchat.PD@idents)   # 查看 PD 中的细胞亚群标签

# Bubble plot: selected source and target clusters
netVisual_bubble(cellchat.PD, sources.use = 11, targets.use = c(4,1,16,7,8), remove.isolate = FALSE)

netVisual_bubble(cellchat.PD, sources.use = c(3,7), targets.use = c(1:5), remove.isolate = FALSE)

# Bubble plot for specific signaling pathways (WNT, FASLG, NRG, CHEMERIN)
cellchat.PD@netP$pathways    
netVisual_bubble(cellchat.PD, sources.use = c(1:16), targets.use = c(11), 
                 signaling = c("WNT", "FASLG", "NRG", "CHEMERIN"), 
                 remove.isolate = FALSE)

##---------------------------------------
## Step 4. Chord plots and gene expression violin plots
##---------------------------------------

# Chord diagrams for specific signaling pathways in PD
cellchat.p6.1 <- netVisual_aggregate(cellchat.PD, signaling ="WNT", layout = "chord")
cellchat.p6.2 <- netVisual_aggregate(cellchat.PD, signaling ="FASLG", layout = "chord")
cellchat.p6.3 <- netVisual_aggregate(cellchat.PD, signaling ="NRG", layout = "chord")
cellchat.p6.4 <- netVisual_aggregate(cellchat.PD, signaling ="CHEMERIN", layout = "chord")

# Violin plots showing gene expression of pathway genes
cellchat.p7.1 <- plotGeneExpression(cellchat.PD, signaling = "WNT")
cellchat.p7.2 <- plotGeneExpression(cellchat.PD, signaling = "FASLG")
cellchat.p7.3 <- plotGeneExpression(cellchat.PD, signaling = "NRG")
cellchat.p7.4 <- plotGeneExpression(cellchat.PD, signaling = "CHEMERIN")

# 
cellchat.p8.1 <- netAnalysis_signalingRole_network(cellchat.PD, signaling = "WNT", width = 8, height = 2.5, font.size = 10)
cellchat.p8.2 <- netAnalysis_signalingRole_network(cellchat.PD, signaling = "FASLG", width = 8, height = 2.5, font.size = 10)
cellchat.p8.3 <- netAnalysis_signalingRole_network(cellchat.PD, signaling = "NRG", width = 8, height = 2.5, font.size = 10)
cellchat.p8.4 <- netAnalysis_signalingRole_network(cellchat.PD, signaling = "CHEMERIN", width = 8, height = 2.5, font.size = 10)



##---------------------------------------
## Step 5. Correlation analysis of OAT with other genes
##---------------------------------------
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(ggExtra)

# Load fibroblast object and subset PD group
Idents(obj.fib) <- obj.fib$diseaseStatus
obj.PD <- subset(obj.fib, idents = c("PD"))

# Select genes of interest
genes_of_interest <- c("OAT", "CXCL1", "CXCL13","CXCL6","CXCL2","CXCL5","PDCD1",
                       "FZD4","FZD6","LRP5","WNT2","FAS","DCN","LAMP5",
                       "ARG1","ARG2","SOD2","ITGA8","RARRES2","CMKLR1")

# Extract expression matrix for selected genes
expression_matrix <- GetAssayData(object = obj.fib, slot = "data")
extracted_matrix <- expression_matrix[genes_of_interest, ]
extracted_df <- as.data.frame(t(extracted_matrix))
extracted_df <- FetchData(obj.PD, vars = genes_of_interest, slot = "data")

# Example correlation plot: OAT vs RARRES2
p <- ggplot(extracted_df, aes_string(x = 'OAT', y = 'RARRES2')) +
  geom_point(size = 2, color = '#C1E0DB', alpha = 0.5) +
  theme_bw() +
  theme(axis.title = element_text(size = 16),
        axis.text = element_text(size = 14),
        axis.ticks.length = unit(0.25,'cm'),
        axis.ticks = element_line(size = 1),
        panel.border = element_rect(size = 1.5),
        panel.grid = element_blank()) +
  geom_smooth(method = 'lm', se = T, color = '#F9B208', size = 1.5, fill = '#FEA82F') +
  stat_cor(method = "pearson", digits = 3, size = 6)

# Add marginal density histograms
p2 <- ggMarginal(p, type = "densigram",
                 xparams = list(binwidth = 0.1, fill = "#B3E283", size = .7),
                 yparams = list(binwidth = 0.1, fill = "#8AB6D6", size = .7))
