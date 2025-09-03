library(clusterProfiler) # GO enrichment analysis, KEGG pathway enrichment analysis
library(org.Hs.eg.db)    # Gene annotation database
library(enrichplot)
library(ggplot2)
library(GOplot)
library(dplyr)

# Set working directory
setwd("working directory")
# Read input data
genes_df <- read.csv("core_gene_set.csv")  

# Convert gene symbols to ENTREZ IDs
entrezIDs <- bitr(genes_df$gene, fromType = "SYMBOL", 
                  toType = c("ENTREZID", "SYMBOL"),
                  OrgDb = org.Hs.eg.db) # `OrgDb` specifies the human genome annotation database


# Use entrezIDs 
gene<- entrezIDs$ENTREZID
## GO enrichment analysis
go<- enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, qvalueCutoff = 0.05, ont="all", readable =T)
write.table(go,file="GO.txt",sep="\t",quote=F,row.names = F)

## Visualization
## Barplot
pdf(file="GO-barplot.pdf",width = 8,height = 10)
## `showCategory` controls the number of terms displayed
barplot(go, drop = TRUE, showCategory =6, label_format=100, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

## Bubble plot
pdf(file="GO-bubbleplot.pdf",width = 6,height = 6)
dotplot(go, showCategory = 3, label_format=100, split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(dplyr)

selected_paths <- c(
  "arginine catabolic process",
  "proline metabolic process",
  "glutamine family amino acid biosynthetic process",
  "dipeptidyl-peptidase activity",
  "transaminase activity",
  "transferase activity, transferring nitrogenous groups",
  "mitochondrial matrix"
)


# Filter enrichResult object, retain enrichResult structure
go_selected <- go

# Highlight selected pathways (**only works for BP**)

# **Create y-axis label wrapping**
library(stringr)
go_selected@result <- go@result %>% filter(Description %in% selected_paths)
# Set line wrap width (you can adjust according to plot width, e.g., width = 40)
go_selected@result$Description <- str_wrap(go_selected@result$Description, width = 40)
go_selected@result$Description <- str_to_sentence(go_selected@result$Description)

# Generate GO bubble plot
pdf(file = "GO_selected_pathways.pdf", width = 6, height = 3)

p <- dotplot(go_selected,
             showCategory = length(selected_paths),
             label_format = 100,
             split = "ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scale = 'free') +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 10),
    strip.text = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.spacing = unit(0.3, "lines"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10)
  )

print(p)
dev.off()


# KEGG analysis
kk <- enrichKEGG(gene = gene, keyType = "kegg", organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05, pAdjustMethod = "fdr")   
write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)                         

# Load ReactomePA package
library(ReactomePA)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ReactomePA")

# Reactome enrichment analysis
reactome_result <- enrichPathway(
  gene         = gene,         # Note: input must be Entrez IDs
  organism     = "human",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable     = TRUE
)

## Visualization
## Barplot
pdf(file="KEGG-barplot.pdf",width = 8,height = 6)
barplot(kk, drop = TRUE, showCategory = 8, label_format=100)
dev.off()

## Bubble plot
pdf(file="KEGG-bubbleplot.pdf",width = 6,height = 4)
dotplot(kk, showCategory = 8, label_format=100)
dev.off()


# Load required packages
library(clusterProfiler)
library(ggplot2)
library(dplyr)

# Set KEGG pathways to display
selected_kegg <- c(
  "Arginine and proline metabolism"
)

# Filter enrichResult object, retain enrichResult structure
kk_selected <- kk
kk_selected@result <- kk@result %>% filter(Description %in% selected_kegg)
kk_selected@result$Description <- str_wrap(kk_selected@result$Description, width = 40)
kk_selected@result$Description <- str_to_sentence(kk_selected@result$Description)

# Generate KEGG bubble plot
pdf(file = "KEGG-bubble_selected.pdf", width = 6, height = 3)

p <- dotplot(kk_selected, showCategory = length(selected_kegg), label_format = 100)+ 
  theme_minimal()+  # Make the plot cleaner
  theme( 
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 10),
    strip.text = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    panel.spacing = unit(0.3, "lines"),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10)
  )

print(p)
dev.off()
