##############################################
## pseudo-bulk analysis
##############################################

##----Step 1. Prepare pseudo-bulk counts----
# Create sample-level groups by combining sample ID and OAT status
obj.fib <- readRDS("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/obj.fib.modified.rds")
obj.fib$sample_group <- paste0(obj.fib$orig.ident, "_", obj.fib$annotation2)
table(obj.fib$sample_group)

# Aggregate raw counts at the sample_group level
pb_expr <- AggregateExpression(
  obj.fib,
  group.by = "sample_group",
  assays = "RNA",
  slot = "counts",    # Must use raw counts for DESeq2
  return.seurat = FALSE
)[[1]] # Extract matrix
dim(pb_expr)
pb_expr[1:3, 1:3]


##----Step 2. Build metadata----
library(tibble)
meta <- data.frame(sample_group = colnames(pb_expr)) %>%
  mutate(
    sample = gsub("_.*", "", sample_group),          # Extract sample ID
    OAT_status = gsub(".*_", "", sample_group),      # Extract OAT label
    OAT_status = ifelse(OAT_status == "OAT+ Fib", "OAT_pos", "OAT_neg")
  )
rownames(meta) <- meta$sample_group
head(meta)
table(meta$OAT_status)


##-----Step 3. Differential expression with DESeq2-----
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = pb_expr, colData = meta, design = ~ OAT_status)
# Filter out lowly expressed genes
dds <- dds[rowSums(counts(dds)) >= 10, ]
# Run DESeq2
dds <- DESeq(dds)
res <- results(dds, contrast = c("OAT_status", "OAT_pos", "OAT_neg"))
resOrdered <- res[order(res$padj), ]
DEG_deseq2 <- as.data.frame(resOrdered) %>% na.omit()
# Add DEG type (up / down / stable)
DEG_deseq2 <- DEG_deseq2 %>%
  mutate(Type = if_else(pvalue > 0.05, "stable",
                        if_else(abs(log2FoldChange) < 1, "stable",
                                if_else(log2FoldChange >= 1, "up", "down")))) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  rownames_to_column("Symbol")
head(DEG_deseq2)

##-----Step 4. Volcano plot of DEGs-----
library(ggplot2)
ggplot(DEG_deseq2, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(size = 2.5, alpha = 0.8, aes(color = Type)) +
  scale_color_manual(values = c("#00468B", "gray", "#E64B35")) +
  ylim(0, 15) + xlim(-5, 5) +
  labs(x = "Log2 Fold Change", y = "-log10 adjusted p-value") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = 'black') +
  geom_vline(xintercept = c(-1, 1), linetype = 2, color = 'black') +
  theme_bw()



##----Step 5. GSEA with fgsea----
library(fgsea)
library(msigdbr)

# Rank genes by log2FC
ranks <- DEG_deseq2$log2FoldChange
names(ranks) <- DEG_deseq2$Symbol
ranks <- sort(ranks, decreasing = TRUE)

# Use KEGG/Reactome pathways
genesets_k <- read.gmt("D:/R/R-4.3.2/library/scMetabolism/data//KEGG_metabolism_nc.gmt") 
geneset_list <- split(genesets_k$gene, genesets_k$term)

# Run fgsea
library(BiocParallel)
register(SerialParam())  # avoid parallel issues
fgseaRes <- fgsea(pathways = geneset_list, stats = ranks)
fgseaRes <- fgseaRes[order(fgseaRes$padj), ]
head(fgseaRes)

# Plot top enriched pathway
plotEnrichment(geneset_list[[fgseaRes$pathway[1]]], ranks) +
  labs(title = fgseaRes$pathway[1])

##----Step 6. Barplot of enriched pathways----
library(dplyr)
top_pathways <- fgseaRes %>%
  filter(pval < 0.05) %>%
  arrange(desc(NES)) %>%
  head(15) %>%
  mutate(group = ifelse(NES > 0, "OAT_pos", "OAT_neg"))

ggplot(top_pathways, aes(x = reorder(pathway, NES), y = NES, fill = group)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("OAT_pos" = "firebrick", "OAT_neg" = "steelblue")) +
  labs(x = "Pathway", y = "Normalized Enrichment Score (NES)", fill = "Enriched in") +
  theme_bw()



##----Step 7. ssGSEA validation----
library(GSVA)
pb_expr_mat <- as.matrix(pb_expr)
# Run ssGSEA
ssgsea_res <- gsva(pb_expr_mat, geneset_list, method = "ssgsea", kcdf = "Gaussian", verbose = FALSE)
# Format results
ssgsea_df <- as.data.frame(t(ssgsea_res)) %>%
  rownames_to_column("sample_group") %>%
  left_join(meta %>% rownames_to_column("sample_group"), by = "sample_group")
# Example: boxplot for one pathway
ggplot(ssgsea_df, aes(x = OAT_status, y = `Arginine and proline metabolism`, fill = OAT_status)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2) +
  labs(title = "ssGSEA: Arginine and proline metabolism") +
  scale_fill_manual(values = c("#619CFF", "#F8766D")) +
  theme_bw()

##-----Step 8. Statistical testing & heatmap-----
# Mann-Whitney test for each pathway
library(purrr)
diff_test <- function(pathway_name) {
  group1 <- ssgsea_df %>% filter(OAT_status == "OAT_pos") %>% pull(pathway_name)
  group2 <- ssgsea_df %>% filter(OAT_status == "OAT_neg") %>% pull(pathway_name)
  wilcox.test(group1, group2)$p.value
}

pvals <- map_dbl(colnames(ssgsea_df)[2:(ncol(ssgsea_df)-2)], diff_test)
names(pvals) <- colnames(ssgsea_df)[2:(ncol(ssgsea_df)-2)]
pval_df <- data.frame(pathway = names(pvals), pvalue = pvals) %>%
  arrange(pvalue)
head(pval_df)

# Select top 10 significant pathways
top_pathways <- pval_df %>% filter(pvalue < 0.1) %>% head(10) %>% pull(pathway)
# Heatmap visualization
ss_top <- ssgsea_df[, top_pathways]
rownames(ss_top) <- ssgsea_df$sample_group
anno <- data.frame(Group = ssgsea_df$OAT_status)
rownames(anno) <- ssgsea_df$sample_group

pheatmap(
  t(scale(t(ss_top))),
  annotation_row = anno,
  color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  cluster_rows = F, cluster_cols = T,
  fontsize_row = 10, border_color = NA
)


##----Step 9. fgsea styled dotplot----
library(ggplot2)
library(ggrepel)
library(dplyr)

# ==== 数据预处理 ====
fgseaRes <- fgseaRes[order(fgseaRes$NES, decreasing = TRUE), ]
fgseaRes$xlab <- 1:nrow(fgseaRes)  # 添加排名位置
fgseaRes$logp <- -log10(fgseaRes$pval)  # 计算 -log10(p-value)
fgseaRes$setSize <- fgseaRes$size / 2   # 控制点大小（可自行调整比例）

# ==== 指定要标注的重点通路（可修改） ====
highlight_pathways <- c("Arginine and proline metabolism", "Glycosaminoglycan degradation", 
                        "Lipoic acid metabolism", "Linoleic acid metabolism")
fgseaRes$highlight <- ifelse(fgseaRes$pathway %in% highlight_pathways, "highlight", "none")

# ==== 颜色方案 ====
fill_color <- "#1f77b4"  # 蓝色系

# ==== 主图 ====
p <- ggplot(fgseaRes, aes(x = xlab, y = NES)) +
  geom_point(aes(size = setSize, fill = logp), shape = 21, color = "black", alpha = 0.9) +
  scale_fill_gradient(low = "lightgrey", high = fill_color, name = "-log10(p-value)") +
  scale_size_continuous(range = c(2, 10), name = "Pathway size") +
  theme_classic(base_size = 15) +
  xlab("Ranked pathways") +
  ylab("Normalized Enrichment Score (NES)") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "right"
  )

# ==== 添加标注 ====
label_data <- fgseaRes %>% filter(highlight == "highlight")
p <- p + geom_text_repel(
  data = label_data,
  aes(label = pathway),
  size = 4,
  color = "black",
  nudge_y = 0.5,
  segment.color = "grey40",
  segment.size = 0.3
)

# ==== 保存并展示 ====
ggsave("fgsea_styled_plot.pdf", p, width = 8, height = 6)
p

