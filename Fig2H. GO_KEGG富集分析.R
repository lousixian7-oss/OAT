library(clusterProfiler)#GO富集分析、KEGG通路富集分析
library(org.Hs.eg.db)#基因注释数据库
library(enrichplot)
library(ggplot2)
library(GOplot)
library(dplyr)

#设置工作路径
setwd("工作路径")
#读入数
genes_df <- read.csv("核心基因集.csv")  

# 将基因符号转换为ENTREZ ID
entrezIDs <- bitr(genes_df$gene, fromType = "SYMBOL", 
               toType = c("ENTREZID", "SYMBOL"),
               OrgDb = org.Hs.eg.db) # `OrgDb`参数指定使用的人类基因组注释数据库


#使用entrezIDs 
gene<- entrezIDs$ENTREZID
##GO富集分析
go<- enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, qvalueCutoff = 0.05,ont="all",readable =T)
write.table(go,file="GO.txt",sep="\t",quote=F,row.names = F) #

##可视化
##条形图
pdf(file="GO-柱状图.pdf",width = 8,height = 10)
##showCategory改变自己想要展示的条目数量
barplot(go, drop = TRUE, showCategory =6,label_format=100,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()

##气泡图
pdf(file="GO-气泡图.pdf",width = 6,height = 6)
dotplot(go,showCategory = 3,label_format=100,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
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


# 过滤 enrichResult 类型的对象，保留 enrichResult 结构
go_selected <- go

# 重点标记的通路（**只对 BP 生效**）

# # **创建 y 轴标签颜色映射**

library(stringr)
go_selected@result <- go@result %>% filter(Description %in% selected_paths)
# 设置换行宽度（你可以根据图宽微调，例如 width = 40）
go_selected@result$Description <- str_wrap(go_selected@result$Description, width = 40)
go_selected@result$Description <- str_to_sentence(go_selected@result$Description)
# 生成 GO 气泡图

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

# p <- dotplot(go_selected, showCategory = length(selected_paths), label_format = 100, split="ONTOLOGY") +
#   facet_grid(ONTOLOGY~., scale='free')
  # theme_minimal() +  # 让图表更简洁
  # theme(
  #   axis.text.y = element_blank(),  # **隐藏默认的 y 轴标签**
  #   panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # **添加边框**
  # ) +
  # geom_text(data = go_selected@result,
  #           aes(x = min(go_selected@result$Count) - 1,  # **确保文本出现在 y 轴**
  #               y = reorder(Description, Count),
  #               label = Description,
  #           hjust = 1, size = 4, show.legend = FALSE))  # **确保颜色正确解析**
# 
# print(p)
# dev.off()



#kegg分析
kk <- enrichKEGG(gene = gene,keyType = "kegg",organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05, pAdjustMethod = "fdr")   
write.table(kk,file="KEGG.txt",sep="\t",quote=F,row.names = F)                         
# 加载包
library(ReactomePA)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ReactomePA")

# Reactome富集分析
reactome_result <- enrichPathway(
  gene         = gene,         # 注意是Entrez ID
  organism     = "human",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable     = TRUE
)

##可视化
##条形图
pdf(file="KEGG-柱状图.pdf",width = 8,height = 6)
barplot(kk, drop = TRUE, showCategory = 8,label_format=100)
dev.off()

##气泡图
pdf(file="KEGG-气泡图.pdf",width = 6,height = 4)
dotplot(kk, showCategory = 8,label_format=100)
dev.off()


# 加载必要的 R 包
library(clusterProfiler)
library(ggplot2)
library(dplyr)

# 设定需要展示的 KEGG 通路
selected_kegg <- c(
  "Arginine and proline metabolism"
  
)


# 过滤 enrichResult 类型的对象，保留 enrichResult 结构
kk_selected <- kk
kk_selected@result <- kk@result %>% filter(Description %in% selected_kegg)
kk_selected@result$Description <- str_wrap(kk_selected@result$Description, width = 40)
kk_selected@result$Description <- str_to_sentence(kk_selected@result$Description)
# # 重点标记的 KEGG 通路（**高亮**）
# highlighted_kegg <- c("Inflammatory bowel disease",
#                       "IL-17 signaling pathway",
#                       "Osteoclast differentiation")
# 
# # **创建 y 轴标签颜色映射**
# kk_selected@result <- kk_selected@result %>%
#   mutate(label_color = ifelse(Description %in% highlighted_kegg, "#CC6666", "black"))

# 生成 KEGG 气泡图
pdf(file = "KEGG-气泡图.pdf", width = 6, height = 3)

p <- dotplot(kk_selected, showCategory = length(selected_kegg), label_format = 100)+ 
  theme_minimal()+  # 让图表更简洁
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
 
  #   axis.text.y = element_blank(),  # **隐藏默认的 y 轴标签**
  #  panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # **确保边框显示**
  # )
  # geom_text(data = kk_selected@result, 
  #           aes(x = min(kk_selected@result$Count) - 1,  # **确保文本出现在 y 轴**
  #               y = reorder(Description, Count), 
  #               label = Description, 
  #               color = ifelse(Description %in% highlighted_kegg, "#CC6666", "black")), 
  #           hjust = 1, size = 4, show.legend = FALSE) +  
  # scale_color_manual(values = c("#CC6666" = "#CC6666", "black" = "black"))  # **确保颜色正确解析**

print(p)
dev.off()




library("enrichplot")
geneList<-genes_df$gene
## categorySize can be scaled by 'pvalue' or 'geneNum'
pdf("plot.pdf", width = 12, height = 10)
cnetplot(go, foldChange=geneList)
dev.off()
cnetplot(go, categorySize="pvalue", foldChange=geneList)
pdf("plot.pdf", width = 12, height = 10)
cnetplot(kk, foldChange=geneList, circular = TRUE, colorEdge = TRUE)


library(tidyverse)
 devtools::install_github("davidsjoberg/ggsankey")
library(ggsankey)
library(ggplot2)
install.packages("cols4all")
library(cols4all)
df <- read.csv("df.csv", header = T)
df1 <- df[,-1]
colnames(df1)
df1_trans <- df1 %>%make_long(gene, pathway)
df1_trans$node <- factor(df1_trans$node,levels = c(df1$pathway %>% unique()%>% rev(),
                                                   df1$gene %>% unique() %>% rev()))

colnames(df1_trans)
ggplot(df1_trans, aes(x = x,
                      next_x= next_x,
                      node= node,
                      next_node= next_node,
                      fill= node,
                      label= node)) +
  geom_sankey(flow.fill="#DFDFDF",
              flow.color="grey60",
              node.fill=dittoColors()[1:44],
              width=0.15) + 
  geom_sankey_text(size = 3,
                   color= "black",
                   hjust=1) + 
  theme_void()


df <- read.table("KEGG.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE)
library(tidyr)
df <- df %>%
  separate_rows(geneID, sep = "/")  # 按 "/" 拆分基因列
ggplot(df, aes(axis1 = geneID, axis2 = Description)) +
  geom_alluvium(aes(fill = Description), width = 1/12) +
  geom_stratum(width = 1/12, fill = "gray") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Gene", "KEGG Pathway"), expand = c(0.2, 0.2)) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "KEGG Pathway Sankey Diagram", x = "", y = "")
df <- read_excel("1.xlsx")  # 请替换文件路径
df <- df %>%
  separate_rows(geneID, sep = "/")  # 按 "/" 拆分基因列
df<-df[1:123,]

ggplot(df, aes(axis1 = geneID, axis2 = Description)) +
  geom_alluvium(aes(fill = Description), width = 1/6) +  # **填充颜色由 Description 决定**
  geom_stratum(aes(fill = Description), width = 1/6) +   # **修改灰色，使其与流颜色一致**
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4) +
  scale_x_discrete(limits = c("Gene", "GO Term"), expand = c(0.3, 0.3)) +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "GO Enrichment Sankey Diagram", x = "", y = "")



library(ggplot2)
library(ggtree)
library(tidyverse)
library(clusterProfiler)
library(DOSE)
library(enrichplot)
library(readxl)
library(readxl)
library(ggplot2)
library(ggalluvial)
library(dplyr)
library(readr)

# 读取数据（请替换成你的文件路径）
df <- read_excel("2.xlsx")
df <-df[1:20,]
# 选择需要的列并整理数据
df_plot <- df %>%
  select(category, subcategory, Description, Count) %>%
  rename(Category = category, Subcategory = subcategory, Pathway = Description) %>%
  arrange(Category, Subcategory)

my.color <- c("#FDDED7","#A8D3A0","#FFDD8E",'#C1E0DB',"#C9A1CA" ,'#70CDBE',"#7AC3DF","#F6AA61","#AC99D2")

ggplot(df_plot, aes(axis1 = Category, axis2 = Subcategory, axis3 = Pathway, y = Count)) +
  geom_flow(aes(fill = Category), curve_type = "cubic", width = 0.3) +
  geom_stratum(width = 0.3, fill = "gray80", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4) +
  scale_x_discrete(limits = c("Category", "Subcategory", "Pathway")) +
  labs(title = "KEGG Pathway Classification", x = "", y = "Gene Count") +
  theme_minimal()

ggplot(df_plot, aes(axis1 = Category, axis2 = Subcategory, axis3 = Pathway, y = Count)) +
  # **绘制流**
  geom_flow(aes(fill = Category), curve_type = "cubic", width = 0.3) +
  
  # **单独修改第一列（Category）色块的长度**
  geom_stratum(aes(fill = Category, width = ifelse(after_stat(x) == 1, 0.3,ifelse(after_stat(x) == 2, 0.3, 0.3))), 
               color = "black") +  
  scale_fill_manual(values = my.color) +
  # **添加文本**
  #geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 5) +
  
  # **X 轴顺序**
  scale_x_discrete(limits = c("Category", "Subcategory", "Pathway")) +
  
  # **标签**
  labs(title = "KEGG Pathway Classification", x = "", y = "Gene Count") +
  
  # **美化主题**
  theme_minimal()




