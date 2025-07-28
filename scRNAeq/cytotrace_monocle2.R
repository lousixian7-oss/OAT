obj.fib <- readRDS("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/obj.fib.rds")
# ==== 脚本说明 ====
# 目的：CytoTRACE辅助monocle确定分化起点 
# 包括：加载数据 → 运行模型 → 可视化 → 输出结果
# 数据：obj.fib cds_DGT
# ---- monocle2分析 ----
library(Seurat)
library(monocle)
sce_fib_new_monocle_fun <- seurat_to_monocle(obj.fib, assay = "RNA", slot = "counts")
#Estimate size factors and dispersions
sce_fib_new_monocle <- estimateSizeFactors(sce_fib_new_monocle_fun)
sce_fib_new_monocle <- estimateDispersions(sce_fib_new_monocle)
#质量控制
sce_fib_new_monocle <- detectGenes(sce_fib_new_monocle, min_expr = 0.1)#至少10%的表达，这一步完成之后，会出现num_cells_expressed这一列
print(head(fData(sce_fib_new_monocle)))
print(head(pData(sce_fib_new_monocle)))
expressed_genes <- row.names(subset(fData(sce_fib_new_monocle), num_cells_expressed >= 10))
pData(sce_fib_new_monocle)$Total_mRNAs <- Matrix::colSums(exprs(sce_fib_new_monocle))#将UMI加入cds
sce_fib_new_monocle <- sce_fib_new_monocle[,pData(sce_fib_new_monocle)$Total_mRNAs < 1e6]#阈值按照自己实际情况
#differentialGeneTest
cds_DGT <- sce_fib_new_monocle
diff_test_res <- differentialGeneTest(cds_DGT,fullModelFormulaStr = "~annotation3")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
# ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval)][1:2000]#或者选择前多少的基因
cds_DGT <- setOrderingFilter(cds_DGT, ordering_genes)
plot_ordering_genes(cds_DGT)
getwd()

#接下来我们以第三种方法进行
#降维,排序
cds_DGT <- reduceDimension(cds_DGT, max_components = 2,reduction_method = 'DDRTree')
cds_DGT <- orderCells(cds_DGT)
cds_DGT <- orderCells(cds_DGT, root_state = NULL, num_paths = NULL, reverse = T) 
saveRDS(cds_DGT,"D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H10. Final/Figure2/cds_DGT_monocle.rds")
#初步可视化
plot_cell_trajectory(cds_DGT, cell_size = 2.2, color_by = "celltype") +
  facet_wrap(~celltype, nrow = 2)
plot_cell_trajectory(cds_DGT, color_by = "annotation3")
plot_cell_trajectory(cds_DGT, color_by = "State")
plot_cell_trajectory(cds_DGT, color_by = "Pseudotime")
plot_cell_trajectory(cds_DGT, color_by = "diseaseStatus")

#基因表达拟时图
plot_cell_trajectory(cds_DGT, markers = "OAT",use_color_gradient=T,cell_size = 1,cell_link_size = 1.5)
#各个分组拟时图
plot_cell_trajectory(cds_DGT, color_by = "annotation3") + facet_wrap(~diseaseStatus, nrow = 2)
cds_DGT <- orderCells(cds_DGT, root_state = 1)  # 假设State 3是最stem-like的

pdf("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H10. Final/Figure2/Pserdotime_umap.pdf", width = 6, height = 6)
plot_cell_trajectory(cds_DGT, color_by = "Pseudotime")
dev.off()
pdf("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H10. Final/Figure2/Pserdotime_OAT_umap.pdf", width = 6, height = 6)
plot_cell_trajectory(cds_DGT, markers = "OAT",use_color_gradient=T,cell_size = 1,cell_link_size = 1.5)
dev.off()
pdf("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H10. Final/Figure2/Pserdotime_OAT_time.pdf", width = 5, height = 4)
plot_genes_in_pseudotime(cds_DGT["OAT", ],color_by = "annotation3") 
dev.off()

# 找出打分最高的 cell 所属 state
pData(cds_DGT)$cytotrace <- results$CytoTRACE[rownames(pData(cds_DGT))]
# 可视化染色
plot_cell_trajectory(cds_DGT, color_by = "cytotrace") + scale_color_viridis_c()
# 查看哪个 State 中 CytoTRACE 平均值最高
state_scores <- aggregate(pData(cds_DGT)$cytotrace, by = list(State = pData(cds_DGT)$State), mean)
root_state <- state_scores$State[which.max(state_scores$x)]
root_state
# 明确指定轨迹方向起点
cds_DGT <- orderCells(cds_DGT, root_state = root_state)

library(RColorBrewer)
cols <- rev(brewer.pal(11, "Spectral"))  # 蓝→红渐变
plot_cell_trajectory(cds_DGT, color_by = "cytotrace", size = 1, show_backbone = TRUE) +
  scale_colour_gradientn(colours = cols)
pdf("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H10. Final/Figure2/cytotrace_Spectral.pdf", width = 6, height = 6)
plot_cell_trajectory(cds_DGT, color_by = "cytotrace", size = 1, show_backbone = TRUE) +
  scale_colour_gradientn(colours = cols)
dev.off()

# ==== 脚本说明 ====
# 目的：CytoTRACE 分析流程,补充分析，辅助确定起点
# 数据：obj.fib
# ---- CytoTRACE分析 ----
library(CytoTRACE)
library(Seurat)
library(RColorBrewer)
mat <- as.matrix(obj.fib@assays$RNA@counts)
results <- CytoTRACE(mat = mat)# 得到list对象      
obj.fib$CytoTRACE <- results$CytoTRACE #将得分加入obj.fib
NM <- SCP::FeatureDimPlot(obj.fib, 
                          features = "CytoTRACE", 
                          reduction = "umap", 
                          pt.size = 1.2,
                          theme_use = ggplot2::theme_classic, 
                          theme_args = list(base_size = 16)
)
pdf("D:/scRNAseq/06.Proj_processing/Sen_yazhouyan/H10. Final/Figure2/cytotrace_umap.pdf", width = 5, height = 5)
SCP::FeatureDimPlot(obj.fib, 
                    features = "CytoTRACE", 
                    reduction = "umap", 
                    pt.size = 1.2,
                    theme_use = ggplot2::theme_classic, 
                    theme_args = list(base_size = 16)
)
dev.off()


# ==== 脚本说明 ====
# 目的：CytoTRACE2 分析流程,补充分析，辅助确定起点
# 数据：obj.fib
# ---- CytoTRACE2分析 ----
# Step1：运行CytoTRACE2模型
library(CytoTRACE2) 
sce
save.dir <- ""
cytotrace2_result <- cytotrace2(
  obj.fib,
  species = "human",           # 如果是小鼠就写 "mouse"
  slot_type = "counts",
  is_seurat = TRUE,
  batch_size = 10000,
  smooth_batch_size = 1000,
  parallelize_models = TRUE,
  parallelize_smoothing = TRUE,
  seed = 1234
)
# Step2：将打分结果加入obj.fib
obj.fib$CytoTRACE2_Score     <- cytotrace2_result$CytoTRACE2_Score
obj.fib$CytoTRACE2_Potency   <- cytotrace2_result$CytoTRACE2_Potency
obj.fib$CytoTRACE2_Relative  <- cytotrace2_result$CytoTRACE2_Relative
# Step 3: 构建分组注释 annotation 信息
annotation_df <- data.frame(celltype = obj.fib$annotation3)
rownames(annotation_df) <- colnames(obj.fib)
plots <- plotData(
  cytotrace2_result = cytotrace2_result,
  expression_data = NULL,
  annotation = annotation_df,
  is_seurat = TRUE,
  slot_type = "counts",
  pc_dims = 50  # 维数如有PCA可保留为 50
)
# Step 4: 可视化 CytoTRACE2 结果
plots$CytoTRACE2_UMAP
plots$CytoTRACE2_Potency_UMAP
plots$CytoTRACE2_Relative_UMAP
plots$CytoTRACE2_Boxplot_byPheno
plots$Phenotype_UMAP
