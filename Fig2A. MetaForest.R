#meta 森林图
#能够调整的地方都是有注释的 根据含义来调整 看看具体出来的图有什么变化
#这样可以完全搞懂修改的地方具体是做什么的
library(forestploter)
library(readr)


data = read_csv("res.csv")
#这里的15可以改变森林图那些竖线的宽度
data$' ' <- paste(rep(" ", 15), collapse = " ") 

#增加置信区间显示列
data$'OR(95% CI)' <- ifelse(is.na(data$or), "",
                           sprintf("%.2f (%.2f to %.2f)",
                                   data$or,
                                   data$or_lci95, data$or_uci95))

#这里的%.3f是为了保留3位小数
#data$P = ifelse(is.na(data$P), "", sprintf("%.3f", data$P))
data$P = ifelse(is.na(data$P), "", data$P)
#这里的表头ifesle 是为了让一些NA值显示为空
#如果对应的表头修改了，需要同时改这里的表头名，删掉或者注释都行
#新增了表头，需要这个显示为空的功能，按照格式添加就行
data$Exposure = ifelse(is.na(data$Exposure), "", data$Exposure)
data$`No.of SNP` = ifelse(is.na(data$`No.of SNP`), "", data$`No.of SNP`)
#data$`MR-Egger  intercept P value` = ifelse(is.na(data$`MR-Egger  intercept P value`), "", data$`MR-Egger  intercept P value`)

#森林图主题  可以调各种样式
tm <- forest_theme(base_size = 10, # 基础大小
                   base_family = "serif", #字体可以执行 windowsFonts()来查看
                   # 可信区间点的形状，线型、颜色、宽度
                   ci_pch = 16,
                   ci_col = "#4575b4", # #762a83
                   ci_lty = 1,
                   ci_lwd = 1.5,
                   ci_Theight = 0.2, # 可信区间两端加短竖线

                   # 参考线宽度、形状、颜色
                   refline_lwd = 1,
                   refline_lty = "dashed",
                   refline_col = "grey20",

                   # 汇总菱形的填充色和边框色
                   summary_fill = "#4575b4",
                   summary_col = "#4575b4",

                   # 脚注大小、字体、颜色
                   footnote_cex = 0.6,
                   footnote_fontface = "italic",
                   footnote_col = "blue",
                   )

plot <- forest(
	#需要的列 默认不变 调这里可以显示自己需要的数据 具体看右侧的环境变量的data 展开填入对应的列号
	data[, c(1,2,3,4,9,8)], 
	#or值 以及对应的置信区间
               est = data$or, 
               lower = data$or_lci95,
               upper = data$or_uci95,
               ci_column = 5,  #前面定义了一个空的列 用来画森林图的 按需修改
               ref_line = 1,  #对应的参考线，如果是beta值 请改成0
	#区间大小 这里控制区间的显示范围
                xlim = c(0.80, 1.30), 
	#刻度的粒度 这里控制刻度的显示 需要删掉ticks_at前面的#注释符号
               ticks_at = c(0.80,0.90, 1, 1.10,1.20,1.30), 
               theme = tm,
               )
#画图 这里可以改p1.pdf 输出的文件名
#也可以更改width宽和height高 调整图的宽高
pdf("p1.pdf", family = "GB1", width = 10, height = 8)
plot
dev.off()

