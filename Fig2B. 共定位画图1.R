
library(locuscomparer)
library(data.table)
gwas_fn = fread('periodontal.out.csv')######需要改
head(gwas_fn)
#cleaned_data <- na.omit(gwas_fn)
gwas_fn <- gwas_fn[,c('SNP','pval')]######需要改
colnames(gwas_fn) <- c('rsid','pval')
fwrite(gwas_fn,file="gwas.txt",sep="\t")
gwas_fn=fread("gwas.txt")

workdir <- getwd()  # 获取当前工作目录的路径
# 
# workdir <- paste0(getwd(), "/exposure")  
# eqtl_files <- list.files(path = workdir, pattern = ".txt", full.names = TRUE) 
eqtl="OAT.txt"
eqtl_fn  <- fread(eqtl)######需要改
leadsnppos="126428674"######共定位topsnp的位置
leadsnppos <- as.numeric(leadsnppos)
eqtl_fn <- eqtl_fn[eqtl_fn$BP> (leadsnppos - 20000) & eqtl_fn$BP < (leadsnppos + 20000), ]
# 选取需要的列，重命名列名
eqtl_fn <- eqtl_fn[,c('SNP','p')]
colnames(eqtl_fn) <- c('rsid','pval')
# filename <- unlist(strsplit(basename(eqtl), "_"))[3]
# 调用 locuscompare 函数并生成PDF图像
# pdf(file.path(workdir, paste0(filename, "OAT.pdf")), width = 6, height = 8)
pdf(file.path(workdir, paste0("OAT.pdf")), width = 6, height = 6)
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Periodontitis', title2 = 'eQTL')
dev.off() 

eqtl="DPP8.txt"
eqtl_fn  <- fread(eqtl)######需要改
leadsnppos="66196997"######共定位topsnp的位置65734801	65810042
leadsnppos <- as.numeric(leadsnppos)
eqtl_fn <- eqtl_fn[eqtl_fn$BP> (65734801) & eqtl_fn$BP < (65810042), ]
# 选取需要的列，重命名列名
eqtl_fn <- eqtl_fn[,c('SNP','p')]
colnames(eqtl_fn) <- c('rsid','pval')
# filename <- unlist(strsplit(basename(eqtl), "_"))[3]
# 调用 locuscompare 函数并生成PDF图像
# pdf(file.path(workdir, paste0(filename, "OAT.pdf")), width = 6, height = 8)
pdf(file.path(workdir, paste0("DPP8.pdf")), width = 6, height = 6)
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Periodontitis', title2 = 'eQTL')
dev.off()   

