library(locuscomparer)
library(data.table)
gwas_fn = fread('periodontal.out.csv')
head(gwas_fn)
# cleaned_data <- na.omit(gwas_fn)
gwas_fn <- gwas_fn[,c('SNP','pval')]
colnames(gwas_fn) <- c('rsid','pval')
fwrite(gwas_fn,file="gwas.txt",sep="\t")
gwas_fn=fread("gwas.txt")

workdir <- getwd()  # Get the path of the current working directory
# 
# workdir <- paste0(getwd(), "/exposure")  
# eqtl_files <- list.files(path = workdir, pattern = ".txt", full.names = TRUE) 
eqtl="OAT.txt"
eqtl_fn  <- fread(eqtl)
leadsnppos="126428674"  ###### Position of the top SNP for colocalization
leadsnppos <- as.numeric(leadsnppos)
eqtl_fn <- eqtl_fn[eqtl_fn$BP> (leadsnppos - 20000) & eqtl_fn$BP < (leadsnppos + 20000), ]
# Select required columns and rename
eqtl_fn <- eqtl_fn[,c('SNP','p')]
colnames(eqtl_fn) <- c('rsid','pval')
# filename <- unlist(strsplit(basename(eqtl), "_"))[3]
# Call the locuscompare function and generate PDF plots
# pdf(file.path(workdir, paste0(filename, "OAT.pdf")), width = 6, height = 8)
pdf(file.path(workdir, paste0("OAT.pdf")), width = 6, height = 6)
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Periodontitis', title2 = 'eQTL')
dev.off() 

eqtl="DPP8.txt"
eqtl_fn  <- fread(eqtl)
leadsnppos="66196997"  ###### Position of the top SNP for colocalization 65734801–65810042
leadsnppos <- as.numeric(leadsnppos)
eqtl_fn <- eqtl_fn[eqtl_fn$BP> (65734801) & eqtl_fn$BP < (65810042), ]
# Select required columns and rename
eqtl_fn <- eqtl_fn[,c('SNP','p')]
colnames(eqtl_fn) <- c('rsid','pval')
# filename <- unlist(strsplit(basename(eqtl), "_"))[3]
# Call the locuscompare function and generate PDF plots
# pdf(file.path(workdir, paste0(filename, "OAT.pdf")), width = 6, height = 8)
pdf(file.path(workdir, paste0("DPP8.pdf")), width = 6, height = 6)
locuscompare(in_fn1 = gwas_fn, in_fn2 = eqtl_fn, title = 'Periodontitis', title2 = 'eQTL')
dev.off()   
