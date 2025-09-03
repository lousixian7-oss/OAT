library(dplyr)
library(data.table)
library(coloc)
library(vroom)
library(locuscomparer)
library(magrittr)

# UKB_PPP
IVs <- fread("./OAT.txt",data.table = F)
IVs <- IVs %>%
  filter(Chr == 10 &
           BP >= 126085872-1000000 &
           BP <= 126107545+1000000 )
colnames(IVs)
data <- IVs %>% dplyr::select("SNP", "Chr", "BP", "A1",
                              "A2", "Freq", "b", "SE",
                              "p")

# Rename columns
colnames(data) <- c('SNP', 'chr', "pos", 'A1', 'A2', "eaf", "beta", "se", "p")
data$eaf <- as.numeric(data$eaf)
# Prepare data required for colocalization
data$MAF <- ifelse(data$eaf<0.5,data$eaf,1-data$eaf)
data$n <- "31199"
data$MAF <- as.numeric(data$MAF)
data$se<-as.numeric(data$se)
data$varbeta <- data$se^2
data$beta<-as.numeric(data$beta)
data$z <- data$beta / data$se
data <- subset(data, !duplicated(SNP))
lead <- data %>% dplyr::arrange(p)
leadPos <- lead$pos[1]
leadPos <- as.numeric(leadPos)
lead$pos<-as.numeric(lead$pos)
QTLdata <- lead %>% filter(pos > leadPos - 500000, pos < leadPos + 500000) %>% na.omit()

### Outcome data
GWASdata <- read_outcome_data(snps = QTLdata$SNP, 
                              filename="finngen_R11_K11_GINGIVITIS_PERIODONTAL.gz",
                              sep = "\t",
                              beta_col = "beta",
                              se_col = "sebeta",
                              effect_allele_col = "alt",
                              other_allele_col = "ref",
                              snp_col = "rsids",
                              chr_col = "#chrom",
                              pos_col = "pos",
                              pval_col = "pval",
                              eaf_col = "af_alt")
column_types <- lapply(QTLdata, class)
column_types
GWASdata$chr.outcome<- as.numeric(GWASdata$chr.outcome)
GWASdata$pos.outcome<- as.numeric(GWASdata$pos.outcome)
GWASdata <- GWASdata %>% dplyr::select("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome",
                                       "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome","pval.outcome")

# Rename columns
colnames(GWASdata) <- c('SNP', 'chr', "pos", 'A1', 'A2', "eaf", "beta", "se", "p")
GWASdata$n <- "34206"
GWASdata$n<-as.numeric(GWASdata$n)
GWASdata$MAF <- ifelse(GWASdata$eaf<0.5,GWASdata$eaf,1-GWASdata$eaf)
GWASdata$MAF <- as.numeric(GWASdata$MAF)
GWASdata$se<-as.numeric(GWASdata$se)
GWASdata$varbeta <- GWASdata$se^2
GWASdata$beta<-as.numeric(GWASdata$beta)
GWASdata$z <- GWASdata$beta/GWASdata$se

### Process data
sameSNP <- intersect(QTLdata$SNP,GWASdata$SNP)
QTLdata <- QTLdata[QTLdata$SNP %in% sameSNP, ] %>% dplyr::arrange(SNP) %>% na.omit()
GWASdata1 <- GWASdata[GWASdata$SNP %in% sameSNP, ] %>% dplyr::arrange(SNP) %>% na.omit()
QTLdata$n<-as.numeric(QTLdata$n)
QTLdata$p<-as.numeric(QTLdata$p)
# Colocalization analysis
# GWASdata$MAF <- as.numeric(GWASdata$MAF)
coloc_data <- list(dataset1=list(snp=QTLdata$SNP,beta=QTLdata$beta,varbeta=QTLdata$varbeta,
                                 N=QTLdata$n,MAF=QTLdata$MAF,z = QTLdata$z,
                                 pvalues=QTLdata$p,type="quant"), 
                   dataset2=list(snp=GWASdata1$SNP,beta=GWASdata1$beta,varbeta=GWASdata1$varbeta,
                                 N=GWASdata1$samplesize,MAF=GWASdata1$MAF,z = GWASdata1$z,
                                 pvalues=GWASdata1$p,type="cc"))

result <- coloc.abf(dataset1=coloc_data$dataset1, dataset2=coloc_data$dataset2)
result1<- result$summary
result2<- result$results %>% dplyr::arrange(desc(SNP.PP.H4))

write.csv(result1,"coloc_data_result1.csv",row.names = F)
write.csv(result2,"coloc_data_result2.csv",row.names = F)
