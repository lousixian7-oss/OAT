# Set your own working directory

# - Load required packages, install them first if not already installed
library(TwoSampleMR)
library(data.table)
library(tidyverse)
library(readxl)
library(writexl)
library(mr.raps)
library(ieugwasr)
library(dplyr)
library(MendelianRandomization)


# Exposure: put the data in your working directory
FileNames <-list.files(paste0(getwd()),pattern=".csv")
exp_dat_ids <- FileNames
exps <- FileNames


# Outcome
out<-read.csv("periodontal.out.csv",header = T) # Read local outcome file (replace with your own outcome file name)
out$trait <- 'periodontal'  # Change to the trait name used in your outcome data
outcomeid <- out
rm(out)
head(outcomeid)
outcomeid<-select(outcomeid,-c(X.1,X,nearest_genes,af_alt_controls,af_alt_cases,mlogp))
colnames(outcomeid)<-c("chr","pos","other_allele","effect_allele","SNP","pval","beta","se","eaf","ncase","ncontrol","trait")#
outcomeid$samplesize=outcomeid$ncase+outcomeid$ncontrol
head(outcomeid)
#
dir.create(path = "mendelian test") # Create a new folder in the working directory to save results

####### Loop code below ######################
for (qaq in 1484:length(exp_dat_ids)) { # 
  exp_dat_id <- exp_dat_ids[qaq]
  exp <- exps[qaq]
  
  d3<- try(fread(paste0(getwd(),"/",FileNames[qaq]),fill=TRUE),silent = F)

  
  if (nrow(d3) > 0) {
    d3<- as.data.frame(d3)

    
    d3<-format_data(
      d3,
      type = "exposure",
      header = TRUE,
      snp_col = "SNP",
      beta_col = "beta.exposure",
      se_col = "se.exposure",
      eaf_col = "eaf.exposure",
      effect_allele_col = "effect_allele.exposure",
      other_allele_col = "other_allele.exposure",
      pval_col = "pval.exposure",
      samplesize_col = "samplesize.exposure",
      gene_col = "gene",
      id_col = "id",
      chr_col = "chr",
      pos_col = "pos",
    )
    exp_data<-d3
    # exp_data <- clump_data(d3,clump_kb = 500,clump_r2 = 0.01)
    
    library(ieugwasr)     
    
    # d4<- ld_clump(
    #   # dat = X1,
    #   clump_kb = 10000,
    #   clump_r2 = 0.2,    # LD threshold
    #   pop = "EUR",
    #   dplyr::tibble(rsid=d3$SNP, pval=d3$pval.exposure, id=d3$id.exposure),
    #   # Path to plink.exe
    #   plink_bin = "E:/medelian randomization/MR analysis/local plink/plink_win64_20230116/plink.exe",
    #   # European reference panel path
    #   bfile = "E:/medelian randomization/MR analysis/local plink/EUR/EUR"
    # )
    # 
    # exp_data<-subset(d3,SNP %in% d4$rsid)
    
    
    
    if(length(exp_data[,1])>0){
      outcome_dat<-merge(exp_data,outcomeid,by.x = "SNP",by.y = "SNP")
      if(length(outcome_dat[,1])>0){  
        write.csv(outcome_dat,file = "d.csv")
        out_data <- read_outcome_data(
          snps = exp_data$SNP, 
          filename = "d.csv",
          sep = ",")
        
        ### Exclude SNPs strongly associated with outcome
        out_data <- subset(out_data,pval.outcome>1e-5)
        
        dat <- TwoSampleMR::harmonise_data(
          exposure_dat = exp_data,
          outcome_dat = out_data)
        
        #### Directly remove palindromic SNPs
        dat <-subset(dat,mr_keep==TRUE)
        
        
        ##### New code to calculate R2
        if(is.null(dat$beta.exposure[1])==T || is.na(dat$beta.exposure[1])==T){print("Data does not contain beta, cannot calculate F statistic")
          qaq<-qaq+2
        }else{
          dat <- get_f(dat, F_value = 10)
          
          ############## MR-PRESSO ######################    
          if (nrow(dat) > 3) {
            mr_presso_res <- mr_Presso(dat, num = 1000)
            mr_presso_main <- mr_presso_pval(mr_presso_res)
            dat <- mr_presso_snp(mr_presso_res, mr_presso_main, dat, type = "data")
            dat <-subset(dat,mr_keep==TRUE)
            
            resMRPRESSO=mr_presso_res[["Main MR results"]]
            resMRPRESSO
            global_test_p <- mr_presso_res[["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]
            se1=sqrt(((resMRPRESSO[1,3])^2)/qchisq(Pval_raw <- resMRPRESSO[1,6],1,lower.tail=F))
            se2=sqrt(((beta_cor <- resMRPRESSO[2,3])^2)/qchisq(Pval_cor <- resMRPRESSO[2,6],1,lower.tail=F))
            resMRPRESSO <- resMRPRESSO %>%
              dplyr::mutate(se = c(se1,se2))
            
            # Generate global test P and outliers
            outliers <- dat$SNP[mr_presso_res[["MR-PRESSO results"]][["Distortion Test"]][["Outliers Indices"]]]
            outliers = as.data.frame(outliers)
            outliers <- paste(outliers$outliers, collapse = ", ")
            global_test_p = as.data.frame(global_test_p)
            resMRPRESSO
            TTT <- as.data.frame(resMRPRESSO)
            TTT
            openxlsx::write.xlsx(TTT, paste0("mendelian test/",exp,"-MR-PRESSO.xlsx"))
          } 
          
          res <-choose_MR(dat=dat)
          res_hete <-purrr::map(.x=seq_along(res),
                                .f=~res[[.x]][[1]])
          res_hete <-do.call(rbind,res_hete)
          res_hete
          res1 <- generate_odds_ratios(res[[1]][[2]])
          res1
          res1$estimate <- paste0(
            format(round(res1$or, 2), nsmall = 2), " (", 
            format(round(res1$or_lci95, 2), nsmall = 2), "-",
            format(round(res1$or_uci95, 2), nsmall = 2), ")")
          res1
          print(paste0(exp,"_SNP count_",res1$nsnp[1]))
          resdata <- dat
          openxlsx::write.xlsx(dat,file = paste0("mendelian test/",exp,"-dat.xlsx"), rowNames = FALSE)
          
          openxlsx::write.xlsx(res1,paste0("mendelian test/",exp,"-res.xlsx"))
          
          
          ###### Steiger test ######
          res_steiger <- steiger_test(dat) 
          
          
          # ## Power calculation for binary outcome
          N <- dat$samplesize.outcome[1]
          alpha <- 0.05
          R2xz <- sum(dat$R2)
          K <- (dat$ncase.outcome[1] / dat$ncontrol.outcome[1])
          if (length(dat[, 1]) == 1) {
            OR <- res1 %>% filter(method == "Wald ratio") %>% pull(or)
          } else if (length(dat[, 1]) > 1) {
            OR <-
              res1 %>%filter(grepl("Inverse variance weighted", method)) %>%pull(or)
          }
          epower = NA
          power <- results_binary(N, alpha, R2xz, K, OR, epower)
          
          library(magrittr)
          # Main result 
          res3 <- res1#[1:5,]
          res3 <- res3[,-c(10:14)]
          res4 <- tidyr::pivot_wider(
            res3,names_from ="method",names_vary = "slowest",
            values_from = c("b","se","pval","estimate") )
          
          col_names <- colnames(res4)
          
          new_col_names <- gsub("\\(.*\\)", "", col_names)
          
          colnames(res4) <- new_col_names
          
          ## Steiger
          res_steiger2 <- dplyr::select(res_steiger,
                                        correct_causal_direction,steiger_pval)
          # Power
          power2 <- tidyr::pivot_wider(
            power,names_from ="Parameter",names_vary = "slowest",
            values_from = c("Value","Description") )
          power2 <- power2[,1]
          
          
          # Merge
          res_ALL <- cbind(res4, res_steiger2,power2
          )
          res_ALL$F <- mean(dat$F,na.rm=TRUE)
          res_ALL$R2 <- sum(dat$R2)
          
          
          
          if (length(dat[, 1]) > 0 && length(dat[, 1]) <= 2) {
            write.csv(res_ALL, file = paste0("mendelian test/", exp, "1.csv"), row.names = FALSE)
          } else {
            # res_hete <- TwoSampleMR::mr_heterogeneity(dat)
            res_plei <- TwoSampleMR::mr_pleiotropy_test(dat)
            res_leaveone <- mr_leaveoneout(dat)
            
            p1 <- mr_scatter_plot(res[[1]][[2]], dat)
            pdf(paste0("mendelian test/", exp, "_scatter.pdf"))
            print(p1[[1]])
            dev.off()
            
            res_single <- mr_singlesnp(dat,all_method)
            p2 <- mr_forest_plot(res_single)
            pdf(paste0("mendelian test/", exp, "_forest.pdf"))
            print(p2[[1]])
            dev.off()
            
            p3 <- mr_funnel_plot(res_single)
            pdf(paste0("mendelian test/", exp, "_funnel.pdf"))
            print(p3[[1]])
            dev.off()
            
            res_loo <- mr_leaveoneout(dat)
            pdf(paste0("mendelian test/", exp, "_leave_one_out.pdf"))
            print(mr_leaveoneout_plot(res_loo))
            dev.off()
            
            
            
            res_hete2 <- tidyr::pivot_wider(
              res_hete, names_from = "method", names_vary = "slowest",
              values_from = c("Q", "Q_df", "Q_pval")
            ) %>% 
              dplyr::select(-id.exposure, -id.outcome, -outcome, -exposure)
            res_hete2 <- res_hete2[, 4:6]
            
            res_plei2 <- dplyr::select(res_plei, egger_intercept, se, pval)
            
            res_ALL <- cbind(res_ALL, res_hete2, res_plei2)
            
            if (nrow(dat) > 3) {
              
              res_ALL$outliers <- outliers
              res_ALL <- cbind(res_ALL, global_test_p)
            }
            write.csv(res_ALL, file = paste0("mendelian test/", exp, ".csv"), row.names = FALSE)
          }
        }
      }
    }
  }
}



# Change the exported path to your own working directory  
csv_files <- list.files("E:/medelian randomization/mendelian test", pattern = "csv.csv", full.names = TRUE)      
combined_df <- read.csv(csv_files[1], stringsAsFactors = FALSE)      

# Convert the 'exposure' column to character    
if (!is.null(combined_df$exposure)) {      
  combined_df$exposure <- as.character(combined_df$exposure)      
}   

# Convert the 'global_test_p' column to character  
if (!is.null(combined_df$global_test_p)) {      
  combined_df$global_test_p <- as.character(combined_df$global_test_p)      
}  

for (i in 2:length(csv_files)) {      
  temp_df <- read.csv(csv_files[i], stringsAsFactors = FALSE)      
  
  # Convert the 'exposure' column to character    
  if (!is.null(temp_df$exposure)) {      
    temp_df$exposure <- as.character(temp_df$exposure)      
  }    
  
  # Convert the 'global_test_p' column to character  
  if (!is.null(temp_df$global_test_p)) {      
    temp_df$global_test_p <- as.character(temp_df$global_test_p)      
    temp_df$global_test_p[temp_df$global_test_p == "<"] <- ""      
  }      
  
  combined_df <- bind_rows(combined_df, temp_df)      
}


combined_df <- combined_df %>%
  mutate(fdr_Inverse.variance.weighted. = p.adjust(pval_Inverse.variance.weighted., method = "fdr")) %>%
  relocate(fdr_Inverse.variance.weighted., .after = pval_Inverse.variance.weighted.)
write.csv(combined_df,"mendelian testres.csv")
