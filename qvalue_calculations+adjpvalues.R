#library import
library(MASS)
library(tidyverse)
library(dplyr)
library(readxl)
library(qvalue)

################## project part 1 ##################
####################################################
#Determine how q-value varies by PTV

####################################################
#q-value installation
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#
#BiocManager::install("qvalue")


####################################################
#combine variantinfo into a single file
variantinfo <- NULL
columns_variantinfo = c("SNP", "Position")
for (i in 1:22) {
  
  variantinfo_chromosome = read_delim(paste("variantinfo/snps.genotyped.chr", i, ".csv", sep = ""), col_select = columns_variantinfo)
  variantinfo_chromosome <- variantinfo_chromosome %>% mutate(chromosome = i)
  
  variantinfo <- rbind(variantinfo, variantinfo_chromosome) 
  write.csv(variantinfo, "variantinfo/variantinfofile.csv", row.names = FALSE)
}

####################################################
#trait data import
importdata <- function(trait) {
  data = read_excel(paste("GWAS_Data_Reweighted/", "Schizophrenia", "gwas_reweighted.xlsx", sep = ""))
  head(data)
  
  #merge positions and chromosome into dataset
  variantinfo = read_delim("variantinfo/variantinfofile.csv")
  
  #switched to left_join from merge so that all values remain in the initial file
  merged_data <- left_join(x = data, y = variantinfo)
  #head(merged_data)
  
  #Calculate QValues and add q value column to merged data
  p_values <- merged_data$`P`
  q_values_results <- qvalue(p_values)
  q_values <- q_values_results$qvalues
  merged_data$q_values <- q_values
  
  
  #sort dataset by qvalue
  sorted_data <- merged_data[order(merged_data$q_values), ]
  write.csv(sorted_data, paste("modified_datasets/", trait, "dataset.csv", sep = ""), row.names = FALSE)
  
}

####################################################
#prints PTV values for groups of q values
covariate <- function(trait, num) {
  data <- read_delim(paste("modified_datasets/", trait, "dataset.csv", sep = ""))
  
  # Calculate the number of rows per group
  rows_per_group <- ceiling(nrow(data) / num)
  
  result <- data %>%
    mutate(Group = ceiling(row_number() / rows_per_group)) %>%
    group_by(Group) %>%
    summarize(Average_PTV = mean(PTV, na.rm = TRUE))
  
  return(result)
}


################## project part 2 ##################
####################################################
library("CAMT")
library("OPWeight")
source("blockweighting.R")

adjpvalues <- function(trait, alpha, opw_EffectType) {
  data = read_delim(paste("modified_datasets/", trait, "dataset.csv", sep = ""))
  
  #################################
  #Block Weighting
  tol <- 0.001
  nblock <- 10
  m <- length(data)
  
  set=floor(seq(from=m/nblock,to=m,length.out=nblock))
  BW <- BlockWeight(set, alpha, tol)
  
  print("Block Weighting:")
  print(sum(data$P/BW<alpha))   
  
  #################################
  #CAMT
  
  camt.obj.fwer <- camt.fwer(pvals = data$P, pi0.var = data$PTV)
  print("CAMT fwer:")
  print(length(camt.obj.fwer$rejection))
  
  camt.obj.fdr <- camt.fdr(pvals = data$P, pi0.var = data$PTV)
  print("CAMT fdr:")
  print(sum(camt.obj.fdr$fdr<alpha))
  
  ##############################
  #Optimal Covariate Weighting
  data <- data[data$PTV > 0, ]
  opw_results <- opw(pvalue = data$P, filter = data$PTV, alpha = alpha, effectType = opw_EffectType, method = "BON")
  print("Optimal Covariate Weighting:")
  print(sum(opw_results$rejections_list[,1]<alpha))
  
}

schizophrenia = adjpvalues("Schizophrenia", 5e-8, "binary")
save(schizophrenia, file = "schizophrenia_results")
load("schizophrenia_results")

