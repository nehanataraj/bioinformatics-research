#library import
library(MASS)
library(tidyverse)
library(dplyr)
library(readxl)
library(qvalue)


####################################################
# Importing dataset with chromosome and position markers
importdata <- function(trait) {
  
  data = read_excel(paste("GWAS_Data_Reweighted/", trait, "gwas_reweighted.xlsx", sep = ""))
  
  # Merge positions and chromosome into dataset
  # variantinfofile.csv contains UK Biobank data for all chromosomes
  variantinfo = read_delim("variantinfo/variantinfofile.csv")
  
  # Switched to left_join from merge so that all values remain in the initial file
  colnames(data)
  colnames(variantinfo)
  merged_data <- left_join(x = data, y = variantinfo)
  
  write.csv(merged_data, paste("modified_datasets/", trait, "dataset.csv", sep = ""), row.names = FALSE)
  
}

importdata("BMI1")
importdata("Bipolar_Disorder")
importdata("Height1")
importdata("Ever_Smoked")
importdata("Years_of_Education1")
importdata("Crohns_Disease")
importdata("Coronary_Artery_Disease")
importdata("Ulcerative_Colitis")
importdata("Type_2_Diabetes")
importdata("Triglycerides")
importdata("Rheumatoid_Arthritis")
importdata("LDL")
importdata("HDL")
importdata("FastingGlucose_Manning")
importdata("Anorexia")

################## project part 2 ##################
####################################################
library(CAMT)
library(OPWeight)
#source("../../../Neha to read/Block Weighting Code/CalculateBlockWeight.R") #for paul computer 
source("R_files/blockweighting.R")  #for neha's computer

#changes: 
#1) Change multiple testing adjustment.  OPW and CAMT make their own adjustments. 
# while block weight does not. Input alpha should be 0.05. If you make it 5e-08 then 
#opw and CAMT will adjust of top of that and be too stringent. 

#2) The block weight p-values should be compared to alpha/m

#3) Added in results for FINDOR. This is just the p_weighted column in the data

#4) Add code to save the objects created by OPW and CAMT. These take a long time
# to run and I want to save the created objects (which include the weighted p-values)
# so that don't have to rerun if we want to do something different with the output. You must
# create a directory called output_save. Then it will create a directory inside of that 
#with the same name as the trait and write the objects to there. 

#for testing without using function:
trait="Schizophrenia"
opw_EffectType="binary"
alpha=0.05


adjpvalues <- function(trait, alpha, opw_EffectType) {

  #this line for reading on Paul's computer: 
  #data = read_delim(paste("../../../Neha to read/Data/GWAS_Data_Reweighted/", trait, "gwas_reweighted", sep = ""))
  #This line for reading on Neha's computer:
  data = read_delim(paste("modified_datasets/", trait, "dataset.csv", sep = ""))
  
  data <- data[data$PTV > 0, ]   #remove rows where PTV is negative
  m=nrow(data)
  
   
  #########################################
  (num_rej_noweight=sum(data$P<alpha/m))  #number rejected with unweighted bonferroni
  (num_rej_noweight_BH=sum(p.adjust(data$P, method="BH")<alpha)) #number rejected with unweighted FDR
  qvalsNW=qvalue(data$P)
  (num_rej_noweight_qval=sum(qvalsNW$qvalues<0.05))
  
  
  #################################
  #Block Weighting
  tol <- 0.001
  nblock <- 10
#  m <- length(data)  CORRECTION
  
  set=floor(seq(from=m/nblock,to=m,length.out=nblock))
  BW <- BlockWeight(set, alpha, tol)
  
  pval_weighted_BW=data$P/BW
  num_rej_BW_BON=sum(pval_weighted_BW<alpha/m)
  
  num_rej_BW_BH=sum(p.adjust(pval_weighted_BW, method="BH")<alpha)
  
  qval_BW=qvalue(ifelse(pval_weighted_BW<1,pval_weighted_BW,1),pi0=qvalsNW$pi0)
  num_rej_BW_qval=sum(qval_BW$qvalues<alpha)
  
  print("Block Weighting:")
  print(sum(data$P/BW<(alpha/m)))     #CHANGED to make alpha/m
  
  
  
  #################################
  #CAMT
                                                                                ##NOTE CAMT DEFAULT alpha=0.1
  camt.obj.fwer <- camt.fwer(pvals = data$P, pi0.var = data$PTV,alpha=alpha)  #changed to add argument alpha
 
  num_rej_camt_bon=length(camt.obj.fwer$rejection)
  
  print("CAMT fwer:")
  print(num_rej_camt_bon)
  
  
  camt.obj.fdr <- camt.fdr(pvals = data$P, pi0.var = data$PTV)
  
  num_rej_camt_FDR=sum(camt.obj.fdr$fdr<alpha)
  print("CAMT fdr:")
  print(num_rej_camt_FDR)
  
  
  
  ##############################
  #Optimal Covariate Weighting

   
  opw_results <- opw(pvalue = data$P, filter = data$PTV, tail=2L, alpha = alpha, effectType = opw_EffectType, method = "BON")
 
  num_rej_opw_bon=nrow(opw_results$rejections_list[,1])
  
  print("Optimal Covariate Weighting:")
  print(num_rej_opw_bon)
  
  #below I calculate weighed p-values from opw
  #the weights outputted from OPw are sorted by covariate rank. This sorts the in the same way
  # so that can calculate correct weighted p-values. orginal order column added for sorting back to original order
  
  data_sort_filter=mutate(data,"original order"=1:nrow(data))
  data_sort_filter=data_sort_filter[order(data_sort_filter$PTV,decreasing = T),]
  
  pval_weighted_opw=data_sort_filter$P/opw_results$weight
  
  
  #smoothed opw weights. opw weights flucuate presumably due to numerical error
  #this fits a smooth curve to the weights
  ssw=smooth.spline(x=1:length(opw_results$weight),y=opw_results$weight)
  pval_weighted_opw_smooth=data_sort_filter$P/ssw$y
  num_rej_opw_bon_smooth=sum(pval_weighted_opw_smooth<alpha/m)
  
  
  num_rej_opw_bon_p=sum(pval_weighted_opw<alpha/m)
  if(num_rej_opw_bon!=num_rej_opw_bon_p) print("error: two versions of opw bonferroni calc don't match!")
  
  #FDR via BH method
  padj_opw_BH=p.adjust(pval_weighted_opw,method = "BH")
  num_rej_opw_BH=sum(padj_opw_BH<0.05)  
  
  padj_opw_BH_smooth=p.adjust(pval_weighted_opw_smooth,method = "BH")
  num_rej_opw_BH_smooth=sum(padj_opw_BH_smooth<0.05)  
  
  
  #FDR via q-values
  #qvalues:
  #calc using q-values rather than BH method
  #uses pi0 calculated from unweighted p-values. 
  qvals_weighted=qvalue(ifelse(pval_weighted_opw<1,pval_weighted_opw,1),pi0=qvalsNW$pi0)
  num_rej_opw_qval=sum(qvals_weighted$qvalues<0.05)
  #is using the original pi0 valid? I think? 
  
  qvals_weighted_smooth=qvalue(ifelse(pval_weighted_opw_smooth<1,pval_weighted_opw_smooth,1),pi0=qvalsNW$pi0)
  num_rej_opw_qval_smooth=sum(qvals_weighted_smooth$qvalues<0.05)
  
  
  #This code added 2/15: Markers were previosly not sorted back into the correct order. 
  #this data frame has all markers sorted by covariate, with opw weighted p-values and qvalues)
  data_sort_filter=mutate(data_sort_filter,"opw weighted p-values"=pval_weighted_opw,"opw smoothed p-values"=pval_weighted_opw_smooth,"opw smoothed q-values"=qvals_weighted_smooth$qvalues)
  
  #this data frame is sorted back into the original order: 
  data_opw=data_sort_filter[,]
  data_opw=data_sort_filter[order(data_sort_filter$`original order`),]
  
####################################
#  FINDOR
  
  num_rej_FINDOR_bon=sum(data$P_weighted<alpha/m)
  pFIN_adj=p.adjust(data$P_weighted,method = "BH")
  num_rej_FINDOR_BH=sum(pFIN_adj<alpha)
  
  num_rej_FINDOR_qval=sum(qvalue(data$P_weighted)$qvalues<alpha)
####################################
###Save output objects to folder
  
output_dir=paste("output_save/",trait,sep="") #name of output directory

#check if directory already exists. If not, create it
if (!dir.exists(output_dir)) {dir.create(output_dir)}

save(camt.obj.fdr,camt.obj.fwer,opw_results,file=paste(trait,"_out.RData"))    
  

###################################
###Save all p-values to file
allpvals=data.frame(data$SNP, data$P,data$P_weighted,pval_weighted_BW,data_opw$`opw weighted p-values`,data_opw$`opw smoothed p-values`)
names(allpvals)=c("SNP","p-value","FINDOR","Block Weight","opw","opw smooth")
write.csv(allpvals,file=paste(output_dir,"/",trait,"_pvals.csv",sep = ""),row.names = F)

###################################
#make table with numbers of rejections:

errortype=c("FWER","BH","q-value")
methods=c("no weight","block weighting","CAMT", "OPW", "OPW smoothed","FINDOR")
rej_bon=c(num_rej_noweight,num_rej_BW_BON,num_rej_camt_bon,num_rej_opw_bon,num_rej_opw_bon_smooth, num_rej_FINDOR_bon)
rej_BH=c(num_rej_noweight_BH,num_rej_BW_BH,NA,num_rej_opw_BH,num_rej_opw_BH_smooth,num_rej_FINDOR_BH)
rej_qval=c(num_rej_noweight_qval,num_rej_BW_qval,num_rej_camt_FDR,num_rej_opw_qval,num_rej_opw_qval_smooth,num_rej_FINDOR_qval)

numrej=data.frame(methods,rej_bon,rej_BH,rej_qval) 
names(numrej)=c("Method","FWER","BH","q-value")
write.csv(numrej,file=paste(output_dir,"/",trait,"_number_rejected.csv",sep = ""),row.names = F)



########################################################################
########################################################################
# New code, turns significant values into a dataframe and stores rather than adding them

# Function to select relevant columns
select_columns <- function(df, weighted_p = NULL) {
  base_cols <- c("SNP", "chromosome", "Position", "P")
  if (!is.null(weighted_p)) {
    df$P_weighted <- weighted_p
    return(df[, c(base_cols, "P_weighted")])
  }
  return(df[, base_cols])
}

# No weighting
noweight_df <- select_columns(data[data$P < alpha/m, ])
noweight_BH_df <- select_columns(data[p.adjust(data$P, method="BH") < alpha, ])
noweight_qval_df <- select_columns(data[qvalue(data$P)$qvalues < 0.05, ])

# Block weighting
blockweight_df <- select_columns(data[pval_weighted_BW < alpha/m, ], pval_weighted_BW[pval_weighted_BW < alpha/m])
blockweight_BH_df <- select_columns(data[p.adjust(pval_weighted_BW, method="BH") < alpha, ], 
                                    pval_weighted_BW[p.adjust(pval_weighted_BW, method="BH") < alpha])
blockweight_qval_df <- select_columns(data[qvalue(ifelse(pval_weighted_BW < 1, pval_weighted_BW, 1))$qvalues < alpha, ],
                                      pval_weighted_BW[qvalue(ifelse(pval_weighted_BW < 1, pval_weighted_BW, 1))$qvalues < alpha])

# CAMT
camt_fwer_df <- select_columns(data[camt.obj.fwer$rejection, ])
camt_fdr_df <- select_columns(data[camt.obj.fdr$fdr < alpha, ])

# For OPW, sorted data
sorted_data <- data[order(data$PTV, decreasing = TRUE), ]
opw_df <- select_columns(sorted_data[pval_weighted_opw < alpha/m, ], pval_weighted_opw[pval_weighted_opw < alpha/m])
opw_BH_df <- select_columns(sorted_data[p.adjust(pval_weighted_opw, method="BH") < alpha, ],
                            pval_weighted_opw[p.adjust(pval_weighted_opw, method="BH") < alpha])
opw_qval_df <- select_columns(sorted_data[qvalue(ifelse(pval_weighted_opw < 1, pval_weighted_opw, 1))$qvalues < alpha, ],
                              pval_weighted_opw[qvalue(ifelse(pval_weighted_opw < 1, pval_weighted_opw, 1))$qvalues < alpha])

opw_smooth_df <- select_columns(sorted_data[pval_weighted_opw_smooth < alpha/m, ], 
                                pval_weighted_opw_smooth[pval_weighted_opw_smooth < alpha/m])
opw_smooth_BH_df <- select_columns(sorted_data[p.adjust(pval_weighted_opw_smooth, method="BH") < alpha, ],
                                   pval_weighted_opw_smooth[p.adjust(pval_weighted_opw_smooth, method="BH") < alpha])
opw_smooth_qval_df <- select_columns(sorted_data[qvalue(ifelse(pval_weighted_opw_smooth < 1, pval_weighted_opw_smooth, 1))$qvalues < alpha, ],
                                     pval_weighted_opw_smooth[qvalue(ifelse(pval_weighted_opw_smooth < 1, pval_weighted_opw_smooth, 1))$qvalues < alpha])

# FINDOR (using P_weighted from original data)
findor_df <- select_columns(data[data$P_weighted < alpha/m, ], data$P_weighted[data$P_weighted < alpha/m])
findor_BH_df <- select_columns(data[p.adjust(data$P_weighted, method="BH") < alpha, ],
                               data$P_weighted[p.adjust(data$P_weighted, method="BH") < alpha])
findor_qval_df <- select_columns(data[qvalue(data$P_weighted)$qvalues < alpha, ],
                                 data$P_weighted[qvalue(data$P_weighted)$qvalues < alpha])

# Save to CSV files
write.csv(noweight_df, file=paste(output_dir,"/",trait,"_noweight_FWER.csv",sep=""), row.names=FALSE)
write.csv(noweight_BH_df, file=paste(output_dir,"/",trait,"_noweight_BH.csv",sep=""), row.names=FALSE)
write.csv(noweight_qval_df, file=paste(output_dir,"/",trait,"_noweight_qval.csv",sep=""), row.names=FALSE)

write.csv(blockweight_df, file=paste(output_dir,"/",trait,"_blockweight_FWER.csv",sep=""), row.names=FALSE)
write.csv(blockweight_BH_df, file=paste(output_dir,"/",trait,"_blockweight_BH.csv",sep=""), row.names=FALSE)
write.csv(blockweight_qval_df, file=paste(output_dir,"/",trait,"_blockweight_qval.csv",sep=""), row.names=FALSE)

write.csv(camt_fwer_df, file=paste(output_dir,"/",trait,"_camt_FWER.csv",sep=""), row.names=FALSE)
write.csv(camt_fdr_df, file=paste(output_dir,"/",trait,"_camt_FDR.csv",sep=""), row.names=FALSE)

write.csv(opw_df, file=paste(output_dir,"/",trait,"_opw_FWER.csv",sep=""), row.names=FALSE)
write.csv(opw_BH_df, file=paste(output_dir,"/",trait,"_opw_BH.csv",sep=""), row.names=FALSE)
write.csv(opw_qval_df, file=paste(output_dir,"/",trait,"_opw_qval.csv",sep=""), row.names=FALSE)

write.csv(opw_smooth_df, file=paste(output_dir,"/",trait,"_opw_smooth_FWER.csv",sep=""), row.names=FALSE)
write.csv(opw_smooth_BH_df, file=paste(output_dir,"/",trait,"_opw_smooth_BH.csv",sep=""), row.names=FALSE)
write.csv(opw_smooth_qval_df, file=paste(output_dir,"/",trait,"_opw_smooth_qval.csv",sep=""), row.names=FALSE)

write.csv(findor_df, file=paste(output_dir,"/",trait,"_findor_FWER.csv",sep=""), row.names=FALSE)
write.csv(findor_BH_df, file=paste(output_dir,"/",trait,"_findor_BH.csv",sep=""), row.names=FALSE)
write.csv(findor_qval_df, file=paste(output_dir,"/",trait,"_findor_qval.csv",sep=""), row.names=FALSE)

# Added code 2/21, moves to subdirectories
# Define subdirectories based on methods
methods <- c("blockweight", "camt", "findor", "noweight", "opw", "opw_smooth")

# Create subdirectories if they do not exist
for (method in methods) {
  method_dir <- file.path(output_dir, method)
  if (!dir.exists(method_dir)) {
    dir.create(method_dir)
  }
}

# Move files to appropriate subdirectories
move_files <- function(pattern, from_dir, to_dir) {
  files <- list.files(from_dir, pattern = pattern, full.names = TRUE)
  for (file in files) {
    new_path <- file.path(to_dir, basename(file))
    file.rename(file, new_path)
  }
}

# Move files for each method
move_files("noweight_.*\\.csv$", output_dir, file.path(output_dir, "noweight"))
move_files("blockweight_.*\\.csv$", output_dir, file.path(output_dir, "blockweight"))
move_files("camt_.*\\.csv$", output_dir, file.path(output_dir, "camt"))
move_files("findor_.*\\.csv$", output_dir, file.path(output_dir, "findor"))
move_files("opw_smooth_.*\\.csv$", output_dir, file.path(output_dir, "opw_smooth"))
move_files("opw_[^s].*\\.csv$", output_dir, file.path(output_dir, "opw"))  # Move non-smooth OPW files

#############################################

#markerdat=read.table("../../Data/UK Biobank/snps.genotyped.chr1.csv",header=T)

#markerdat=read_delim("../../../../1000Genomes/Data/SNPChrPosAllele_b129/SNPChrPosAllele_b129.txt",col_names = c("chr","chrPosFrom", "chrPosTo", "rs","ChrAllele","variantAllele","snpAlleleChrOrien", "snp2chrOrien", "snpClassAbbrev","snpClassCode","mapLocType","mapLocCnt","mapWeight","contigLabel","uplacedContig"),skip=4,delim="|")



#sig_noweight=allpvals[allpvals$`p-value`<alpha/m,]
#sig_noweight[3,1]%in%markerdat[,4]

#markermatch=match(sig_noweight[,1],as.vector(markerdat[,4]))


#sig_opw=allpvals[allpvals$`opw smooth`<alpha/m,]



return(list(compare_num_rej=numrej,all_pvalues=allpvals))
} #end adjpvalues function


# Running analysis
schizophrenia = adjpvalues("Schizophrenia", 0.05, "binary")
save(schizophrenia, file = "schizophrenia_results")
load("schizophrenia_results")

bmi1 = adjpvalues("BMI1", 0.05, "continous")
save(bmi1, file = "bmi1_results")
load("bmi1_results")

bipolar_disorder = adjpvalues("Bipolar_Disorder", 0.05, "binary")
save(bipolar_disorder, file = "bipolar_disorder_results")
load("bipolar_disorder_results")

height1 = adjpvalues("Height1", 0.05, "continuous")
save(height1, file = "height1_results")
load("height1_results")

ever_smoked = adjpvalues("Ever_Smoked", 0.05, "binary")
save(ever_smoked, file = "ever_smoked_results")
load("ever_smoked_results")
#many values returned as zero

years_of_education1 = adjpvalues("Years_of_Education1", 0.05, "binary")
save(years_of_education1, file = "years_of_education1_results")
load("years_of_education1_results")

crohns_disease = adjpvalues("Crohns_Disease", 0.05, "binary")
save(crohns_disease, file = "crohns_disease_results")
load("crohns_disease_results")

coronary_artery_disease = adjpvalues("Coronary_Artery_Disease", 0.05, "binary")
save(coronary_artery_disease, file = "coronary_artery_disease_results")
load("coronary_artery_disease_results")

ulcerative_colitis = adjpvalues("Ulcerative_Colitis", 0.05, "binary")
save(ulcerative_colitis, file = "ulcerative_colitis_results")
load("ulcerative_colitis_results")

type_2_diabetes = adjpvalues("Type_2_Diabetes", 0.05, "binary")
save(type_2_diabetes, file = "type_2_diabetes_results")
load("type_2_diabetes_results")

triglycerides = adjpvalues("Triglycerides", 0.05, "continuous")
save(triglycerides, file = "triglycerides_results")
load("triglycerides_results")

rheumatoid_arthritis = adjpvalues("Rheumatoid_Arthritis", 0.05, "binary")
save(rheumatoid_arthritis, file = "rheumatoid_arthritis_results")
load("rheumatoid_arthritis_results")

LDL = adjpvalues("LDL", 0.05, "binary")
save(LDL, file = "LDL_results")
load("LDL_results ")

fastingglucose_manning = adjpvalues("FastingGlucose_Manning", 0.05, "continuous")
save(fastingglucose_manning, file = "fastingglucose_manning_results")
load("fastingglucose_manning_results")

anorexia = adjpvalues("Anorexia", 0.05, "binary")
save(anorexia, file = "anorexia_results")
load("anorexia_results")
#many tests returned zero


