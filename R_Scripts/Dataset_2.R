##### Install required libraries
GetPackages <- function(required.packages) {
  packages.not.installed <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if(length(packages.not.installed)){install.packages(packages.not.installed, dependencies = T)}
  suppressMessages(lapply(required.packages, require, character.only = T))
}

GetPackages(c("tidyverse", "randomForest", "dplyr", "igraph", "caret", "reshape", "rsample",
              "reshape2", "devtools", "PerformanceAnalytics", "ggplot2", "lubridate",
              "bootstrap", "corrplot", "ggraph", "doParallel", "ranger", "data.table", "h2o",
              "sparsio"))

devtools::install_github("krlmlr/ulimit")

##### Load in the genotype matrices #####

SNP_list <- fread(input = "~/Stephan_SNP_HbSS_Dataset/Dataset_2.raw")  # Load the plink data
SNP_list <- SNP_list[, grep("HET", colnames(SNP_list)):=NULL]  # Remove the HET columns (cut dataframe size in half)
SNP_list$IID <- as.numeric(gsub("ID-([0-9]+)", "\\1", SNP_list$IID)) # Convert ID-001 to 001 for matching with FACS
SNP_Fcell <- fread(input = "~/Twins/Data/FID_ID_FCellFACS.csv")  # Load fcell level information
SNP_merged <- merge(x = SNP_list, y = SNP_Fcell, by = "IID")  # Add the fcell information to the SNP data

SNP_merged <- SNP_merged[, grep("AT", colnames(SNP_merged)):=NULL]
SNP_merged <- SNP_merged[, grep("SEX", colnames(SNP_merged)):=NULL]
SNP_merged <- SNP_merged[, grep("PHENOTYPE", colnames(SNP_merged)):=NULL]
SNP_merged$FID.x <- NULL  # Remove FID.x

saveRDS(SNP_merged, file = "~/Twins/Data/Dataset_2.rds")
SNP_merged <- readRDS(file = "~/Twins/Data/Dataset_2.rds")



# SNP_merged[,1:5][1,] # Checking columns names
# SNP_merged[,dim(SNP_merged_1)[2]-5:dim(SNP_merged_1)][1,] # Checking columns names

# IDs <- SNP_merged[,1:2]
# FACS_info <- SNP_merged[,dim(SNP_merged)[2]:dim(SNP_merged)[2]]
# SNP_merged_1 <- SNP_merged[,3:round(dim(SNP_merged)[2]/3)]
# SNP_merged_2 <- SNP_merged[,(round(dim(SNP_merged)[2]/3)+1):(2*round(dim(SNP_merged)[2]/3))]
# SNP_merged_3 <- SNP_merged[,((2*round(dim(SNP_merged)[2]/3))+1):((dim(SNP_merged)[2])-1)]
# 
# # Part 1
# SNP_merged_1 <- cbind(IDs, SNP_merged_1, FACS_info)
# 
# colnames(SNP_merged)


