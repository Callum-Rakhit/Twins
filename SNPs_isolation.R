install.packages("bigmemory")
library(bigmemory)

SNP_list_1 <- read.delim("~/Documents/twins_ML_project/plink/1-4-T23-gwa610K.output.raw", sep = " ", nrows = 1)
SNP_list_2 <- read.delim("~/Documents/twins_ML_project/plink/11-22-T23-gwa610K.output.raw", sep = " ", nrows = 1)
SNP_list_3 <- read.delim("~/Documents/twins_ML_project/plink/5-10-T23-gwa610K.output.raw", sep = " ", nrows = 1)

bigmemory::as.big.matrix(x = )  # Can only be used on matrices, which can only have one type of data

