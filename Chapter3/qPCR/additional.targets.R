######### This file calculates relative expression values from qPCR runs.
##### 1. data needs to be organsised
##### 2. individual reactions that failed need to be kicked out
##### 3. calculation of Cq values (2nd derivative of amplification curve)

### 1. organising sheets

# import of useful packages
library(qpcR)
setwd("//file/Usersm$/mmi99/Home/Desktop/qPCR")
## works from my MacBook, but not from work computer
setwd("~/Google Drive/qPCR")

### normalising SYBR green data against ROX data and then transpose data, so it can be used for qPCR script (qPCR_May2020_...R)

pm46_96 <- read.csv("additional_targets/processed_data/csv_files/PM46_96hpi.csv")
pm47_96 <- read.csv("additional_targets/processed_data/csv_files/PM47_96hpi.csv")
pm48_96 <- read.csv("additional_targets/processed_data/csv_files/PM48_96hpi.csv")
pm51_96 <- read.csv("additional_targets/processed_data/csv_files/PM51_96hpi.csv")
pm52_96 <- read.csv("additional_targets/processed_data/csv_files/PM52_96hpi.csv")
pm53_96 <- read.csv("additional_targets/processed_data/csv_files/PM53_96hpi.csv")
pm54_96 <- read.csv("additional_targets/processed_data/csv_files/PM54_96hpi.csv")
pm55_96 <- read.csv("additional_targets/processed_data/csv_files/PM55_96hpi.csv")
pm60_96 <- read.csv("additional_targets/processed_data/csv_files/PM60_96hpi.csv")
pm62_96 <- read.csv("additional_targets/processed_data/csv_files/PM62_96hpi.csv")

##### nt function to normalise and transpose
nt <- function(x){
  sybr <- x[2:73, 3:42]
  rox <- x[81:152, 3:42]
  norm <- sybr/rox
  t <- t(norm)
  #write.table(x, file = "", sep="\t",row.names=FALSE, col.names = FALSE)
  return(t)
}

setwd("~/Google Drive/qPCR/additional_targets/processed_data/txt_files")

write.table(nt(pm46_96), file = "pm46_96.txt", sep="\t",row.names=FALSE, col.names = FALSE)
write.table(nt(pm47_96), file = "pm47_96.txt", sep="\t",row.names=FALSE, col.names = FALSE)
write.table(nt(pm48_96), file = "pm48_96.txt", sep="\t",row.names=FALSE, col.names = FALSE)
write.table(nt(pm51_96), file = "pm51_96.txt", sep="\t",row.names=FALSE, col.names = FALSE)
write.table(nt(pm52_96), file = "pm52_96.txt", sep="\t",row.names=FALSE, col.names = FALSE)
write.table(nt(pm53_96), file = "pm53_96.txt", sep="\t",row.names=FALSE, col.names = FALSE)
write.table(nt(pm54_96), file = "pm54_96.txt", sep="\t",row.names=FALSE, col.names = FALSE)
write.table(nt(pm55_96), file = "pm55_96.txt", sep="\t",row.names=FALSE, col.names = FALSE)
write.table(nt(pm60_96), file = "pm60_96.txt", sep="\t",row.names=FALSE, col.names = FALSE)
write.table(nt(pm62_96), file = "pm62_96.txt", sep="\t",row.names=FALSE, col.names = FALSE)

### 2. getting rid of individual reactions that failed
### see spreadsheet 'additional_targets_qPCR_outlier_check'

pm46_96 <- nt(pm46_96)
pm47_96 <- nt(pm47_96)
pm48_96 <- nt(pm48_96)
pm51_96 <- nt(pm51_96)
pm52_96 <- nt(pm52_96)
pm53_96 <- nt(pm53_96)
pm54_96 <- nt(pm54_96)
pm55_96 <- nt(pm55_96)
pm60_96 <- nt(pm60_96)
pm62_96 <- nt(pm62_96)


### sometimes the algorithm fails to detect reactions that don't amplify or amplify late as outliers
### --> 
### change reactions with final fluorescence value below 5 to a vector of 1s

amp5 <- function(x){
  for (i in 1:72){
    if (x[40,i] < 5) {
      for (j in 1:40){
        x[j,i] = 1
      }
    } 
  }
  return(x)
}

pm46_96 <- amp5(pm46_96)
pm47_96 <- amp5(pm47_96)
pm48_96 <- amp5(pm48_96)
pm51_96 <- amp5(pm51_96)
pm52_96 <- amp5(pm52_96)
pm53_96 <- amp5(pm53_96)
pm54_96 <- amp5(pm54_96)
pm55_96 <- amp5(pm55_96)
pm60_96 <- amp5(pm60_96)
pm62_96 <- amp5(pm62_96)

#################################################################################################
### getting rid of individual reactions that were either gDNA contaminated or had wrong melt peak 
#################################################################################################
###function to change reaction x to a vector of 1s
cr40 <- function(x, y){
  for (i in 1:40){
    x[i,y] = 1
  }
  return(x)
}



setwd("~/Google Drive/qPCR/additional_targets/processed_data/txt_files/ind_react_change")

pm46_96 <- cr40(pm46_96, 4)
pm46_96 <- cr40(pm46_96, 8)
pm46_96 <- cr40(pm46_96, 12)
pm46_96 <- cr40(pm46_96, 16)
pm46_96 <- cr40(pm46_96, 20)
pm46_96 <- cr40(pm46_96, 24)
pm46_96 <- cr40(pm46_96, 28)
pm46_96 <- cr40(pm46_96, 32)
pm46_96 <- cr40(pm46_96, 36)
pm46_96 <- cr40(pm46_96, 40)
pm46_96 <- cr40(pm46_96, 44)
pm46_96 <- cr40(pm46_96, 48)
pm46_96 <- cr40(pm46_96, 52)
pm46_96 <- cr40(pm46_96, 56)
pm46_96 <- cr40(pm46_96, 60)
pm46_96 <- cr40(pm46_96, 64)
pm46_96 <- cr40(pm46_96, 68)
pm46_96 <- cr40(pm46_96, 69)
pm46_96 <- cr40(pm46_96, 70)
pm46_96 <- cr40(pm46_96, 71)
pm46_96 <- cr40(pm46_96, 72)
write.table(pm46_96, file = "pm46_96.txt", sep="\t",row.names=FALSE, col.names = FALSE)

write.table(pm47_96, file = "pm47_96.txt", sep="\t",row.names=FALSE, col.names = FALSE)

pm48_96 <- cr40(pm48_96, 3)
pm48_96 <- cr40(pm48_96, 4)
pm48_96 <- cr40(pm48_96, 8)
pm48_96 <- cr40(pm48_96, 12)
pm48_96 <- cr40(pm48_96, 16)
pm48_96 <- cr40(pm48_96, 20)
pm48_96 <- cr40(pm48_96, 24)
pm48_96 <- cr40(pm48_96, 28)
pm48_96 <- cr40(pm48_96, 32)
pm48_96 <- cr40(pm48_96, 36)
pm48_96 <- cr40(pm48_96, 40)
pm48_96 <- cr40(pm48_96, 44)
pm48_96 <- cr40(pm48_96, 48)
pm48_96 <- cr40(pm48_96, 52)

pm48_96 <- cr40(pm48_96, 60)

pm48_96 <- cr40(pm48_96, 68)
pm48_96 <- cr40(pm48_96, 69)
pm48_96 <- cr40(pm48_96, 70)
pm48_96 <- cr40(pm48_96, 71)
pm48_96 <- cr40(pm48_96, 72)
write.table(pm48_96, file = "pm48_96.txt", sep="\t",row.names=FALSE, col.names = FALSE)

write.table(pm51_96, file = "pm51_96.txt", sep="\t",row.names=FALSE, col.names = FALSE)

write.table(pm52_96, file = "pm52_96.txt", sep="\t",row.names=FALSE, col.names = FALSE)

pm53_96 <- cr40(pm53_96, 5)
pm53_96 <- cr40(pm53_96, 6)
pm53_96 <- cr40(pm53_96, 7)
pm53_96 <- cr40(pm53_96, 8)
pm53_96 <- cr40(pm53_96, 13)
pm53_96 <- cr40(pm53_96, 14)
pm53_96 <- cr40(pm53_96, 15)
pm53_96 <- cr40(pm53_96, 16)
pm53_96 <- cr40(pm53_96, 24)
pm53_96 <- cr40(pm53_96, 29)
pm53_96 <- cr40(pm53_96, 30)
pm53_96 <- cr40(pm53_96, 31)
pm53_96 <- cr40(pm53_96, 32)
pm53_96 <- cr40(pm53_96, 37)
pm53_96 <- cr40(pm53_96, 38)
pm53_96 <- cr40(pm53_96, 39)
pm53_96 <- cr40(pm53_96, 40)
write.table(pm53_96, file = "pm53_96.txt", sep="\t",row.names=FALSE, col.names = FALSE)

pm54_96 <- cr40(pm54_96, 24)
pm54_96 <- cr40(pm54_96, 56)
pm54_96 <- cr40(pm54_96, 64)
write.table(pm54_96, file = "pm54_96.txt", sep="\t",row.names=FALSE, col.names = FALSE)

pm55_96 <- cr40(pm55_96, 4)
pm55_96 <- cr40(pm55_96, 8)
pm55_96 <- cr40(pm55_96, 12)
pm55_96 <- cr40(pm55_96, 16)
pm55_96 <- cr40(pm55_96, 20)
pm55_96 <- cr40(pm55_96, 24)
pm55_96 <- cr40(pm55_96, 28)
pm55_96 <- cr40(pm55_96, 32)

pm55_96 <- cr40(pm55_96, 40)
pm55_96 <- cr40(pm55_96, 44)
pm55_96 <- cr40(pm55_96, 48)
pm55_96 <- cr40(pm55_96, 52)
pm55_96 <- cr40(pm55_96, 56)
pm55_96 <- cr40(pm55_96, 60)
pm55_96 <- cr40(pm55_96, 64)
pm55_96 <- cr40(pm55_96, 68)
pm55_96 <- cr40(pm55_96, 69)
pm55_96 <- cr40(pm55_96, 70)
pm55_96 <- cr40(pm55_96, 71)
pm55_96 <- cr40(pm55_96, 72)
write.table(pm55_96, file = "pm55_96.txt", sep="\t",row.names=FALSE, col.names = FALSE)

pm60_96 <- cr40(pm60_96, 32)
pm60_96 <- cr40(pm60_96, 36)
write.table(pm60_96, file = "pm60_96.txt", sep="\t",row.names=FALSE, col.names = FALSE)

pm62_96 <- cr40(pm62_96, 5)
pm62_96 <- cr40(pm62_96, 6)
pm62_96 <- cr40(pm62_96, 7)
pm62_96 <- cr40(pm62_96, 8)
pm62_96 <- cr40(pm62_96, 12)
pm62_96 <- cr40(pm62_96, 24)
pm62_96 <- cr40(pm62_96, 25)
pm62_96 <- cr40(pm62_96, 26)
pm62_96 <- cr40(pm62_96, 27)
pm62_96 <- cr40(pm62_96, 28)
pm62_96 <- cr40(pm62_96, 32)
pm62_96 <- cr40(pm62_96, 36)
pm62_96 <- cr40(pm62_96, 40)
pm62_96 <- cr40(pm62_96, 44)
pm62_96 <- cr40(pm62_96, 48)
pm62_96 <- cr40(pm62_96, 52)
pm62_96 <- cr40(pm62_96, 60)
write.table(pm62_96, file = "pm62_96.txt", sep="\t",row.names=FALSE, col.names = FALSE)

### 3. This part calculates relative expression values from qPCR runs.
## Cq value is the 2nd derivative of the amplification curve.
## cDNA amount is calculated using the primer efficiency.
## gDNA amount is subtracted from cDNA amount (RT - noRT).
## cDNA amount is normalised to  interrun calibrator.

# import of useful packages
library(qpcR)
library(tidyverse)
library(ggplot2)
library(tidyr)
library(ComplexHeatmap)

######################################################################################
# all individual functions are listed here. They are then used in the masterfunction #
######################################################################################

### function that calculates the Cq value (in this specific case the 2nd derivative of the amplification curve)
getCpD2 <- function (x, type = c("fit", "curve"), cp = "cpD2", eff = "sigfit", 
                     ...) 
{
  type <- match.arg(type)
  options(expressions = 20000)
  if (class(x)[1] == "pcrfit") 
    x <- modlist(x)
  if (class(x)[1] != "modlist") 
    stop("'x' must be either of class 'pcrfit' or 'modlist'!")
  modNAMES <- sapply(x, function(x) x$MODEL$name)
  modNAMES <- unique(modNAMES[!is.na(modNAMES)])
  if (modNAMES %in% c("mak2", "mak3", "chag", "cm3")) 
    type <- "fit"
  if (type == "fit") {
    RN <- lapply(x, function(x) names(coef(x)))
    RN <- unique(unlist(RN))
    RN <- na.omit(RN)
    NR = length(RN)
  }
  else {
    RN <- "cpD2"
    NR = 1
  }
  RES <- matrix(nrow = NR, ncol = length(x))
  NAMES <- sapply(x, function(a) a$name)
  for (i in 1:length(x)) {
    #counter(i)
    flush.console()
    tempMOD <- x[[i]]
    if (type == "fit") {
      COEF <- try(coef(tempMOD), silent = TRUE)
      if (inherits(COEF, "try-error")) 
        COEF <- NA
      if (length(COEF) > 0) 
        RES[1:length(COEF), i] <- COEF
    }
    if (type == "curve") {
      outNAME <- switch(cp, cpD2 = "cpD2", cpD1 = "cpD1", 
                        maxE = "cpE", expR = "cpR", Cy0 = "Cy0", CQ = "cpCQ", 
                        maxRatio = "cpMR", stop())
      tempRES <- tryCatch(efficiency(tempMOD, plot = FALSE, 
                                     type = cp, ...), error = function(e) NA)
      tempCT <- tryCatch(tempRES[[outNAME]], error = function(e) NA)
      RES[1, i] <- tempCT
      
    }
  }
  colnames(RES) <- NAMES
  rownames(RES) <- RN
  cat("\n")
  return(RES)
}

### function that organises replicates
orgRep <- function(x){
  mock_RT <- c(mean(c(x[1,1], x[1,2], x[1,3]), na.rm = TRUE), 
               mean(c(x[1,9], x[1,10], x[1,11]), na.rm = TRUE),
               mean(c(x[1,17], x[1,18], x[1,19]), na.rm = TRUE))
  
  mock_noRT <- c(x[1,4], x[1,12], x[1,20])
  
  Sphingomonas_leaf34_RT <- c(mean(c(x[1,25], x[1,26], x[1,27]), na.rm = TRUE), 
                              mean(c(x[1,33], x[1,34], x[1,35]), na.rm = TRUE),
                              mean(c(x[1,41], x[1,42], x[1,43]), na.rm = TRUE))
  
  Sphingomonas_leaf34_noRT <- c(x[1,28], x[1,36], x[1,44])
  
  Microbacterium_leaf347_RT <- c(mean(c(x[1,49], x[1,50], x[1,51]), na.rm = TRUE), 
                                 mean(c(x[1,57], x[1,58], x[1,59]), na.rm = TRUE),
                                 mean(c(x[1,65], x[1,66], x[1,67]), na.rm = TRUE))
  
  Microbacterium_leaf347_noRT <- c(x[1,52], x[1,60], x[1,68])
  
  Williamsia_leaf354_RT <- c(mean(c(x[1,5], x[1,6], x[1,7]), na.rm = TRUE), 
                             mean(c(x[1,13], x[1,14], x[1,15]), na.rm = TRUE),
                             mean(c(x[1,21], x[1,22], x[1,23]), na.rm = TRUE))
  
  Williamsia_leaf354_noRT <- c(x[1,8], x[1,16], x[1,24])
  
  Pseudomonas_syringae_RT <- c(mean(c(x[1,29], x[1,30], x[1,31]), na.rm = TRUE), 
                               mean(c(x[1,37], x[1,38], x[1,39]), na.rm = TRUE),
                               mean(c(x[1,45], x[1,46], x[1,47]), na.rm = TRUE))
  
  Pseudomonas_syringae_noRT <- c(x[1,32], x[1,40], x[1,48])
  
  interrun_RT <- c(mean(c(x[1,53], x[1,54], x[1,55], x[1,61], x[1,62], x[1,63]), na.rm = TRUE), NA, NA)
  interrun_noRT <- c(mean(c(x[1,56], x[1,64])), NA, NA)
  
  RT<- data.frame(mock_RT = mock_RT, Sphingomonas_leaf34_RT = Sphingomonas_leaf34_RT, Microbacterium_leaf347_RT = Microbacterium_leaf347_RT, Williamsia_leaf354_RT = Williamsia_leaf354_RT, Pseudomonas_syringae_RT = Pseudomonas_syringae_RT,
                  mock_noRT = mock_noRT, Sphingomonas_leaf34_noRT = Sphingomonas_leaf34_noRT, Microbacterium_leaf347_noRT = Microbacterium_leaf347_noRT, Williamsia_leaf354_noRT = Williamsia_leaf354_noRT, Pseudomonas_syringae_noRT = Pseudomonas_syringae_noRT,
                  interrun_RT = interrun_RT, interrun_noRT = interrun_noRT)
  return(RT)
}

### fucntion to calculate template amount, subtract no_RT and normalise to interrun calibrator
pcrnorm <- function(df, primer_eff){
  
  ## calculating template amount
  df_t <- 1/primer_eff^df
  
  
  ## minus noRT
  df_t[is.na(df_t)] <- 0 # repalcing NA for 0s
  df_rt <- df_t[,1:5]
  df_nt <- df_t[,6:10]
  df_ct <- df_rt-df_nt
  
  
  ## norm to interrun calibrator
  df_nct <- df_ct/df_t[1,11]
  
  return(df_nct)
}

######################################################################################
# masterfunction #
######################################################################################

masterfunction <- function(x, y){
  
  ## create a modlist
  ml <- modlist(x[[1]], model = b7, norm = TRUE)
  
  ## calculate Cq using 'getCpD2'
  cpD2 <- getCpD2(ml, type = "curve")
  
  ## get rid of Cq values above 39 and below 2
  cpD2[cpD2 > 39] <- NA
  cpD2[cpD2 <2] <- NA
  
  ## organise replicates using 'orgRep' function
  org <- orgRep(cpD2) 
  
  ## calculate template amount, subtract gDNA and normalise to interrun calibrator using 'pcrnorm' function
  dfp <- pcrnorm(org, y)
  
  ## return dfp
  return(dfp)
}

######################################################################################
# small version of 'masterfunction' to retireve Cq-values #
######################################################################################

ma_Cq <- function(x, y){
  
  ## create a modlist
  ml <- modlist(x[[1]], model = b7, norm = TRUE)
  
  ## calculate Cq using 'getCpD2'
  cpD2 <- getCpD2(ml, type = "curve")
  
  ## get rid of Cq values above 39 and below 2
  cpD2[cpD2 > 39] <- NA
  cpD2[cpD2 <2] <- NA
  
  ## organise replicates using 'orgRep' function
  org <- orgRep(cpD2) 
  
  ## return dfp
  return(org)
}

######################################################################################
# loading and processing data #
######################################################################################
setwd("~/Google Drive/qPCR")
# load data into R
pm46_96 <- read.delim("additional_targets/processed_data/txt_files/ind_react_change/pm46_96.txt", header = FALSE)
pm47_96 <- read.delim("additional_targets/processed_data/txt_files/ind_react_change/pm47_96.txt", header = FALSE)
pm48_96 <- read.delim("additional_targets/processed_data/txt_files/ind_react_change/pm48_96.txt", header = FALSE)
pm51_96 <- read.delim("additional_targets/processed_data/txt_files/ind_react_change/pm51_96.txt", header = FALSE)
pm52_96 <- read.delim("additional_targets/processed_data/txt_files/ind_react_change/pm52_96.txt", header = FALSE)
pm53_96 <- read.delim("additional_targets/processed_data/txt_files/ind_react_change/pm53_96.txt", header = FALSE)
pm54_96 <- read.delim("additional_targets/processed_data/txt_files/ind_react_change/pm54_96.txt", header = FALSE)
pm55_96 <- read.delim("additional_targets/processed_data/txt_files/ind_react_change/pm55_96.txt", header = FALSE)
pm60_96 <- read.delim("additional_targets/processed_data/txt_files/ind_react_change/pm60_96.txt", header = FALSE)
pm62_96 <- read.delim("additional_targets/processed_data/txt_files/ind_react_change/pm62_96.txt", header = FALSE)

pm15_96 <- read.delim("Big_exp/amp_change_r/pm15_96.txt", header = FALSE)
pm16_96 <- read.delim("Big_exp/amp_change_r/pm16_96.txt", header = FALSE)

pm11_96 <- read.delim("Big_exp/amp_change_r/pm11_96.txt", header = FALSE)
pm23_96 <- read.delim("Big_exp/amp_change_r/pm23_96.txt", header = FALSE)
pm42_96 <- read.delim("Big_exp/amp_change_r/pm42_96.txt", header = FALSE)

#create suitable data format
dat_pm46_96 <- pcrimport(file = "pm46_96", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)
dat_pm47_96 <- pcrimport(file = "pm47_96", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)
dat_pm48_96 <- pcrimport(file = "pm48_96", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)
dat_pm51_96 <- pcrimport(file = "pm51_96", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)
dat_pm52_96 <- pcrimport(file = "pm52_96", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)
dat_pm53_96 <- pcrimport(file = "pm53_96", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)
dat_pm54_96 <- pcrimport(file = "pm54_96", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)
dat_pm55_96 <- pcrimport(file = "pm55_96", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)
dat_pm60_96 <- pcrimport(file = "pm60_96", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)
dat_pm62_96 <- pcrimport(file = "pm62_96", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)

dat_pm15_96 <- pcrimport(file = "pm15_96", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)

dat_pm16_96 <- pcrimport(file = "pm16_96", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)

dat_pm11_96 <- pcrimport(file = "pm11_96", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)
dat_pm23_96 <- pcrimport(file = "pm23_96", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)
dat_pm42_96 <- pcrimport(file = "pm42_96", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)

## primer efficiencies ##

primer_eff_pm46 <- 1.953892
primer_eff_pm47 <- 1.904817
primer_eff_pm48 <- 2.012300827
primer_eff_pm51 <- 2.043823145
primer_eff_pm52 <- 1.981850457
primer_eff_pm53 <- 1.937201401
primer_eff_pm54 <- 1.875446
primer_eff_pm55 <- 1.901463
primer_eff_pm60 <- 1.873389731
primer_eff_pm62 <- 1.915484105

primer_eff_pm15 <- 1.89641
primer_eff_pm16 <- 2.01136

primer_eff_pm11 <- 1.91429
primer_eff_pm23 <- 1.91604674
primer_eff_pm42 <- 2.00500197


#### small version of masterfunction to get Cq-values #####

c46_96 <- ma_Cq(dat_pm46_96, primer_eff_pm46)
c47_96 <- ma_Cq(dat_pm47_96, primer_eff_pm47)
c48_96 <- ma_Cq(dat_pm48_96, primer_eff_pm48)
c51_96 <- ma_Cq(dat_pm51_96, primer_eff_pm51)
c52_96 <- ma_Cq(dat_pm52_96, primer_eff_pm52)
c53_96 <- ma_Cq(dat_pm53_96, primer_eff_pm53)
c54_96 <- ma_Cq(dat_pm54_96, primer_eff_pm54)
c55_96 <- ma_Cq(dat_pm55_96, primer_eff_pm55)
c60_96 <- ma_Cq(dat_pm60_96, primer_eff_pm60)
c62_96 <- ma_Cq(dat_pm62_96, primer_eff_pm62)

c11_96 <- ma_Cq(dat_pm11_96, primer_eff_pm11)
c23_96 <- ma_Cq(dat_pm23_96, primer_eff_pm23)
c42_96 <- ma_Cq(dat_pm42_96, primer_eff_pm42)

c46_96 <- gather(c46_96[,1:5], 'treatment', 'norm_exp', 1:5)
meanc46_96 <- mean(c46_96$norm_exp)
c47_96 <- gather(c47_96[,1:5], 'treatment', 'norm_exp', 1:5)
meanc47_96 <- mean(c47_96$norm_exp)
c48_96 <- gather(c48_96[,1:5], 'treatment', 'norm_exp', 1:5)
meanc48_96 <- mean(c48_96$norm_exp)
c51_96 <- gather(c51_96[,1:5], 'treatment', 'norm_exp', 1:5)
meanc51_96 <- mean(c51_96$norm_exp)
c52_96 <- gather(c52_96[,1:5], 'treatment', 'norm_exp', 1:5)
meanc52_96 <- mean(c52_96$norm_exp)
c53_96 <- gather(c53_96[,1:5], 'treatment', 'norm_exp', 1:5)
meanc53_96 <- mean(c53_96$norm_exp, na.rm = TRUE)
c54_96 <- gather(c54_96[,1:5], 'treatment', 'norm_exp', 1:5)
meanc54_96 <- mean(c54_96$norm_exp)
c55_96 <- gather(c55_96[,1:5], 'treatment', 'norm_exp', 1:5)
meanc55_96 <- mean(c55_96$norm_exp)
c60_96 <- gather(c60_96[,1:5], 'treatment', 'norm_exp', 1:5)
meanc60_96 <- mean(c60_96$norm_exp)
c62_96 <- gather(c62_96[,1:5], 'treatment', 'norm_exp', 1:5)
meanc62_96 <- mean(c62_96$norm_exp, na.rm = TRUE)

c11_96 <- gather(c11_96[,1:5], 'treatment', 'norm_exp', 1:5)
meanc11_96 <- mean(c11_96$norm_exp, na.rm = TRUE)
c23_96 <- gather(c23_96[,1:5], 'treatment', 'norm_exp', 1:5)
meanc23_96 <- mean(c23_96$norm_exp, na.rm = TRUE)
c42_96 <- gather(c42_96[,1:5], 'treatment', 'norm_exp', 1:5)
meanc42_96 <- mean(c42_96$norm_exp, na.rm = TRUE)

mean_vector <- c(meanc11_96, meanc23_96, meanc42_96, meanc46_96, meanc47_96, meanc48_96, meanc51_96, meanc52_96, meanc53_96, meanc54_96, meanc55_96, meanc60_96, meanc62_96)


######################################################
##################    masterfunction     #############
######################################################
m46_96 <- masterfunction(dat_pm46_96, primer_eff_pm46)
m47_96 <- masterfunction(dat_pm47_96, primer_eff_pm47)
m48_96 <- masterfunction(dat_pm48_96, primer_eff_pm48)
m51_96 <- masterfunction(dat_pm51_96, primer_eff_pm51)
m52_96 <- masterfunction(dat_pm52_96, primer_eff_pm52)
m53_96 <- masterfunction(dat_pm53_96, primer_eff_pm53)
m54_96 <- masterfunction(dat_pm54_96, primer_eff_pm54)
m55_96 <- masterfunction(dat_pm55_96, primer_eff_pm55)
m60_96 <- masterfunction(dat_pm60_96, primer_eff_pm60)
m62_96 <- masterfunction(dat_pm62_96, primer_eff_pm62)

m15_96 <- masterfunction(dat_pm15_96, primer_eff_pm15)

### the pm16_96 run had an issue in the last 3rd which was repeated and stitched together manually
setwd("~/Google Drive/qPCR")
m16_96 <- read.delim("Big_exp/m16_96_nfc.txt", header = TRUE)

m11_96 <- masterfunction(dat_pm11_96, primer_eff_pm11)
m23_96 <- masterfunction(dat_pm23_96, primer_eff_pm23)
m42_96 <- masterfunction(dat_pm42_96, primer_eff_pm42)

######################################################

ref_96 <- (m15_96 + m16_96)/2

############## normalising to ref ####################

nr46_96 <- m46_96/ref_96
nr47_96 <- m47_96/ref_96
nr48_96 <- m48_96/ref_96
nr51_96 <- m51_96/ref_96
nr52_96 <- m52_96/ref_96
nr53_96 <- m53_96/ref_96
nr54_96 <- m54_96/ref_96
nr55_96 <- m55_96/ref_96
nr60_96 <- m60_96/ref_96
nr62_96 <- m62_96/ref_96

nr11_96 <- m11_96/ref_96
nr23_96 <- m23_96/ref_96
nr42_96 <- m42_96/ref_96

####### creating a big dataframe for plotting ########

nr46_96$primer_pair <- 'VSP1'
nr47_96$primer_pair <- 'LOX2'
nr48_96$primer_pair <- 'ERF1'

nr51_96$primer_pair <- 'PAD3'
nr52_96$primer_pair <- 'CRK5'
nr53_96$primer_pair <- 'CRK23'
nr54_96$primer_pair <- 'AT1G51890'
nr55_96$primer_pair <- 'SWEET4'
nr60_96$primer_pair <- 'SWEET12'
nr62_96$primer_pair <- 'STP13'

nr11_96$primer_pair <- 'ETmarker'
nr23_96$primer_pair <- 'PR1'
nr42_96$primer_pair <- 'JAmarker'

marker_set <- rbind(nr46_96, nr47_96, nr48_96, nr51_96, nr52_96, nr53_96, nr54_96, nr55_96, nr60_96, nr62_96, nr11_96, nr23_96, nr42_96)

marker_set <- gather(marker_set, 'treatment', 'ref_norm', 1:5)
mock_set <- subset(marker_set, marker_set$treatment == 'mock_RT')

############################################
#### graphing log2 fold change to mock #####
############################################
marker_set[marker_set == 0] <- NA
marker_set[is.nan(marker_set)] <- NA
marker_set[is.infinite(marker_set)] <- NA
## The is.finite works on vector and not on data.frame object. So, we can loop through the data.frame using lapply and get only the 'finite' values.
#marker_set[Reduce(`&`, lapply(marker_set, is.finite)),] 
#doesn't really work

marker_set <- marker_set %>% 
  group_by(treatment, primer_pair) %>% 
  mutate(mean_ref_norm = mean(ref_norm, na.rm = TRUE, nan.rm = TRUE, infinite.rm = TRUE))

mock_set <- mock_set %>% 
  group_by(treatment, primer_pair) %>% 
  mutate(mean_ref_norm = mean(ref_norm, na.rm = TRUE, nan.rm = TRUE, infinite.rm = TRUE))

###################################
### log2 of fold change to ref ###
###################################
marker_set$log2_ref_norm <- log2(marker_set$ref_norm)

marker_set <- marker_set %>% 
  group_by(treatment, primer_pair) %>% 
  mutate(mean_log2_ref_norm = mean(log2_ref_norm, na.rm = TRUE, nan.rm = TRUE, infinite.rm = TRUE),
         se_log2_ref_norm = sd(log2_ref_norm, na.rm = TRUE)/sqrt(n()))

#####################################
### add biological replicates to data
#####################################
marker_set <- marker_set %>% ungroup() %>% add_count()
marker_set$biol_rep <- rep(1:3, times = marker_set$n[1]/3)

### change observation names in marker_set

marker_set$treatment <- gsub("mock_RT", "mock", marker_set$treatment)
marker_set$treatment <- gsub("Microbacterium_leaf347_RT", "Micro347", marker_set$treatment)
marker_set$treatment <- gsub("Pseudomonas_syringae_RT", "Pst", marker_set$treatment)
marker_set$treatment <- gsub("Sphingomonas_leaf34_RT", "Sphingo34", marker_set$treatment)
marker_set$treatment <- gsub("Williamsia_leaf354_RT", "Willi354", marker_set$treatment)

marker_set$treatment <- factor(marker_set$treatment, levels = c("mock", "Sphingo34", "Micro347", "Willi354", "Pst"))

##############################
### normalise against mock ###
##############################
marker_set <- marker_set %>%
  dplyr::select(primer_pair, treatment, mean_ref_norm, biol_rep) %>% 
  group_by(primer_pair) %>% 
  filter(treatment == 'mock' & biol_rep == 1) %>%
  right_join(marker_set, ., by =c('primer_pair'))


marker_set <-marker_set[,c(-10, -12)] 

names(marker_set)[10] <- 'mean_ref_norm_mock'
names(marker_set)[9] <- 'biol_rep'
names(marker_set)[4] <- 'mean_ref_norm'
names(marker_set)[2] <- 'treatment'

marker_set$mock_norm <- marker_set$ref_norm/marker_set$mean_ref_norm_mock

###################################
### log2 of fold change to mock ###
###################################
marker_set$log2_mock_norm <- log2(marker_set$mock_norm)

marker_set <- marker_set %>% 
  group_by(treatment, primer_pair) %>% 
  mutate(mean_log2_mock_norm = mean(log2_mock_norm, na.rm = TRUE, nan.rm = TRUE),
         se_log2_mock_norm = sd(log2_mock_norm, na.rm = TRUE)/sqrt(n()))


### which targets correlate in their gene expression profiles?

correlation_set <- marker_set %>% select(primer_pair, treatment, biol_rep, log2_ref_norm, log2_mock_norm)

correlation_set_ref <- correlation_set %>%
  select(-log2_mock_norm) %>%
  group_by(biol_rep, treatment) %>% 
  spread(key=primer_pair, value = log2_ref_norm)

ggplot(correlation_set_ref, aes(x=ETmarker, y=LOX2))+
  geom_point(size = 2.5)+
  theme_bw()+
  geom_smooth(method = 'lm')

library(nlme)
lm1 = lm(formula = STP13 ~ ETmarker, data = correlation_set_ref)
summary(lm1)
summary(lm1)$adj.r.squared # returns the adjusted r square
summary(lm1)$coefficients[2,4] # returns the p-value
summary(lm1)$terms[[2]] # returns first target gene
summary(lm1)$terms[[3]] # returns second target gene

str(summary(lm1)) # helps to find out how to call individual parts of the summary function

lm2 = lm(formula = STP13 ~ LOX2, data = correlation_set_ref)
summary(lm2)

lm3 = lm(formula = VSP1 ~ PAD3, data = correlation_set_ref)
summary(lm3)

lm4 = lm(formula = ETmarker ~ LOX2, data = correlation_set_ref)
summary(lm4)

lm.list <- list(lm1, lm2, lm3, lm4)
lm.frame <- data.frame()

for (i in 1:4){
  lm.frame[1,i] <- as.character(lm.list[[i]]$terms[[2]])
  lm.frame[2,i] <- as.character(lm.list[[i]]$terms[[3]])
  lm.frame[3,i] <- summary(lm.list[[i]])$adj.r.squared
  lm.frame[4,i] <- summary(lm.list[[i]])$coefficients[2,4]
}

####### now lets correlate each target gene against itself and each other target gene
library(corrplot)
library(RColorBrewer)
# change order of target genes
level_order <- c('PR1', 'VSP1', 'JAmarker', 'LOX2', 'ETmarker', 'ERF1', 'CRK5', 'CRK23', 'AT1G51890', 'SWEET4', 'SWEET12', 'STP13', 'PAD3') #this vector might be useful for other plots/analyses

counter <- length(level_order)
correlation_set2_ref <- correlation_set_ref[,3:15]
correlation_set2_ref <- subset(correlation_set2_ref, select = level_order) # changes order of df
df.correlation.set2_ref <- as.data.frame(correlation_set2_ref)

#create a df containing all adj. Rsquares
df.r.squared <-data.frame()

for (i in 1:counter){
  for(j in 1:counter){
    lm <- lm(formula = df.correlation.set2_ref[,i] ~ df.correlation.set2_ref[,j], data = df.correlation.set2_ref)
    df.r.squared[i,j] <- summary(lm)$adj.r.squared
  }
}
row.names(df.r.squared) <- level_order
colnames(df.r.squared) <- level_order
setwd("~/Google Drive/qPCR/Big_exp")

write.csv(df.r.squared, file = "correlation.matrix.r.squared.ref.csv", sep="\t",row.names=T, col.names = T)

#create a df containing all p.values
df.p.value <-data.frame()

for (i in 1:counter){
  for(j in 1:counter){
    lm <- lm(formula = df.correlation.set2_ref[,i] ~ df.correlation.set2_ref[,j], data = df.correlation.set2_ref)
    df.p.value[i,j] <- summary(lm)$coefficient[2,4]
  }
}
row.names(df.p.value) <- level_order
colnames(df.p.value) <- level_order
setwd("~/Google Drive/qPCR/Big_exp")

write.csv(df.p.value, file = "correlation.matrix.p.value.ref.csv", sep="\t",row.names=T, col.names = T)

# create a df containing all slopes
df.slope <-data.frame()

for (i in 1:counter){
  for(j in 1:counter){
    lm <- lm(formula = df.correlation.set2_ref[,i] ~ df.correlation.set2_ref[,j], data = df.correlation.set2_ref)
    df.slope[i,j] <- summary(lm)$coefficient[2,1]
  }
}
row.names(df.slope) <- level_order
colnames(df.slope) <- level_order
setwd("~/Google Drive/qPCR/Big_exp")

write.csv(df.slope, file = "correlation.matrix.slope.ref.csv", sep="\t",row.names=T, col.names = T)

## plot the adj. R square values for each combination using 'corrplot' with asterisks depicting the p.value sig. levels
heat_order <- c('ETmarker', 'LOX2', 'STP13', 'ERF1', 'CRK5', 'JAmarker',
                'SWEET4', 'AT1G51890', 'PR1', 'SWEET12', 'CRK23',
                'PAD3', 'VSP1')
df.r.squared <- df.r.squared[heat_order, heat_order]
df.p.value <- df.p.value[heat_order, heat_order]

corrplot(as.matrix(df.r.squared), p.mat = as.matrix(df.p.value),
         insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = 1.7, pch.col = "grey61",
         method = "color", type = "upper",
         #order = "hclust", addrect = 2,
         col = brewer.pal(n = 10, name = "BrBG"),
         tl.col = "black")



#### plot correlations based on data normalised against mock
correlation_set_mock <- correlation_set %>%
  select(-log2_ref_norm) %>%
  group_by(biol_rep, treatment) %>% 
  spread(key=primer_pair, value = log2_mock_norm)

counter <- length(level_order)
correlation_set2_mock <- correlation_set_mock[,3:15]
correlation_set2_mock <- subset(correlation_set2_mock, select = level_order) # changes order of df
df.correlation.set2_mock <- as.data.frame(correlation_set2_mock)

#create a df containing all adj. Rsquares
df.r.squared.mock <-data.frame()

for (i in 1:counter){
  for(j in 1:counter){
    lm <- lm(formula = df.correlation.set2_mock[,i] ~ df.correlation.set2_mock[,j], data = df.correlation.set2_mock)
    df.r.squared.mock[i,j] <- summary(lm)$adj.r.squared
  }
}
row.names(df.r.squared.mock) <- level_order
colnames(df.r.squared.mock) <- level_order
setwd("~/Google Drive/qPCR/Big_exp")

write.csv(df.r.squared.mock, file = "correlation.matrix.r.squared.mock.csv", sep="\t",row.names=T, col.names = T)

#create a df containing all p.values
df.p.value.mock <-data.frame()

for (i in 1:counter){
  for(j in 1:counter){
    lm <- lm(formula = df.correlation.set2_mock[,i] ~ df.correlation.set2_mock[,j], data = df.correlation.set2_mock)
    df.p.value.mock[i,j] <- summary(lm)$coefficient[2,4]
  }
}
row.names(df.p.value.mock) <- level_order
colnames(df.p.value.mock) <- level_order
setwd("~/Google Drive/qPCR/Big_exp")

write.csv(df.p.value.mock, file = "correlation.matrix.p.value.mock.csv", sep="\t",row.names=T, col.names = T)

# create a df containing all slopes
df.slope.mock <-data.frame()

for (i in 1:counter){
  for(j in 1:counter){
    lm <- lm(formula = df.correlation.set2_mock[,i] ~ df.correlation.set2_mock[,j], data = df.correlation.set2_mock)
    df.slope.mock[i,j] <- summary(lm)$coefficient[2,1]
  }
}
row.names(df.slope.mock) <- level_order
colnames(df.slope.mock) <- level_order
setwd("~/Google Drive/qPCR/Big_exp")

write.csv(df.slope.mock, file = "correlation.matrix.slope.mock.csv", sep="\t",row.names=T, col.names = T)

## plot the adj. R square values for each combination using 'corrplot' with asterisks depicting the p.value sig. levels

corrplot(as.matrix(df.r.squared.mock), p.mat = as.matrix(df.p.value.mock), insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = 1.7, pch.col = "black",
         method = "color", type = "upper",
         #order = "hclust", addrect = 2,
         col = brewer.pal(n = 10, name = "PiYG"),
         tl.col = "black")

###### towards a heatmap
#### heatmap based on mean expression
heat5 <- marker_set %>%
  filter(treatment != 'Sphingo34') %>%
  dplyr::select(primer_pair, treatment, mean_log2_mock_norm) %>% distinct()
heat5 <- spread(heat5, treatment, mean_log2_mock_norm)
heat5 <- heat5 %>% remove_rownames %>% column_to_rownames(var="primer_pair")
heat5 <- as.matrix(heat5)

mean_expression <- marker_set %>% dplyr::select(primer_pair, treatment, mean_ref_norm) %>%
  group_by(primer_pair) %>% summarise(mean_expression = mean(mean_ref_norm))

se5 <- marker_set %>% 
  filter(treatment != 'Sphingo34') %>%
  dplyr::select(primer_pair, treatment, se_log2_mock_norm) %>% distinct()
se5 <- spread(se5, treatment, se_log2_mock_norm)
se5 <- se5 %>% remove_rownames %>% column_to_rownames((var="primer_pair"))
se5 <- as.matrix(se5)

col_ep <- colorRamp2(c(-2, 2), c("white", "darkblue"))
min(m$m)
max(m$m)
m <- as.matrix(m %>% remove_rownames %>% column_to_rownames(var="primer_pair"))
heat5 <- heat5[heat_order,]
m <- m[heat_order,]
se5 <- se5[heat_order,]
##
higherlevel_order <- c('ethylene', 'jasmonic acid', 'sugar transporter', 'ethylene', 'receptor-like kinase', 'jasmonic acid',
                       'sugar transporter', 'receptor-like kinase', 'salicylic acid', 'sugar transporter', 'receptor-like kinase',
                       'secondary metabolism', 'jasmonic acid')
hl_order <- as.data.frame(higherlevel_order)
row.names(hl_order) <- heat_order
#
#global heatmap setting
ht_opt(legend_border = 'black',
       heatmap_border = TRUE,
       annotation_border = TRUE)

Heatmap(heat5, cluster_columns = FALSE, row_km = 3, row_km_repeats = 1000, border = TRUE, name = 'log2(fold change)',
        rect_gp = gpar(col = "black", lwd = 0.69), 
        right_annotation = rowAnnotation('log2(norm. expression)' = m,
                                         'id' = as.matrix(hl_order),
                                         gp = gpar(col = "black"),
                                         col = list('log2(norm. expression)' = col_ep,
                                                    'id' = c('ethylene' = '#fed439ff', 'jasmonic acid' = '#709ae1ff',
                                                             'salicylic acid' = '#fd7466ff', 'sugar transporter' = '#fd8cc1ff',
                                                             'receptor-like kinase' = '#8a9197ff', 'secondary metabolism' = '#1a9993ff')),
                                         show_annotation_name = FALSE,
                                         gap = unit(1, "mm")),
        row_gap = unit(3, "mm"),
        width = ncol(heat5)*unit(10, "mm"), 
        height = nrow(heat5)*unit(10, "mm"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", se5[i, j]), x, y, gp = gpar(fontsize = se5[i,j]*5+5))
        })

Heatmap(heat5, cluster_columns = FALSE, row_km = 2, row_km_repeats = 1000, border = TRUE, name = 'log2(fold change)',
        rect_gp = gpar(col = "black", lwd = 0.69), 
        right_annotation = rowAnnotation(
                                         'id' = as.matrix(hl_order),
                                         gp = gpar(col = "black"),
                                         col = list(
                                                    'id' = c('ethylene' = '#fed439ff', 'jasmonic acid' = '#709ae1ff',
                                                             'salicylic acid' = '#fd7466ff', 'sugar transporter' = '#fd8cc1ff',
                                                             'receptor-like kinase' = '#8a9197ff', 'secondary metabolism' = '#1a9993ff')),
                                         show_annotation_name = FALSE,
                                         gap = unit(1, "mm")),
        row_gap = unit(3, "mm"),
        width = ncol(heat5)*unit(10, "mm"), 
        height = nrow(heat5)*unit(10, "mm"),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%.1f", se5[i, j]), x, y, gp = gpar(fontsize = se5[i,j]*5+5))
        })






#### correlation plot only on mean
#create a df containing all adj. Rsquares
df.r.squared.heat5 <-data.frame()
theat5 <- as.data.frame(t(as.data.frame(heat5)))
level_order <- colnames(theat5)
for (i in 1:counter){
  for(j in 1:counter){
    lm <- lm(formula = theat5[,i] ~ theat5[,j], data = theat5)
    df.r.squared.heat5[i,j] <- summary(lm)$adj.r.squared
  }
}
row.names(df.r.squared.heat5) <- level_order
colnames(df.r.squared.heat5) <- level_order

#create a df containing all p.values
df.p.value.heat5 <-data.frame()

for (i in 1:counter){
  for(j in 1:counter){
    lm <- lm(formula = theat5[,i] ~ theat5[,j], data = theat5)
    df.p.value.heat5[i,j] <- summary(lm)$coefficient[2,4]
  }
}
row.names(df.p.value.heat5) <- level_order
colnames(df.p.value.heat5) <- level_order

heat_order <- c('ETmarker', 'LOX2', 'STP13', 'ERF1', 'CRK5', 'JAmarker',
                'SWEET4', 'AT1G51890', 'PR1', 'SWEET12', 'CRK23',
                'PAD3', 'VSP1')
df.r.squared.heat5 <- df.r.squared.heat5[heat_order, heat_order]
df.p.value.heat5 <- df.p.value.heat5[heat_order, heat_order]

corrplot(as.matrix(df.r.squared.heat5), p.mat = as.matrix(df.p.value.heat5),
         insig = "label_sig", sig.level = c(.001, .01, .05), pch.cex = 1.7, pch.col = "grey74",
         method = "color", type = "upper",
         #order = "hclust", addrect = 2,
         col = brewer.pal(n = 10, name = "BrBG"),
         tl.col = "black")





