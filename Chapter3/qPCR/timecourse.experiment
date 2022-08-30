#### Timecourse experiment

### This file is to calculate relative expression values from qPCR runs.
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
library(ggpubr)
library(dunn.test)

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
  mock_RT <- c(mean(c(x[1,1], x[1,2], x[1,3])), 
               mean(c(x[1,9], x[1,10], x[1,11])),
               mean(c(x[1,17], x[1,18], x[1,19])))
  
  mock_noRT <- c(x[1,4], x[1,12], x[1,20])
  
  Sphingomonas_leaf34_RT <- c(mean(c(x[1,25], x[1,26], x[1,27])), 
                              mean(c(x[1,33], x[1,34], x[1,35])),
                              mean(c(x[1,41], x[1,42], x[1,43])))
  
  Sphingomonas_leaf34_noRT <- c(x[1,28], x[1,36], x[1,44])
  
  Microbacterium_leaf347_RT <- c(mean(c(x[1,49], x[1,50], x[1,51])), 
                                 mean(c(x[1,57], x[1,58], x[1,59])),
                                 mean(c(x[1,65], x[1,66], x[1,67])))
  
  Microbacterium_leaf347_noRT <- c(x[1,52], x[1,60], x[1,68])
  
  Williamsia_leaf354_RT <- c(mean(c(x[1,5], x[1,6], x[1,7])), 
                             mean(c(x[1,13], x[1,14], x[1,15])),
                             mean(c(x[1,21], x[1,22], x[1,23])))
  
  Williamsia_leaf354_noRT <- c(x[1,8], x[1,16], x[1,24])
  
  Pseudomonas_syringae_RT <- c(mean(c(x[1,29], x[1,30], x[1,31])), 
                               mean(c(x[1,37], x[1,38], x[1,39])),
                               mean(c(x[1,45], x[1,46], x[1,47])))
  
  Pseudomonas_syringae_noRT <- c(x[1,32], x[1,40], x[1,48])
  
  interrun_RT <- c(mean(c(x[1,53], x[1,54], x[1,55], x[1,61], x[1,62], x[1,63]), na.rm = TRUE), NA, NA) ## in run ja09 some interrun calibrator samples ended up as 'NA'. the remaining ones seem very consitent though --> na.rm
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
# loading and processing data #
######################################################################################

# load data into R
setwd("~/Google Drive/qPCR/Big_exp")
setwd("/Volumes/GoogleDrive/My Drive/qPCR/Big_exp")

pm11_0 <- read.delim("amp_change_r/pm11_0.txt", header = FALSE)
pm11_3 <- read.delim("amp_change_r/pm11_3.txt", header = FALSE)
pm11_6 <- read.delim("amp_change_r/pm11_6.txt", header = FALSE)
pm11_9 <- read.delim("amp_change_r/pm11_9.txt", header = FALSE)
pm11_12 <- read.delim("amp_change_r/pm11_12.txt", header = FALSE)
pm11_24 <- read.delim("amp_change_r/pm11_24.txt", header = FALSE)
pm11_48 <- read.delim("amp_change_r/pm11_48.txt", header = FALSE)
pm11_96 <- read.delim("amp_change_r/pm11_96.txt", header = FALSE)

pm15_0 <- read.delim("amp_change_r/pm15_0.txt", header = FALSE)
pm15_3 <- read.delim("amp_change_r/pm15_3.txt", header = FALSE)
pm15_6 <- read.delim("amp_change_r/pm15_6.txt", header = FALSE)
pm15_9 <- read.delim("amp_change_r/pm15_9.txt", header = FALSE)
pm15_12 <- read.delim("amp_change_r/pm15_12.txt", header = FALSE)
pm15_24 <- read.delim("amp_change_r/pm15_24.txt", header = FALSE)
pm15_48 <- read.delim("amp_change_r/pm15_48.txt", header = FALSE)
pm15_96 <- read.delim("amp_change_r/pm15_96.txt", header = FALSE)

pm16_0 <- read.delim("amp_change_r/pm16_0.txt", header = FALSE)
pm16_3 <- read.delim("amp_change_r/pm16_3.txt", header = FALSE)
pm16_6 <- read.delim("amp_change_r/pm16_6.txt", header = FALSE)
pm16_9 <- read.delim("amp_change_r/pm16_9.txt", header = FALSE)
pm16_12 <- read.delim("amp_change_r/pm16_12.txt", header = FALSE)
pm16_24 <- read.delim("amp_change_r/pm16_24.txt", header = FALSE)
pm16_48 <- read.delim("amp_change_r/pm16_48.txt", header = FALSE)
pm16_96 <- read.delim("amp_change_r/pm16_96.txt", header = FALSE)

pm23_0 <- read.delim("amp_change_r/pm23_0.txt", header = FALSE)
pm23_3 <- read.delim("amp_change_r/pm23_3.txt", header = FALSE)
pm23_6 <- read.delim("amp_change_r/pm23_6.txt", header = FALSE)
pm23_9 <- read.delim("amp_change_r/pm23_9.txt", header = FALSE)
pm23_12 <- read.delim("amp_change_r/pm23_12.txt", header = FALSE)
pm23_24 <- read.delim("amp_change_r/pm23_24.txt", header = FALSE)
pm23_48 <- read.delim("amp_change_r/pm23_48.txt", header = FALSE)
pm23_96 <- read.delim("amp_change_r/pm23_96.txt", header = FALSE)

pm34_0 <- read.delim("amp_change_r/pm34_0.txt", header = FALSE)
pm34_12 <- read.delim("amp_change_r/pm34_12.txt", header = FALSE)

pm42_3 <- read.delim("amp_change_r/pm42_3.txt", header = FALSE)
pm42_6 <- read.delim("amp_change_r/pm42_6.txt", header = FALSE)
pm42_9 <- read.delim("amp_change_r/pm42_9.txt", header = FALSE)
pm42_24 <- read.delim("amp_change_r/pm42_24.txt", header = FALSE)
pm42_48 <- read.delim("amp_change_r/pm42_48.txt", header = FALSE)
pm42_96 <- read.delim("amp_change_r/pm42_96.txt", header = FALSE)

#create suitable data format
dat_pm11_0 <- pcrimport(file = "pm11_0", sep = "\t", dec = ".", delCol = 0,
                        delRow = 0, format = "col", sampleDat = 0, names = 0,
                        sampleLen = 1:40, check = FALSE)
dat_pm11_3 <- pcrimport(file = "pm11_3", sep = "\t", dec = ".", delCol = 0,
                        delRow = 0, format = "col", sampleDat = 0, names = 0,
                        sampleLen = 1:40, check = FALSE)
dat_pm11_6 <- pcrimport(file = "pm11_6", sep = "\t", dec = ".", delCol = 0,
                        delRow = 0, format = "col", sampleDat = 0, names = 0,
                        sampleLen = 1:40, check = FALSE)
dat_pm11_9 <- pcrimport(file = "pm11_9", sep = "\t", dec = ".", delCol = 0,
                        delRow = 0, format = "col", sampleDat = 0, names = 0,
                        sampleLen = 1:40, check = FALSE)
dat_pm11_12 <- pcrimport(file = "pm11_12", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)
dat_pm11_24 <- pcrimport(file = "pm11_24", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)
dat_pm11_48 <- pcrimport(file = "pm11_48", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)
dat_pm11_96 <- pcrimport(file = "pm11_96", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)

dat_pm15_0 <- pcrimport(file = "pm15_0", sep = "\t", dec = ".", delCol = 0,
                        delRow = 0, format = "col", sampleDat = 0, names = 0,
                        sampleLen = 1:40, check = FALSE)
dat_pm15_3 <- pcrimport(file = "pm15_3", sep = "\t", dec = ".", delCol = 0,
                        delRow = 0, format = "col", sampleDat = 0, names = 0,
                        sampleLen = 1:40, check = FALSE)
dat_pm15_6 <- pcrimport(file = "pm15_6", sep = "\t", dec = ".", delCol = 0,
                        delRow = 0, format = "col", sampleDat = 0, names = 0,
                        sampleLen = 1:40, check = FALSE)
dat_pm15_9 <- pcrimport(file = "pm15_9", sep = "\t", dec = ".", delCol = 0,
                        delRow = 0, format = "col", sampleDat = 0, names = 0,
                        sampleLen = 1:40, check = FALSE)
dat_pm15_12 <- pcrimport(file = "pm15_12", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)
dat_pm15_24 <- pcrimport(file = "pm15_24", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)
dat_pm15_48 <- pcrimport(file = "pm15_48", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)
dat_pm15_96 <- pcrimport(file = "pm15_96", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)

dat_pm16_0 <- pcrimport(file = "pm16_0", sep = "\t", dec = ".", delCol = 0,
                        delRow = 0, format = "col", sampleDat = 0, names = 0,
                        sampleLen = 1:40, check = FALSE)
dat_pm16_3 <- pcrimport(file = "pm16_3", sep = "\t", dec = ".", delCol = 0,
                        delRow = 0, format = "col", sampleDat = 0, names = 0,
                        sampleLen = 1:40, check = FALSE)
dat_pm16_6 <- pcrimport(file = "pm16_6", sep = "\t", dec = ".", delCol = 0,
                        delRow = 0, format = "col", sampleDat = 0, names = 0,
                        sampleLen = 1:40, check = FALSE)
dat_pm16_9 <- pcrimport(file = "pm16_9", sep = "\t", dec = ".", delCol = 0,
                        delRow = 0, format = "col", sampleDat = 0, names = 0,
                        sampleLen = 1:40, check = FALSE)
dat_pm16_12 <- pcrimport(file = "pm16_12", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)
dat_pm16_24 <- pcrimport(file = "pm16_24", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)
dat_pm16_48 <- pcrimport(file = "pm16_48", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)
dat_pm16_96 <- pcrimport(file = "pm16_96", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)

dat_pm23_0 <- pcrimport(file = "pm23_0", sep = "\t", dec = ".", delCol = 0,
                        delRow = 0, format = "col", sampleDat = 0, names = 0,
                        sampleLen = 1:40, check = FALSE)
dat_pm23_3 <- pcrimport(file = "pm23_3", sep = "\t", dec = ".", delCol = 0,
                        delRow = 0, format = "col", sampleDat = 0, names = 0,
                        sampleLen = 1:40, check = FALSE)
dat_pm23_6 <- pcrimport(file = "pm23_6", sep = "\t", dec = ".", delCol = 0,
                        delRow = 0, format = "col", sampleDat = 0, names = 0,
                        sampleLen = 1:40, check = FALSE)
dat_pm23_9 <- pcrimport(file = "pm23_9", sep = "\t", dec = ".", delCol = 0,
                        delRow = 0, format = "col", sampleDat = 0, names = 0,
                        sampleLen = 1:40, check = FALSE)
dat_pm23_12 <- pcrimport(file = "pm23_12", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)
dat_pm23_24 <- pcrimport(file = "pm23_24", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)
dat_pm23_48 <- pcrimport(file = "pm23_48", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)
dat_pm23_96 <- pcrimport(file = "pm23_96", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)

dat_pm42_3 <- pcrimport(file = "pm42_3", sep = "\t", dec = ".", delCol = 0,
                        delRow = 0, format = "col", sampleDat = 0, names = 0,
                        sampleLen = 1:40, check = FALSE)
dat_pm42_6 <- pcrimport(file = "pm42_6", sep = "\t", dec = ".", delCol = 0,
                        delRow = 0, format = "col", sampleDat = 0, names = 0,
                        sampleLen = 1:40, check = FALSE)
dat_pm42_9 <- pcrimport(file = "pm42_9", sep = "\t", dec = ".", delCol = 0,
                        delRow = 0, format = "col", sampleDat = 0, names = 0,
                        sampleLen = 1:40, check = FALSE)
dat_pm42_24 <- pcrimport(file = "pm42_24", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)
dat_pm42_48 <- pcrimport(file = "pm42_48", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)
dat_pm42_96 <- pcrimport(file = "pm42_96", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)

dat_pm34_0 <- pcrimport(file = "pm34_0", sep = "\t", dec = ".", delCol = 0,
                        delRow = 0, format = "col", sampleDat = 0, names = 0,
                        sampleLen = 1:40, check = FALSE)
dat_pm34_12 <- pcrimport(file = "pm34_12", sep = "\t", dec = ".", delCol = 0,
                         delRow = 0, format = "col", sampleDat = 0, names = 0,
                         sampleLen = 1:40, check = FALSE)

## primer efficiencies

primer_eff_pm11 <- 1.91429
primer_eff_pm15 <- 1.89641
primer_eff_pm16 <- 2.01136
primer_eff_pm23 <- 1.91604674
primer_eff_pm34 <- 1.87531757
primer_eff_pm42 <- 2.00500197



####################################################
m11_0 <- masterfunction(dat_pm11_0, primer_eff_pm11)
m11_3 <- masterfunction(dat_pm11_3, primer_eff_pm11)
m11_6 <- masterfunction(dat_pm11_6, primer_eff_pm11)
m11_9 <- masterfunction(dat_pm11_9, primer_eff_pm11)
m11_12 <- masterfunction(dat_pm11_12, primer_eff_pm11)
m11_24 <- masterfunction(dat_pm11_24, primer_eff_pm11)
m11_48 <- masterfunction(dat_pm11_48, primer_eff_pm11)
m11_96 <- masterfunction(dat_pm11_96, primer_eff_pm11)

m15_0 <- masterfunction(dat_pm15_0, primer_eff_pm15)
m15_3 <- masterfunction(dat_pm15_3, primer_eff_pm15)
m15_6 <- masterfunction(dat_pm15_6, primer_eff_pm15)
m15_9 <- masterfunction(dat_pm15_9, primer_eff_pm15)
m15_12 <- masterfunction(dat_pm15_12, primer_eff_pm15)
m15_24 <- masterfunction(dat_pm15_24, primer_eff_pm15)
m15_48 <- masterfunction(dat_pm15_48, primer_eff_pm15)
m15_96 <- masterfunction(dat_pm15_96, primer_eff_pm15)

m16_0 <- masterfunction(dat_pm16_0, primer_eff_pm16)
m16_3 <- masterfunction(dat_pm16_3, primer_eff_pm16)
m16_6 <- masterfunction(dat_pm16_6, primer_eff_pm16)
m16_9 <- masterfunction(dat_pm16_9, primer_eff_pm16)
m16_12 <- masterfunction(dat_pm16_12, primer_eff_pm16)
m16_24 <- masterfunction(dat_pm16_24, primer_eff_pm16)
m16_48 <- masterfunction(dat_pm16_48, primer_eff_pm16)
#m16_96 <- masterfunction(dat_pm16_96, primer_eff_pm16) ### last3rd was repeated as there were problems in that run
##setwd("~/Desktop/qPCR/a_qPCR_plant_response_timecourse/txt_files")
##write.table(m16_96, file = "m16_96_nfc.txt", sep="\t",row.names=FALSE, col.names = TRUE) ### nfc = need for correction
#above was used to write a table. then the values for Microbacterium were edited by hand. see file 'qPCR_June2020_cpD2_pm16_96_last3rd'
setwd("~/Google Drive/qPCR")
m16_96 <- read.delim("Big_exp/m16_96_nfc.txt", header = TRUE)

m23_0 <- masterfunction(dat_pm23_0, primer_eff_pm23)
m23_3 <- masterfunction(dat_pm23_3, primer_eff_pm23)
m23_6 <- masterfunction(dat_pm23_6, primer_eff_pm23)
m23_9 <- masterfunction(dat_pm23_9, primer_eff_pm23)
m23_12 <- masterfunction(dat_pm23_12, primer_eff_pm23)
m23_24 <- masterfunction(dat_pm23_24, primer_eff_pm23)
m23_48 <- masterfunction(dat_pm23_48, primer_eff_pm23)
m23_96 <- masterfunction(dat_pm23_96, primer_eff_pm23)

m42_3 <- masterfunction(dat_pm42_3, primer_eff_pm42)
m42_6 <- masterfunction(dat_pm42_6, primer_eff_pm42)
m42_9 <- masterfunction(dat_pm42_9, primer_eff_pm42)
m42_24 <- masterfunction(dat_pm42_24, primer_eff_pm42)
m42_48 <- masterfunction(dat_pm42_48, primer_eff_pm42)
m42_96 <- masterfunction(dat_pm42_96, primer_eff_pm42)

m34_0 <- masterfunction(dat_pm34_0, primer_eff_pm34)
m34_12 <- masterfunction(dat_pm34_12, primer_eff_pm34)
######################################################

ref_0 <- (m15_0 + m16_0)/2
ref_3 <- (m15_3 + m16_3)/2
ref_6 <- (m15_6 + m16_6)/2
ref_9 <- (m15_9 + m16_9)/2
ref_12 <- (m15_12 + m16_12)/2
ref_24 <- (m15_24 + m16_24)/2
ref_48 <- (m15_48 + m16_48)/2
ref_96 <- (m15_96 + m16_96)/2

######################################################
setwd("~/Google Drive/qPCR/Big_exp/result_norm_ref")

write.table(m11_0/ref_0, file = "pm11_0_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(m23_0/ref_0, file = "pm23_0_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(m34_0/ref_0, file = "pm34_0_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)

write.table(m11_3/ref_3, file = "pm11_3_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(m23_3/ref_3, file = "pm23_3_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(m42_3/ref_3, file = "pm42_3_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)

write.table(m11_6/ref_6, file = "pm11_6_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(m23_6/ref_6, file = "pm23_6_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(m42_6/ref_6, file = "pm42_6_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)

write.table(m11_9/ref_9, file = "pm11_9_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(m23_9/ref_9, file = "pm23_9_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(m42_9/ref_9, file = "pm42_9_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)

write.table(m11_12/ref_12, file = "pm11_12_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(m23_12/ref_12, file = "pm23_12_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(m34_12/ref_12, file = "pm34_12_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)

write.table(m11_24/ref_24, file = "pm11_24_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(m23_24/ref_24, file = "pm23_24_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(m42_24/ref_24, file = "pm42_24_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)

write.table(m11_48/ref_48, file = "pm11_48_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(m23_48/ref_48, file = "pm23_48_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(m42_48/ref_48, file = "pm42_48_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)

write.table(m11_96/ref_96, file = "pm11_96_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(m23_96/ref_96, file = "pm23_96_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(m42_96/ref_96, file = "pm42_96_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)

#################################################################################################

setwd("~/Google Drive/qPCR/Big_exp/result_norm_ref_mock")

####### normalise to ref and then to mock (-> nfunc 'normalising function') 
nfunc <- function(x, y){
  nr <- x/y
  nr[nr < 0.000000001] <- NA
  mean_mock <- mean(c(nr[1, 1], nr[2, 1], nr[3, 1]), na.rm = TRUE)
  nr/mean_mock
}


write.table(nfunc(m11_0, ref_0), file = "pm11_0_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(nfunc(m23_0, ref_0), file = "pm23_0_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(nfunc(m34_0, ref_0), file = "pm34_0_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)

write.table(nfunc(m11_3, ref_3), file = "pm11_3_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(nfunc(m23_3, ref_3), file = "pm23_3_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(nfunc(m42_3, ref_3), file = "pm42_3_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)

write.table(nfunc(m11_6, ref_6), file = "pm11_6_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(nfunc(m23_6, ref_6), file = "pm23_6_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(nfunc(m42_6, ref_6), file = "pm42_6_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)

write.table(nfunc(m11_9, ref_9), file = "pm11_9_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(nfunc(m23_9, ref_9), file = "pm23_9_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(nfunc(m42_9, ref_9), file = "pm42_9_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)

write.table(nfunc(m11_12, ref_12), file = "pm11_12_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(nfunc(m23_12, ref_12), file = "pm23_12_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(nfunc(m34_12, ref_12), file = "pm34_12_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)

write.table(nfunc(m11_24, ref_24), file = "pm11_24_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(nfunc(m23_24, ref_24), file = "pm23_24_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(nfunc(m42_24, ref_24), file = "pm42_24_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)

write.table(nfunc(m11_48, ref_48), file = "pm11_48_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(nfunc(m23_48, ref_48), file = "pm23_48_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(nfunc(m42_48, ref_48), file = "pm42_48_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)

write.table(nfunc(m11_96, ref_96), file = "pm11_96_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(nfunc(m23_96, ref_96), file = "pm23_96_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)
write.table(nfunc(m42_96, ref_96), file = "pm42_96_nref.txt", sep="\t",row.names=FALSE, col.names = TRUE)

##################################################
#######           graphing ref genes      ########
##################################################
m15_0$timepoint <- '01'
m15_3$timepoint <- '03'
m15_6$timepoint <- '06'
m15_9$timepoint <- '09'
m15_12$timepoint <- '12'
m15_24$timepoint <- '24'
m15_48$timepoint <- '48'
m15_96$timepoint <- '96'

m16_0$timepoint <- '01'
m16_3$timepoint <- '03'
m16_6$timepoint <- '06'
m16_9$timepoint <- '09'
m16_12$timepoint <- '12'
m16_24$timepoint <- '24'
m16_48$timepoint <- '48'
m16_96$timepoint <- '96'

m15_0$time <- 1
m15_3$time <- 3
m15_6$time <- 6
m15_9$time <- 9
m15_12$time <- 12
m15_24$time <- 24
m15_48$time <- 48
m15_96$time <- 96

m16_0$time <- 1
m16_3$time <- 3
m16_6$time <- 6
m16_9$time <- 9
m16_12$time <- 12
m16_24$time <- 24
m16_48$time <- 48
m16_96$time <- 96

m15_0$primer_pair <- 'ref15'
m15_3$primer_pair <- 'ref15'
m15_6$primer_pair <- 'ref15'
m15_9$primer_pair <- 'ref15'
m15_12$primer_pair <- 'ref15'
m15_24$primer_pair <- 'ref15'
m15_48$primer_pair <- 'ref15'
m15_96$primer_pair <- 'ref15'

m16_0$primer_pair <- 'ref16'
m16_3$primer_pair <- 'ref16'
m16_6$primer_pair <- 'ref16'
m16_9$primer_pair <- 'ref16'
m16_12$primer_pair <- 'ref16'
m16_24$primer_pair <- 'ref16'
m16_48$primer_pair <- 'ref16'
m16_96$primer_pair <- 'ref16'

ref_set <- rbind(m15_0, m15_3, m15_6, m15_9, m15_12, m15_24, m15_48, m15_96,
                 m16_0, m16_3, m16_6, m16_9, m16_12, m16_24, m16_48, m16_96)

ref_set <- gather(ref_set, 'treatment', 'norm_exp', 1:5)

theme_cu <- function() {
  theme_bw() +
    theme(panel.border = element_rect(colour = "grey54", linetype = 'solid', size = 1),
          panel.grid.major = element_line(colour = "grey54", linetype = 'solid', size = 0.5),
          panel.grid.minor = element_line(colour = "grey54", linetype = 'solid', size = 0.5),
          axis.text = element_text(colour = "black", family = "sans", size = 15),
          axis.title = element_text(colour = "black", family = "sans", size = 20),
          axis.ticks = element_blank(),
          # legend at the bottom 
          legend.justification = "center",
          legend.position = "bottom",
          legend.title = element_text(size =20),
          legend.text = element_text(colour = "black", family = "sans", size = 20),
          # plot title
          plot.title = element_text(colour = "black", family = "sans", size = 20),
          aspect.ratio = 1
    )
}

ggplot(data = ref_set, aes(timepoint, norm_exp)) +
  geom_boxplot(aes(fill = primer_pair)) +
  labs(title = 'reference gene expression',
       x = 'time (hpi)',
       y = 'norm. template concentration')+ # normalised to the interrun calibrator
  scale_fill_manual(values = c("#eccbae", "#abddde"))+
  theme_cu()

##########################################
### change observation names in ref_set
##########################################
ref_set$treatment <- gsub("mock_RT", "Mock", ref_set$treatment)
ref_set$treatment <- gsub("Sphingomonas_leaf34_RT", "Sphingo34", ref_set$treatment)
ref_set$treatment <- gsub("Microbacterium_leaf347_RT", "Micro347", ref_set$treatment)
ref_set$treatment <- gsub("Williamsia_leaf354_RT", "Willi354", ref_set$treatment)
ref_set$treatment <- gsub("Pseudomonas_syringae_RT", "Pst", ref_set$treatment)

ref_set$treatment <- factor(ref_set$treatment, levels = c("Mock", "Sphingo34", "Micro347", "Willi354", "Pst"))

m11_0$timepoint <- '01'
m11_3$timepoint <- '03'
m11_6$timepoint <- '06'
m11_9$timepoint <- '09'
m11_12$timepoint <- '12'
m11_24$timepoint <- '24'
m11_48$timepoint <- '48'
m11_96$timepoint <- '96'

m23_0$timepoint <- '01'
m23_3$timepoint <- '03'
m23_6$timepoint <- '06'
m23_9$timepoint <- '09'
m23_12$timepoint <- '12'
m23_24$timepoint <- '24'
m23_48$timepoint <- '48'
m23_96$timepoint <- '96'

m34_0$timepoint <- '01'
m42_3$timepoint <- '03'
m42_6$timepoint <- '06'
m42_9$timepoint <- '09'
m34_12$timepoint <- '12'
m42_24$timepoint <- '24'
m42_48$timepoint <- '48'
m42_96$timepoint <- '96'

m11_0$time <- 1
m11_3$time <- 3
m11_6$time <- 6 
m11_9$time <- 9
m11_12$time <- 12
m11_24$time <- 24
m11_48$time <- 48
m11_96$time <- 96

m23_0$time <- 1
m23_3$time <- 3
m23_6$time <- 6 
m23_9$time <- 9
m23_12$time <- 12
m23_24$time <- 24
m23_48$time <- 48
m23_96$time <- 96

m34_0$time <- 1
m42_3$time <- 3
m42_6$time <- 6 
m42_9$time <- 9
m34_12$time <- 12
m42_24$time <- 24
m42_48$time <- 48
m42_96$time <- 96

m11_0$primer_pair <- 'et.marker'
m11_3$primer_pair <- 'et.marker'
m11_6$primer_pair <- 'et.marker'
m11_9$primer_pair <- 'et.marker'
m11_12$primer_pair <- 'et.marker'
m11_24$primer_pair <- 'et.marker'
m11_48$primer_pair <- 'et.marker'
m11_96$primer_pair <- 'et.marker'

m23_0$primer_pair <- 'sa.marker'
m23_3$primer_pair <- 'sa.marker'
m23_6$primer_pair <- 'sa.marker'
m23_9$primer_pair <- 'sa.marker'
m23_12$primer_pair <- 'sa.marker'
m23_24$primer_pair <- 'sa.marker'
m23_48$primer_pair <- 'sa.marker'
m23_96$primer_pair <- 'sa.marker'

m34_0$primer_pair <- 'ja.marker'
m42_3$primer_pair <- 'ja.marker'
m42_6$primer_pair <- 'ja.marker'
m42_9$primer_pair <- 'ja.marker'
m34_12$primer_pair <- 'ja.marker'
m42_24$primer_pair <- 'ja.marker'
m42_48$primer_pair <- 'ja.marker'
m42_96$primer_pair <- 'ja.marker'



marker_set <- rbind(m15_0, m15_3, m15_6, m15_9, m15_12, m15_24, m15_48, m15_96,
                    m16_0, m16_3, m16_6, m16_9, m16_12, m16_24, m16_48, m16_96,
                    m11_0, m11_3, m11_6, m11_9, m11_12, m11_24, m11_48, m11_96,
                    m23_0, m23_3, m23_6, m23_9, m23_12, m23_24, m23_48, m23_96,
                    m34_0, m42_3, m42_6, m42_9, m34_12, m42_24, m42_48, m42_96)


marker_set <- gather(marker_set, 'treatment', 'norm_exp', 1:5)


##########################################
### change observation names in marker_set
##########################################
marker_set$treatment <- gsub("mock_RT", "Mock", marker_set$treatment)
marker_set$treatment <- gsub("Sphingomonas_leaf34_RT", "Sphingo34", marker_set$treatment)
marker_set$treatment <- gsub("Microbacterium_leaf347_RT", "Micro347", marker_set$treatment)
marker_set$treatment <- gsub("Williamsia_leaf354_RT", "Willi354", marker_set$treatment)
marker_set$treatment <- gsub("Pseudomonas_syringae_RT", "Pst", marker_set$treatment)

marker_set$treatment <- factor(marker_set$treatment, levels = c("Mock", "Sphingo34", "Micro347", "Willi354", "Pst"))

marker_set$norm_exp[marker_set$norm_exp == 0] <- NA

### normalise to reference genes

ref <- marker_set %>% filter(primer_pair %in% c('ref15', 'ref16')) %>%
  group_by(timepoint, treatment) %>% mutate(mean_ref = mean(norm_exp, na.rm = TRUE)) %>%
  select(-c(primer_pair, norm_exp, time))

marker_set <- right_join(marker_set, ref, by = c('timepoint', 'treatment')) %>% distinct() %>%
  mutate(ref_norm = norm_exp/mean_ref) %>% select(-c('norm_exp', 'mean_ref'))


marker_set <- marker_set %>% 
  group_by(treatment, primer_pair, timepoint) %>% 
  mutate(mean_ref_norm = mean(ref_norm, na.rm = TRUE, nan.rm = TRUE, infinite.rm = TRUE))

###################################
### log2 of fold change to ref ###
###################################
marker_set$log2_ref_norm <- log2(marker_set$ref_norm)

marker_set <- marker_set %>% 
  group_by(treatment, primer_pair, timepoint) %>% 
  mutate(mean_log2_ref_norm = mean(log2_ref_norm, na.rm = TRUE, nan.rm = TRUE, infinite.rm = TRUE),
         sd_log2_ref_norm = sd(log2_ref_norm, na.rm = TRUE))

#####################################
### add biological replicates to data
#####################################
marker_set <- marker_set %>% ungroup() %>% add_count()
marker_set$biol_rep <- rep(1:3, times = marker_set$n[1]/3)

##############################
### normalise against mock ###
##############################

mock <- marker_set %>% filter(treatment %in% 'Mock') %>%
  select(c('timepoint', 'primer_pair', 'mean_ref_norm')) %>%
  group_by(timepoint, primer_pair) %>%
  summarise(mean_ref_norm_mock = mean_ref_norm)

marker_set <- right_join(marker_set, mock, by = c('timepoint', 'primer_pair'))

marker_set$mock_norm <- marker_set$ref_norm/marker_set$mean_ref_norm_mock

###################################
### log2 of fold change to mock ###
###################################
marker_set$log2_mock_norm <- log2(marker_set$mock_norm)

marker_set <- marker_set %>% 
  group_by(treatment, primer_pair, timepoint) %>% 
  mutate(mean_log2_mock_norm = mean(log2_mock_norm, na.rm = TRUE, nan.rm = TRUE),
         sd_log2_mock_norm = sd(log2_mock_norm, na.rm = TRUE))

#############################################################################
### plotting changes in gene expression in relation to mock control ######### 
#############################################################################

# change order of target genes
level_order <- c('PR1', 'ja.marker', 'et.marker') #this vector might be useful for other plots/analyses

### ggplot figure in heatmap orientation

marker_gg <- marker_set %>% 
  dplyr::select(primer_pair, treatment, timepoint, mean_log2_mock_norm, sd_log2_mock_norm, log2_mock_norm, biol_rep) %>%
  filter(primer_pair != 'ref15') %>% filter(primer_pair != 'ref16') %>%
  distinct()

level_order <- c('Pst', 'Willi354', 'Micro347', 'Sphingo34', 'Mock')
max(marker_gg$mean_log2_mock_norm, na.rm = TRUE)
min(marker_gg$mean_log2_mock_norm, na.rm = TRUE)

marker_gg %>% filter(primer_pair %in% 'et.marker') %>%
  ggplot()+
  geom_vline(xintercept = 0, size = 1)+
  geom_hline(yintercept = 1.5, size = 1)+ geom_hline(yintercept = 2.5, size = 1)+ geom_hline(yintercept = 3.45, size = 1)+
  geom_hline(yintercept = 3.55, size = 1)+ geom_hline(yintercept = 4.5, size = 1)+ geom_hline(yintercept = 5.5, size = 1)+
  geom_hline(yintercept = 6.5, size = 1)+ geom_hline(yintercept = 7.5, size = 1)+
  geom_pointrange(aes(mean_log2_mock_norm, timepoint, shape = factor(treatment, level = level_order), fill = mean_log2_mock_norm,
                      xmin = mean_log2_mock_norm-sd_log2_mock_norm, xmax = mean_log2_mock_norm+sd_log2_mock_norm),
                  position = position_dodge(width = 0.66), size = 2)+
  geom_point(aes(log2_mock_norm, timepoint, shape = factor(treatment, level = level_order)),
             position = position_dodge(width = 0.66), size = 5, fill = 'grey24', alpha = 0.5)+
  scale_shape_manual(values = c(25,24,23,22,21), name = '')+
  scale_fill_gradient2(name = 'log2(fold change)', low = "blue", mid = "white", high = "red", midpoint = 0, limits=c(-2,1.5))+
  scale_y_discrete(limits=rev)+
  theme_cu()+
  theme(legend.position = 'right',
        aspect.ratio = 2.2,
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.minor.y = element_line(size = 2, colour = 'black'),
        legend.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        panel.border = element_rect(size = 2, colour = 'black'))+
  labs(x = '', y = '')


marker_gg %>% 
  filter(primer_pair %in% 'et.marker') %>%
  filter(treatment %in% 'Pst') %>%
  ggplot()+
  geom_vline(xintercept = 0, size = 1)+
  geom_hline(yintercept = 1.5, size = 1)+ geom_hline(yintercept = 2.5, size = 1)+ geom_hline(yintercept = 3.45, size = 1)+
  geom_hline(yintercept = 3.55, size = 1)+ geom_hline(yintercept = 4.5, size = 1)+ geom_hline(yintercept = 5.5, size = 1)+
  geom_hline(yintercept = 6.5, size = 1)+ geom_hline(yintercept = 7.5, size = 1)+
  geom_pointrange(aes(mean_log2_mock_norm, timepoint, fill = mean_log2_mock_norm,
                      xmin = mean_log2_mock_norm-sd_log2_mock_norm, xmax = mean_log2_mock_norm+sd_log2_mock_norm),
                  position = position_dodge(width = 0.66), size = 2, shape = 25)+
  geom_point(aes(log2_mock_norm, timepoint),
             position = position_dodge(width = 0.66), size = 5, fill = 'grey24', alpha = 0.5, shape = 25)+
  scale_fill_gradient2(name = 'log2(fold change)', low = "blue", mid = "white", high = "red", midpoint = 0, limits=c(-2,1.5))+
  scale_y_discrete(limits=rev)+
  theme_cu()+
  theme(legend.position = 'right',
        aspect.ratio = 1.5,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.minor.y = element_line(size = 2, colour = 'black'),
        legend.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        axis.text.x = element_blank(),
        panel.border = element_rect(size = 2, colour = 'black'))+
  labs(x = '', y = '')+
  xlim(-3,2)

marker_gg %>% 
  filter(primer_pair %in% 'et.marker') %>%
  filter(treatment %in% 'Sphingo34') %>%
  ggplot()+
  geom_vline(xintercept = 0, size = 1)+
  geom_hline(yintercept = 1.5, size = 1)+ geom_hline(yintercept = 2.5, size = 1)+ 
  geom_hline(yintercept = 3.5, size = 1)+ geom_hline(yintercept = 4.5, size = 1)+ geom_hline(yintercept = 5.5, size = 1)+
  geom_hline(yintercept = 6.5, size = 1)+ geom_hline(yintercept = 7.5, size = 1)+
  geom_pointrange(aes(mean_log2_mock_norm, timepoint,
                      xmin = mean_log2_mock_norm-sd_log2_mock_norm, xmax = mean_log2_mock_norm+sd_log2_mock_norm),
                  position = position_dodge(width = 0.66), size = 2, shape = 22)+
  geom_point(aes(log2_mock_norm, timepoint),
             position = position_dodge(width = 0.66), size = 5, fill = 'grey24', alpha = 0.5, shape = 22)+
  scale_fill_gradient2(name = 'log2(fold change)', low = "blue", mid = "white", high = "red", midpoint = 0, limits=c(-2,1.5))+
  scale_y_discrete(limits=rev)+
  theme_cu()+
  theme(legend.position = 'right',
        aspect.ratio = 1.5,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.minor.y = element_line(size = 2, colour = 'black'),
        legend.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        #axis.text.x = element_blank(),
        panel.border = element_rect(size = 2, colour = 'black'))+
  labs(x = '', y = '')+
  xlim(-3.5,2.5)

setwd("/Volumes/GoogleDrive/My Drive/qPCR/Big_exp")
cfu <- read_csv('Big_cfu.csv')

cfu.sub <- cfu %>% filter(timepoint %in% c(24,48,96)) %>% 
  group_by(treatment, timepoint) %>% mutate(log_cfu = log10(cfu_g)) %>%
  summarise(mean_log_cfu = mean(log_cfu), se_log_cfu = sd(log_cfu)/sqrt(n())) %>%
  mutate(treatment=ifelse(treatment=="MM1","Pst", 
                      ifelse(treatment=="MM27","Willi354", 
                             ifelse(treatment=="MM22","Micro347", 
                                    ifelse(treatment=="MM15","Sphingo34", "not_interested")))) )

marker.sub <- marker_gg %>% select(treatment, primer_pair, mean_log2_mock_norm, sd_log2_mock_norm) %>% distinct() %>%
  filter(primer_pair %in% 'et.marker') %>% filter(timepoint %in% c('24', '48', '96')) %>%
  mutate(timepoint=ifelse(timepoint=="24",24, 
                          ifelse(timepoint=="48",48,
                               ifelse(timepoint=="96",96,0))))
et.cor <- left_join(marker.sub, cfu.sub, by=c('timepoint', 'treatment')) %>% filter(treatment != 'Mock')

et.cor %>% group_by(timepoint,treatment) %>% mutate(se_log2_mock_norm = sd_log2_mock_norm/sqrt(3)) %>%
  ggplot()+
  stat_smooth(aes(mean_log_cfu, mean_log2_mock_norm),method = "lm", col = "red")+
  geom_errorbarh(aes(xmin=mean_log_cfu-se_log_cfu,xmax=mean_log_cfu+se_log_cfu,y=mean_log2_mock_norm),width=0.1)+
  geom_errorbar(aes(x=mean_log_cfu,
                      ymin=mean_log2_mock_norm-se_log2_mock_norm, ymax=mean_log2_mock_norm+se_log2_mock_norm),width=0.01)+
  geom_point(aes(mean_log_cfu, mean_log2_mock_norm,fill = as.character(timepoint)),
                  size = 5, shape=21)+
  scale_fill_manual(values = c('#abe3ce','#f4d984','#d6aaea'), name = 'hpi')+
  theme_cu()+
  theme(aspect.ratio = 1,
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.grid.minor.y = element_line(size = 2, colour = 'black'),
        #legend.title = element_blank(),
        axis.text = element_text(size = 20),
        panel.border = element_rect(size = 2, colour = 'black'))+
  labs(x = 'bacterial density (log(cfu/g))', y = 'log2(fold change)')

summary(lm(mean_log_cfu ~ mean_log2_mock_norm, et.cor))  

### stats
st <- function(pp,t){
z <- marker_gg %>% 
  filter(primer_pair %in% pp) %>%
  filter(treatment %in% t)
out <- TukeyHSD(aov(log2_mock_norm ~ timepoint, data = z))
return(out)
}
st2 <- function(pp,t){
  z <- marker_gg %>% 
    filter(primer_pair %in% pp) %>%
    filter(treatment %in% t)
  out <- summary(aov(log2_mock_norm ~ timepoint, data = z))
  return(out)
}

st('et.marker', 'Sphingo34') # ns
st('ja.marker', 'Sphingo34') # ns
st('sa.marker', 'Sphingo34') # ns

st('et.marker', 'Micro347') # 12 to everything else
st2('et.marker', 'Micro347')
st('ja.marker', 'Micro347') # ns
st('sa.marker', 'Micro347') # ns

st('et.marker', 'Willi354') #96-1;96-3;96-9;48-12;96-12
st2('et.marker', 'Willi354')
st('ja.marker', 'Willi354') # ns
st('sa.marker', 'Willi354') #48-09
st2('sa.marker', 'Willi354')

st('et.marker', 'Pst') #24-1;48-1;96-1; 48-3;96-3; 48-12; 96-12
st2('et.marker', 'Pst')
st('ja.marker', 'Pst') # ns
st('sa.marker', 'Pst') # ns



m2 <- marker_gg %>% filter(treatment %in% 'Mock') %>%
  group_by(primer_pair, timepoint) %>%
  summarise(log2_mm_norm = log2_mock_norm, biol_rep = biol_rep) %>%
  right_join(marker_gg,., by= c('timepoint','primer_pair')) %>% distinct()
## t.test
tt <- m2 %>% group_by(timepoint, treatment, primer_pair) %>% filter(treatment!= 'Pst') %>%
  summarise(tt = t.test(log2_mock_norm, log2_mm_norm, na.rm = TRUE)$p.value, ttbon = tt*96)

tt.pst <- m2 %>% group_by(timepoint, treatment, primer_pair) %>% filter(treatment%in% 'Pst') %>%
  filter(timepoint != '24') %>%
  summarise(tt = t.test(log2_mock_norm, log2_mm_norm, na.rm = TRUE)$p.value, ttbon = tt*96)
t.test(m2$log2_mock_norm, m2$log2_mm_norm, na.rm = TRUE)$p.value

## wilcox test
# https://www.r-tutor.com/elementary-statistics/non-parametric-methods/wilcoxon-signed-rank-test
wt <- m2 %>% group_by(timepoint, treatment, primer_pair) %>% filter(treatment!= 'Pst') %>%
  summarise(wt = wilcox.test(log2_mock_norm, log2_mm_norm, paired = TRUE)$p.value, wtbon = wt*96)

wt.pst <- m2 %>% group_by(timepoint, treatment, primer_pair) %>% filter(treatment%in% 'Pst') %>%
  filter(timepoint != '24') %>%
  summarise(wt = wilcox.test(log2_mock_norm, log2_mm_norm, paired = TRUE)$p.value, wtbon = wt*96)
wilcox.test(m2$log2_mock_norm, m2$log2_mm_norm, paired = TRUE)$p.value
          
