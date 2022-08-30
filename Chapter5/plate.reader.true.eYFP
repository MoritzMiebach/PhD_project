### Plate reader script for true eFYP expression

library(ggplot2)
library(qqplotr)
library(tidyverse)
library(reshape2)
library(readr)
library(MESS) #needed for 'auc' function
citation(package = "MESS")
library(wesanderson)
library(ggQQunif)
library(corrplot)
library(RColorBrewer)
detach(package:plyr)
library(dplyr)
library(gganimate)
library(av)
library(gifski)
library(ffmpeg)
library(transformr)
library(ragg)
library(dunn.test)
library(ComplexHeatmap)

# master_function that does most of the computation needed on a per run basis
AUC_baselined_data <- function(x){
  ## wide to long data conversion
  dat <- gather(x, well, fluorescence_intensity, A1:H12)
  
  ## split well into plate_column & plate_row
  dat_sample <- as.data.frame(str_split(dat$well,"", simplify = TRUE)) # 'simplify' returns a matrix rather than a list 
  names(dat_sample)[1] <- 'plate_row'
  dat_sample <- dat_sample %>% unite('plate_column', c('V2','V3'), sep = '', remove = TRUE)
  dat <- bind_cols(dat, dat_sample)
  
  ## add plant_type info
  dat <- mutate(dat, plant_type = ifelse(plate_row %in% c('A', 'B', 'C', 'D'), 'pPR1:nls-YFP',
                                         ifelse(plate_row %in% c('E', 'F', 'G', 'H'), 'Col-0', 'error')))
  
  ## add treatment info
  dat <- dat %>% arrange(Time, plate_row) %>% add_column(treatment = rep(c('mock', 'mock', '0.01125mM_SA', '0.01125mM_SA', '0.0225mM_SA', '0.0225mM_SA', '0.045mM_SA', '0.045mM_SA',
                                                                           '0.09mM_SA', '0.09mM_SA', '0.18mM_SA', '0.18mM_SA', 1:36), times = nrow(dat)/48))
  
  
  ## add info of origin
  dat <- dat %>% group_by(plant_type, Time) %>% add_column(epi_endo = rep(c(rep(c('Control'), times = 12), 
                                                                            rep(c('Epiphyte'), times = 18), rep(c('Endophyte'), times = 18)), times = nrow(dat)/48))
  
  ## baseline subtraction
  dat_baseline <- dat %>% group_by(well) %>% filter(Time < 1) %>% summarise(mean_fluo_int_1h = mean(fluorescence_intensity))
  dat <- dat %>% group_by(well) %>% right_join(dat_baseline, dat, by = 'well')
  dat <- dat %>% filter(Time >=1)
  dat <- dat %>% group_by(well) %>% mutate(bs_fluo_int = fluorescence_intensity - mean_fluo_int_1h)
  
  ## spline the data
  dat <- dat %>% group_by(well) %>%
    #filter(well %in% 'A1') %>% 
    mutate(spl = smooth.spline(x = Time, y = bs_fluo_int, nknots = 4)$y)
  
  ## subtract Col-0 spline data from pPR1:eYFP-nls spline data to get true eYFP expression curve
  dat.col <- dat %>% filter(plant_type=='Col-0') %>% group_by(Time,treatment) %>% summarise(spl_col=mean(spl))
  
  dat <- left_join(dat,dat.col, by=c('Time','treatment')) %>% mutate(spl_yfp=spl-spl_col)
  
  ## the area under curve function either subtracts areas in the negative space or adds them.
  ## however the area under the whole curve (from baseline and not from 0) needs to be considered.
  ## get lowest value of curve to be 0 by adding the min value of to all values for each sample.
  dat <- ungroup(dat)
  dat <- dat %>% group_by(plant_type) %>% mutate(pos_min = sqrt(min(spl_yfp, na.rm = TRUE)^2))
  dat <- dat %>% mutate(spl_min_corr = spl_yfp + pos_min)
  
  ## calculate AUC for each well
  dat <- dat %>% group_by(well) %>% mutate(area_under_curve = auc(x = Time, y = spl_min_corr, method = 'linear', absolutearea = TRUE, ties = min)) #'linear' means linear interpolation for the trapezoid method

  ## normalise AUCs against mock !Do this for each plant_type!
  mock_AUC_PR1 <- mean(dat  %>% filter(well %in% c('A1', 'A2')) %>% pull(var = area_under_curve))
  mock_AUC_Col <- mean(dat  %>% filter(well %in% c('E1', 'E2')) %>% pull(var = area_under_curve))
  
  dat_PR1 <- dat %>% filter(plant_type %in% 'pPR1:nls-YFP') %>% group_by(well) %>% mutate(norm_AUC = area_under_curve / mock_AUC_PR1)
  dat_Col <- dat %>% filter(plant_type %in% 'Col-0') %>% group_by(well) %>% mutate(norm_AUC = area_under_curve  / mock_AUC_Col)
  
  dat <- bind_rows(dat_PR1, dat_Col)
  
  return(dat_PR1)
}

# set working directory
setwd("~/My Drive/protoplast_Mo/transpose_clean_forR")

# data import individual 
X2020_12_07_plant1 <- read.csv("2020.12.07_plant1.csv")
X2020_12_08_plant2 <- read.csv("2020.12.08_plant2.csv")
X2020_12_09_plant3 <- read.csv("2020.12.09_plant3.csv")
X2020_12_10_plant4 <- read.csv("2020.12.10_plant4.csv")
X2020_12_11_plant5 <- read.csv("2020.12.11_plant5.csv")

X2020_12_14_plant1 <- read.csv("2020.12.14_plant1.csv")
X2020_12_15_plant2 <- read.csv("2020.12.15_plant2.csv")
X2020_12_16_plant3 <- read.csv("2020.12.16_plant3.csv")
X2020_12_17_plant4 <- read.csv("2020.12.17_plant4.csv")
X2020_12_18_plant5 <- read.csv("2020.12.18_plant5.csv")

X2021_01_18_plant1 <- read.csv("2021.01.18_plant1.csv")
X2021_01_19_plant2 <- read.csv("2021.01.19_plant2.csv")
X2021_01_20_plant3 <- read.csv("2021.01.20_plant3.csv")
X2021_01_21_plant4 <- read.csv("2021.01.21_plant4.csv")
X2021_01_22_plant5 <- read.csv("2021.01.22_plant5.csv")

X2021_01_25_plant1 <- read.csv("2021.01.25_plant1.csv")
X2021_01_26_plant2 <- read.csv("2021.01.26_plant2.csv")
X2021_01_27_plant3 <- read.csv("2021.01.27_plant3.csv")
X2021_01_28_plant4 <- read.csv("2021.01.28_plant4.csv")
X2021_01_29_plant5 <- read.csv("2021.01.29_plant5.csv")

## !!! I accidentally swapped AC235 and AC236 (C7 & C8) in 2020.12.17 --> correct this here
## extract values form each column and then insert them to the right column. --> this way the column order doesn't get fucked up
C7new <- X2020_12_17_plant4$C8
C8new <- X2020_12_17_plant4$C7
X2020_12_17_plant4$C7 <- C7new
X2020_12_17_plant4$C8 <- C8new

## limit data from plant 5 to first 24 h
X2020_12_11_plant5 <- X2020_12_11_plant5[1:288,]
X2020_12_18_plant5 <- X2020_12_18_plant5[1:288,]
X2021_01_22_plant5 <- X2021_01_22_plant5[1:288,]
X2021_01_29_plant5 <- X2021_01_29_plant5[1:288,] 

## perform the main analysis
plant1.1 <- AUC_baselined_data(X2020_12_07_plant1)
plant2.1 <- AUC_baselined_data(X2020_12_08_plant2)
plant3.1 <- AUC_baselined_data(X2020_12_09_plant3)
plant4.1 <- AUC_baselined_data(X2020_12_10_plant4)
plant5.1 <- AUC_baselined_data(X2020_12_11_plant5)

plant1.2 <- AUC_baselined_data(X2020_12_14_plant1)
plant2.2 <- AUC_baselined_data(X2020_12_15_plant2)
plant3.2 <- AUC_baselined_data(X2020_12_16_plant3)
plant4.2 <- AUC_baselined_data(X2020_12_17_plant4)
plant5.2 <- AUC_baselined_data(X2020_12_18_plant5)

plant1.3 <- AUC_baselined_data(X2021_01_18_plant1)
plant2.3 <- AUC_baselined_data(X2021_01_19_plant2)
plant3.3 <- AUC_baselined_data(X2021_01_20_plant3)
plant4.3 <- AUC_baselined_data(X2021_01_21_plant4)
plant5.3 <- AUC_baselined_data(X2021_01_22_plant5)

plant1.4 <- AUC_baselined_data(X2021_01_25_plant1)
plant2.4 <- AUC_baselined_data(X2021_01_26_plant2)
plant3.4 <- AUC_baselined_data(X2021_01_27_plant3)
plant4.4 <- AUC_baselined_data(X2021_01_28_plant4)
plant5.4 <- AUC_baselined_data(X2021_01_29_plant5)

## add strain identity + plant ID
strain_IDs <- read.csv("protoplast_assay_strain_ids_2.csv")
strain_ID_plant1 <- strain_IDs %>% dplyr::select(plate_row, plate_column, id_plant1)
strain_ID_plant1$plate_column <- as.character(strain_ID_plant1$plate_column)
colnames(strain_ID_plant1)[3] <- 'strain_ID'
strain_ID_plant1$plant_ID <- rep('plant1', 96)
strain_ID_plant2 <- strain_IDs %>% dplyr::select(plate_row, plate_column, id_plant2)
strain_ID_plant2$plate_column <- as.character(strain_ID_plant2$plate_column)
colnames(strain_ID_plant2)[3] <- 'strain_ID'
strain_ID_plant2$plant_ID <- rep('plant2', 96)
strain_ID_plant3 <- strain_IDs %>% dplyr::select(plate_row, plate_column, id_plant3)
strain_ID_plant3$plate_column <- as.character(strain_ID_plant3$plate_column)
colnames(strain_ID_plant3)[3] <- 'strain_ID'
strain_ID_plant3$plant_ID <- rep('plant3', 96)
strain_ID_plant4 <- strain_IDs %>% dplyr::select(plate_row, plate_column, id_plant4)
strain_ID_plant4$plate_column <- as.character(strain_ID_plant4$plate_column)
colnames(strain_ID_plant4)[3] <- 'strain_ID'
strain_ID_plant4$plant_ID <- rep('plant4', 96)
strain_ID_plant5 <- strain_IDs %>% dplyr::select(plate_row, plate_column, id_plant5)
strain_ID_plant5$plate_column <- as.character(strain_ID_plant5$plate_column)
colnames(strain_ID_plant5)[3] <- 'strain_ID'
strain_ID_plant5$plant_ID <- rep('plant5', 96)

plant1.1 <- left_join(plant1.1, strain_ID_plant1, by = c("plate_row","plate_column"))
plant1.2 <- left_join(plant1.2, strain_ID_plant1, by = c("plate_row","plate_column"))
plant1.3 <- left_join(plant1.3, strain_ID_plant1, by = c("plate_row","plate_column"))
plant1.4 <- left_join(plant1.4, strain_ID_plant1, by = c("plate_row","plate_column"))

plant2.1 <- left_join(plant2.1, strain_ID_plant2, by = c("plate_row","plate_column"))
plant2.2 <- left_join(plant2.2, strain_ID_plant2, by = c("plate_row","plate_column"))
plant2.3 <- left_join(plant2.3, strain_ID_plant2, by = c("plate_row","plate_column"))
plant2.4 <- left_join(plant2.4, strain_ID_plant2, by = c("plate_row","plate_column"))

plant3.1 <- left_join(plant3.1, strain_ID_plant3, by = c("plate_row","plate_column"))
plant3.2 <- left_join(plant3.2, strain_ID_plant3, by = c("plate_row","plate_column"))
plant3.3 <- left_join(plant3.3, strain_ID_plant3, by = c("plate_row","plate_column"))
plant3.4 <- left_join(plant3.4, strain_ID_plant3, by = c("plate_row","plate_column"))

plant4.1 <- left_join(plant4.1, strain_ID_plant4, by = c("plate_row","plate_column"))
plant4.2 <- left_join(plant4.2, strain_ID_plant4, by = c("plate_row","plate_column"))
plant4.3 <- left_join(plant4.3, strain_ID_plant4, by = c("plate_row","plate_column"))
plant4.4 <- left_join(plant4.4, strain_ID_plant4, by = c("plate_row","plate_column"))

plant5.1 <- left_join(plant5.1, strain_ID_plant5, by = c("plate_row","plate_column"))
plant5.2 <- left_join(plant5.2, strain_ID_plant5, by = c("plate_row","plate_column"))
plant5.3 <- left_join(plant5.3, strain_ID_plant5, by = c("plate_row","plate_column"))
plant5.4 <- left_join(plant5.4, strain_ID_plant5, by = c("plate_row","plate_column"))

## add replicate ID
replicate_ID_1 <- strain_ID_plant1[,1:2]
replicate_ID_1$replicate_ID <- rep('replicate1', 96)
replicate_ID_2 <- strain_ID_plant1[,1:2]
replicate_ID_2$replicate_ID <- rep('replicate2', 96)
replicate_ID_3 <- strain_ID_plant1[,1:2]
replicate_ID_3$replicate_ID <- rep('replicate3', 96)
replicate_ID_4 <- strain_ID_plant1[,1:2]
replicate_ID_4$replicate_ID <- rep('replicate4', 96)
replicate_ID_5 <- strain_ID_plant1[,1:2]
replicate_ID_5$replicate_ID <- rep('replicate5', 96)

plant1.1 <- left_join(plant1.1, replicate_ID_1, by = c("plate_row","plate_column"))
plant2.1 <- left_join(plant2.1, replicate_ID_1, by = c("plate_row","plate_column"))
plant3.1 <- left_join(plant3.1, replicate_ID_1, by = c("plate_row","plate_column"))
plant4.1 <- left_join(plant4.1, replicate_ID_1, by = c("plate_row","plate_column"))
plant5.1 <- left_join(plant5.1, replicate_ID_1, by = c("plate_row","plate_column"))

plant1.2 <- left_join(plant1.2, replicate_ID_2, by = c("plate_row","plate_column"))
plant2.2 <- left_join(plant2.2, replicate_ID_2, by = c("plate_row","plate_column"))
plant3.2 <- left_join(plant3.2, replicate_ID_2, by = c("plate_row","plate_column"))
plant4.2 <- left_join(plant4.2, replicate_ID_2, by = c("plate_row","plate_column"))
plant5.2 <- left_join(plant5.2, replicate_ID_2, by = c("plate_row","plate_column"))

plant1.3 <- left_join(plant1.3, replicate_ID_3, by = c("plate_row","plate_column"))
plant2.3 <- left_join(plant2.3, replicate_ID_3, by = c("plate_row","plate_column"))
plant3.3 <- left_join(plant3.3, replicate_ID_3, by = c("plate_row","plate_column"))
plant4.3 <- left_join(plant4.3, replicate_ID_3, by = c("plate_row","plate_column"))
plant5.3 <- left_join(plant5.3, replicate_ID_3, by = c("plate_row","plate_column"))

plant1.4 <- left_join(plant1.4, replicate_ID_4, by = c("plate_row","plate_column"))
plant2.4 <- left_join(plant2.4, replicate_ID_4, by = c("plate_row","plate_column"))
plant3.4 <- left_join(plant3.4, replicate_ID_4, by = c("plate_row","plate_column"))
plant4.4 <- left_join(plant4.4, replicate_ID_4, by = c("plate_row","plate_column"))
plant5.4 <- left_join(plant5.4, replicate_ID_4, by = c("plate_row","plate_column"))

###############################################################################################################################
###############################################################################################################################  
## custom ggplot theme
theme_cu <- function() {
  theme_minimal() +
    theme(panel.grid.major = element_line(colour = "grey69", linetype = 'solid', size = 0.5),
          axis.title = element_text(colour = "black", family = "sans", size = 16),
          # legend at the bottom 
          legend.position = 'bottom',
          legend.justification = "center",
          legend.title = element_text(size =0),
          legend.text = element_text(colour = "black", family = "sans", size = 16),
          # plot title
          plot.title = element_text(colour = "black", family = "sans", size = 20),
          aspect.ratio = 1,
          panel.grid.minor = element_blank(),
          axis.text = element_text(colour = "black", family = "sans", size = 12),
          panel.border = element_rect(size = 2, colour = 'black', fill = NA)
    )
}

yfp_plot <- function(x,y){
  dat <- x
  ## kick out x treatments as they are inconsistent
  dat <- dat %>% filter(strain_ID != 'x')
  ## calculate threshold as median + sd of all AUC values per run (excluding the kicked out SA treatments)
  dat <- ungroup(dat)
  dat <- dat %>% filter(plant_type %in%'pPR1:nls-YFP') %>% mutate(run_median_AUC = median(norm_AUC))
  dat <- dat %>% filter(plant_type %in%'pPR1:nls-YFP') %>% mutate(run_sd_AUC = sd(norm_AUC))
  
  yfp2 <- dat %>% filter(plant_type %in%'pPR1:nls-YFP') %>% 
    ggplot(aes(strain_ID, spl_min_corr, fill=Time,col=Time))+
    geom_point(alpha=.3,size=2,shape=21)+
    scale_colour_viridis_c(name='time (h)')+
    scale_fill_viridis_c(name='time (h)')+
    theme_cu() %+replace% theme(axis.text.x = element_text(angle=0, hjust=1))+
    theme(axis.ticks.x=element_blank())+
    theme(aspect.ratio = 2,
          legend.position = 'none',
          legend.title = element_text(size = 20),
          legend.title.align = 1,
          panel.grid.major = element_blank())+
    labs(title = paste0('eYFP | ',y), x = '', y = '')+
    coord_flip()+
    ylim(0,26000)
  
  yfp3 <- dat %>% filter(plant_type %in%'pPR1:nls-YFP') %>% 
    ggplot(aes(Time, spl_min_corr, col=strain_ID))+
    geom_line(lwd=1)+
    #scale_colour_viridis_d(name='time (h)')+
    scale_fill_viridis_c(name='time (h)')+
    theme_cu() %+replace% theme(axis.text.x = element_text(angle=0, hjust=1))+
    theme(axis.ticks.x=element_blank())+
    theme(aspect.ratio = 1.5,
          legend.position = 'right',
          legend.title = element_text(size = 20),
          legend.title.align = 1,
          panel.grid.major = element_blank())+
    labs(title = paste0('eYFP | ',y), x = '', y = '')+
    ylim(0,26000)
  
  return(yfp2)
}
setwd("~/My Drive/protoplast_Mo/protoplastJuly2022/range.plots.yfp")

max(plant5.4$spl_min_corr)

yfp_plot(plant1.1,'p1.1')
ggsave('y11.pdf',width=5,height=7)
yfp_plot(plant1.2,'p1.2')
ggsave('y12.pdf',width=5,height=7)
yfp_plot(plant1.3,'p1.3')
ggsave('y13.pdf',width=5,height=7)
yfp_plot(plant1.4,'p1.4')
ggsave('y14.pdf',width=5,height=7)

yfp_plot(plant2.1,'p2.1')
ggsave('y21.pdf',width=5,height=7)
yfp_plot(plant2.2,'p2.2')
ggsave('y22.pdf',width=5,height=7)
yfp_plot(plant2.3,'p2.3')
ggsave('y23.pdf',width=5,height=7)
yfp_plot(plant2.4,'p2.4')
ggsave('y24.pdf',width=5,height=7)

yfp_plot(plant3.1,'p3.1')
ggsave('y31.pdf',width=5,height=7)
yfp_plot(plant3.2,'p3.2')
ggsave('y32.pdf',width=5,height=7)
yfp_plot(plant3.3,'p3.3')
ggsave('y33.pdf',width=5,height=7)
yfp_plot(plant3.4,'p3.4')
ggsave('y34.pdf',width=5,height=7)

yfp_plot(plant4.1,'p4.1')
ggsave('y41.pdf',width=5,height=7)
yfp_plot(plant4.2,'p4.2')
ggsave('y42.pdf',width=5,height=7)
yfp_plot(plant4.3,'p4.3')
ggsave('y43.pdf',width=5,height=7)
yfp_plot(plant4.4,'p4.4')
ggsave('y44.pdf',width=5,height=7)

yfp_plot(plant5.1,'p5.1')
ggsave('y51.pdf',width=5,height=7)
yfp_plot(plant5.2,'p5.2')
ggsave('y52.pdf',width=5,height=7)
yfp_plot(plant5.3,'p5.3')
ggsave('y53.pdf',width=5,height=7)
yfp_plot(plant5.4,'p5.4')
ggsave('y54.pdf',width=5,height=7)

## joining the data -> heavily nested since 'union' only allows for two data to get combined
all_prot <- union(union(union(union(plant1.1, plant1.2), union(plant1.3, plant1.4)), union(union(plant2.1, plant2.2), union(plant2.3, plant2.4))),
                  union(union(union(plant3.1, plant3.2), union(plant3.3, plant3.4)), union(union(plant4.1, plant4.2), union(plant4.3, plant4.4))))
all_prot <- union(all_prot, union(union(plant5.1, plant5.2), union(plant5.3, plant5.4)))

## get rid of fluorescent bacteria (49,132,145,269,886) 
## & potentially fluorescent bacteria (302,463,248,68,4,884)
ap <- all_prot %>% filter(strain_ID != 'x') %>% filter(strain_ID != '49') %>% 
  filter(strain_ID != '132') %>% filter(strain_ID != '145') %>% 
  filter(strain_ID != '269') %>% filter(strain_ID != '886')

## reduce to AUC data
apr <- ap %>% 
  #filter(replicate_ID %in% c('replicate1','replicate2')) %>%
  dplyr::select(well,strain_ID,epi_endo,plant_ID,norm_AUC,area_under_curve,replicate_ID) %>% 
  distinct()
## transfrom data to log2
# AUC raw
#apr <- apr %>% group_by(strain_ID, plant_ID) %>%
#  mutate(log.auc = log2(area_under_curve),m.log.auc = mean(log.auc),sd.log.auc = sd(log.auc)) 
# AUC normalised to mock
mock.auc <- apr %>% group_by(replicate_ID,plant_ID) %>% filter(epi_endo=='Control') %>%
  summarise(mock.auc=mean(area_under_curve))
apr <- left_join(apr,mock.auc,by=c('replicate_ID','plant_ID'))
apr <- apr %>% group_by(strain_ID, plant_ID) %>%
  mutate(n.auc=area_under_curve/mock.auc,log.auc = log2(n.auc),
         m.log.auc = mean(log.auc),sd.log.auc = sd(log.auc)) 
#apr <- apr %>% group_by(strain_ID, plant_ID) %>%
#  mutate(log.auc = log2(norm_AUC),m.log.auc = mean(log.auc),sd.log.auc = sd(log.auc)) 
# AUC normalised to run median
#med.auc <- apr %>% group_by(replicate_ID,plant_ID) %>% summarise(med.auc=median(area_under_curve))
#apr <- left_join(apr,med.auc,by=c('replicate_ID','plant_ID'))
#apr <- apr %>% group_by(strain_ID, plant_ID) %>%
#  mutate(n.auc=area_under_curve/med.auc,log.auc = log2(n.auc),
#        m.log.auc = mean(log.auc),sd.log.auc = sd(log.auc)) 

## get a dataset with only mean of all replicates and map qqplot x values
apr_qq <- apr %>% dplyr::select(strain_ID, plant_ID,epi_endo,m.log.auc) %>% distinct()
apr_qq_epi <- apr_qq %>% filter(epi_endo %in% 'Epiphyte')
apr_qq_endo <- apr_qq %>% filter(epi_endo %in% 'Endophyte')
apr_qq_control <- apr_qq %>% filter(epi_endo %in% 'Control')

apr_qq_epi <- apr_qq_epi %>% arrange(m.log.auc) %>% add_column(qq = qnorm(seq(0.01, 0.99, 1/87)))
apr_qq_endo <- apr_qq_endo %>% arrange(m.log.auc) %>% add_column(qq = qnorm(seq(0.01, 0.99, 1/90)))
apr_qq_control <- apr_qq_control %>% arrange(m.log.auc) %>% add_column(qq = c(0, 0, 0, 0, 0))

apr_qq <- union(union(apr_qq_epi, apr_qq_endo), apr_qq_control)
apr_qq$strain_ID <- with(apr_qq, reorder(strain_ID, m.log.auc))
apr_qq <- apr_qq %>% arrange(m.log.auc) %>% add_column(qq_comb = qnorm(seq(0.01, 0.99, 1/183)))
apr_qq <- apr_qq %>% ungroup() %>% 
  summarise(strain_ID = strain_ID, qq = qq, epi_endo = epi_endo, qq_comb = qq_comb) %>%
  dplyr::select(-c('epi_endo'))


ap <- left_join(ap, apr_qq, by = "strain_ID")
ap %>% filter(replicate_ID == 'replicate1') %>%
  group_by(strain_ID,plant_ID,Time) %>% summarise(spl_min_corr=mean(spl_min_corr),qq_comb) %>%
  ggplot(aes(qq_comb, spl_min_corr, fill=Time,col=Time))+
  geom_point(alpha=.3,size=2,shape=21)+
  scale_colour_viridis_c(name='time (h)')+
  scale_fill_viridis_c(name='time (h)')+
  theme_cu() %+replace% theme(axis.text.x = element_text(angle=0, hjust=1))+
  theme(axis.ticks.x=element_blank())+
  theme(aspect.ratio = 2,
        legend.position = 'none',
        legend.title = element_text(size = 20),
        legend.title.align = 1,
        panel.grid.major = element_blank())+
  labs(x = '', y = '')+
  coord_flip()+
  ylim(0,26000)


apr <- left_join(apr, apr_qq, by = "strain_ID")

scaleFUN <- function(x) sprintf("%.1f", x)

apr %>% filter(epi_endo %in% 'Control') %>%
  group_by(strain_ID) %>% 
  mutate(se = sd(log.auc)/sqrt(4)) %>%
  ggplot()+
  #geom_hline(yintercept = 0, size = 1, color = "grey30")+
  #geom_pointrange(aes(qq_comb, mean_log2_AUC, ymin = mean_log2_AUC-se_log2_AUC, ymax = mean_log2_AUC+se_log2_AUC, col = epi_endo))+
  geom_point(aes(plant_ID, log.auc, fill = replicate_ID), size = 4, shape = 21,
             alpha=.7, position=position_dodge(width=.15))+
  #coord_flip()+
  scale_fill_manual(values = c('#4053d3','#ddb310','#b51d14','#00beff'))+
  scale_color_manual(values = c('#4053d3','#ddb310','#b51d14','#00beff'))+
  labs(y = 'log2(normalised AUC)', x = 'run order')+
  scale_y_continuous(labels=scaleFUN)+
  theme_cu()+
  theme(aspect.ratio = .3,
        axis.ticks.x=element_blank())

apr %>% group_by(strain_ID) %>% 
  mutate(se = sd(log.auc)/sqrt(4)) %>%
  ggplot()+
  #geom_hline(yintercept = 0, size = 1, color = "grey30")+
  geom_hline(yintercept = 0, size = 1, color = "grey24")+
  #geom_smooth(aes(qq_comb, log.auc),method='lm', lwd=2, fill='blue')+
  geom_point(aes(qq_comb, log.auc, fill = replicate_ID), size = 2, shape = 21, alpha=.7)+
  #coord_flip()+
  scale_fill_manual(values = c('#4053d3','#ddb310','#b51d14','#00beff'))+
  scale_y_continuous(labels=scaleFUN)+
  labs(y = 'log2(normalisedAUC)', x = 'theoretical quantiles')+
  theme_cu()+
  theme(aspect.ratio = .3,
        axis.ticks.x=element_blank())

apr %>% group_by(strain_ID) %>% mutate(se = sd(log.auc)/sqrt(4)) %>%
  filter(epi_endo != 'Control') %>%
  #filter(sum_fluo == 4) %>%
  ggplot()+
  #geom_hline(yintercept = 0, size = 1, color = "grey30")+
  geom_hline(yintercept = 0, size = 1, color = 'grey24')+
  #geom_hline(yintercept = 0.343, size = 1, color = '#78B7C5')+
  #geom_hline(yintercept = 0.549, size = 1, color = '#3B9AB2')+
  #geom_rect(data = data.frame(xmin=-Inf, xmax=Inf, ymin=-0.0294, ymax=0.0236),
  #          aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='#f21a00', alpha=0.2, inherit.aes = FALSE)+
  #geom_rect(data = data.frame(xmin=-Inf, xmax=Inf, ymin=0.327, ymax=0.359),
  #          aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='#78B7C5', alpha=0.2, inherit.aes = FALSE)+
  #geom_rect(data = data.frame(xmin=-Inf, xmax=Inf, ymin=0.531, ymax=0.567),
  #          aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='#3B9AB2', alpha=0.2, inherit.aes = FALSE)+
  #geom_pointrange(aes(qq, mean_log2_AUC, ymin = mean_log2_AUC-se_log2_AUC, ymax = mean_log2_AUC+se_log2_AUC,
  #                    col = epi_endo, shape = sum_fluo), size = 0.5)+
  geom_point(aes(qq, m.log.auc, fill = epi_endo), size = 2.5, shape = 21, alpha =.8)+
  geom_smooth(aes(qq, log.auc, col = epi_endo, fill = epi_endo), method = 'lm')+
  #coord_flip()+
  scale_color_manual(values = c('#74a089', '#9a8822'))+
  scale_fill_manual(values = c('#74a089', '#9a8822'))+
  scale_y_continuous(labels=scaleFUN)+
  scale_shape_manual(values = c(25, 24, 23, 21))+
  labs(y = 'log2(normalised AUC)', x = 'theoretical quantiles')+
  theme_cu()+
  theme(aspect.ratio=.3)
  ylim(-1, 0.5)+
  xlim(-2.2,2.2)

apr %>% 
  group_by(strain_ID) %>% mutate(se = sd(log.auc)/sqrt(4)) %>%
  #filter(epi_endo != 'Control') %>%
  ggplot()+
  #geom_hline(yintercept = 0, size = 1, color = "grey30")+
  #geom_hline(yintercept = 0, size = 1, color = 'grey24')+
  #geom_hline(yintercept = 0.343, size = 1, color = '#78B7C5')+
  #geom_hline(yintercept = 0.549, size = 1, color = '#3B9AB2')+
  #geom_rect(data = data.frame(xmin=-Inf, xmax=Inf, ymin=-0.0294, ymax=0.0236),
  #          aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='#f21a00', alpha=0.2, inherit.aes = FALSE)+
  #geom_rect(data = data.frame(xmin=-Inf, xmax=Inf, ymin=0.327, ymax=0.359),
  #          aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='#78B7C5', alpha=0.2, inherit.aes = FALSE)+
  #geom_rect(data = data.frame(xmin=-Inf, xmax=Inf, ymin=0.531, ymax=0.567),
  #          aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='#3B9AB2', alpha=0.2, inherit.aes = FALSE)+
  geom_pointrange(aes(qq_comb, m.log.auc, ymin = m.log.auc-se, ymax = m.log.auc+se,
                      col = epi_endo), size = 0.5)+
  geom_point(aes(qq_comb, m.log.auc, fill = epi_endo), size = 2.5, shape=21)+
  #geom_smooth(aes(qq_comb, log2_AUC, col = epi_endo, fill = epi_endo), method = 'lm')+
  #coord_flip()+
  scale_color_manual(values = c('#f21a00','#74a089', '#9a8822'))+
  scale_fill_manual(values = c('#f21a00','#74a089', '#9a8822'))+
  scale_y_continuous(labels=scaleFUN)+
  scale_shape_manual(values = c(22, 24, 21))+
  labs(y = 'log2(normalised AUC)', x = 'theoretical quantiles')+
  theme_cu()+
  theme(aspect.ratio=.3)
  #ylim(-0.5, 0.75)+
  xlim(-2.5,2.5)

apr %>% group_by(strain_ID) %>% mutate(se = sd(log.auc)/sqrt(4)) %>%
  #filter(epi_endo != 'Control') %>%
  ggplot()+
  #geom_hline(yintercept = 0, size = 1, color = "grey30")+
  #geom_hline(yintercept = 0, size = 1, color = 'grey24')+
  #geom_hline(yintercept = 0.343, size = 1, color = '#78B7C5')+
  #geom_hline(yintercept = 0.549, size = 1, color = '#3B9AB2')+
  #geom_rect(data = data.frame(xmin=-Inf, xmax=Inf, ymin=-0.0294, ymax=0.0236),
  #          aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='#f21a00', alpha=0.2, inherit.aes = FALSE)+
  #geom_rect(data = data.frame(xmin=-Inf, xmax=Inf, ymin=0.327, ymax=0.359),
  #          aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='#78B7C5', alpha=0.2, inherit.aes = FALSE)+
  #geom_rect(data = data.frame(xmin=-Inf, xmax=Inf, ymin=0.531, ymax=0.567),
  #          aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='#3B9AB2', alpha=0.2, inherit.aes = FALSE)+
  geom_pointrange(aes(qq_comb, m.log.auc, ymin = m.log.auc-se, ymax = m.log.auc+se,
                      col = plant_ID), size = 0.5, shape = 21)+
  geom_point(aes(qq_comb, m.log.auc, fill = plant_ID), size = 2.5, shape = 21)+
  #geom_smooth(aes(qq_comb, log2_AUC, col = epi_endo, fill = epi_endo), method = 'lm')+
  #coord_flip()+
  scale_color_manual(values = c( '#c0affb','#e6a176','#00678a','#984464','#5eccab'))+
  scale_fill_manual(values = c('#c0affb','#e6a176','#00678a','#984464','#5eccab'))+
  scale_y_continuous(labels=scaleFUN)+
  scale_shape_manual(values = c(22, 24, 23, 21))+
  labs(y = 'log2(normalised AUC)', x = 'theoretical quantiles')+
  theme_cu()+
  theme(aspect.ratio=.3)
  #ylim(-0.5, 0.75)+
  xlim(-2.5,2.5)
  
################################################################################
## strains on extreme ends of spectrum were selected for further investigation
  
setwd('~/My Drive/protoplast_Mo/transpose_clean_forR') 
cand29_strain_ID <- read.csv("strainID_new29.csv")
cand29_strain_ID$strain_ID <- as.character(cand29_strain_ID$strain_ID)
subset_cand29 <- semi_join(apr, cand29_strain_ID, by = 'strain_ID')

subset_cand29 %>% filter(epi_endo %in% c('Epiphyte', 'Endophyte')) %>%
  group_by(strain_ID) %>% mutate(se = sd(log.auc)/sqrt(4)) %>%
  ggplot()+
  #geom_hline(yintercept = 0, size = 1, color = "grey30")+
  #geom_hline(yintercept = 0, size = 1, color = 'grey24')+
  #geom_hline(yintercept = 0.343, size = 1, color = '#78B7C5')+
  #geom_hline(yintercept = 0.549, size = 1, color = '#3B9AB2')+
  #geom_rect(data = data.frame(xmin=-Inf, xmax=Inf, ymin=-0.0294, ymax=0.0236),
  #          aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='#f21a00', alpha=0.2, inherit.aes = FALSE)+
  #geom_rect(data = data.frame(xmin=-Inf, xmax=Inf, ymin=0.327, ymax=0.359),
  #          aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='#78B7C5', alpha=0.2, inherit.aes = FALSE)+
  #geom_rect(data = data.frame(xmin=-Inf, xmax=Inf, ymin=0.531, ymax=0.567),
  #          aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='#3B9AB2', alpha=0.2, inherit.aes = FALSE)+
  #geom_pointrange(aes(qq_comb, m.log.auc, ymin = m.log.auc-se, ymax = m.log.auc+se, col = epi_endo))+
  geom_point(aes(qq, m.log.auc, fill=epi_endo), size = 3.5, shape = 21)+
  #geom_smooth(aes(qq, log2_AUC, fill = epi_endo, col = epi_endo), method = 'lm')+
  #coord_flip()+
  #scale_color_manual(values = c('#74a089', '#9a8822'))+
  scale_fill_manual(values = c('#74a089', '#9a8822'))+
  scale_y_continuous(labels=scaleFUN)+
  scale_shape_manual(values = c(24, 23, 21))+
  labs(title = '', y = 'log2(normalised AUC)', x = 'theoretical quantiles')+
  theme_cu()+
  theme(aspect.ratio = .3)+
  ylim(-1, 0.5)+
  xlim(-2.2,2.2)
  
subset_cand29 %>% filter(epi_endo %in% c('Epiphyte', 'Endophyte')) %>%
    group_by(strain_ID) %>% mutate(se = sd(log.auc)/sqrt(4)) %>%
    ggplot()+
    #geom_hline(yintercept = 0, size = 1, color = "grey30")+
    #geom_hline(yintercept = 0, size = 1, color = 'grey24')+
    #geom_hline(yintercept = 0.343, size = 1, color = '#78B7C5')+
    #geom_hline(yintercept = 0.549, size = 1, color = '#3B9AB2')+
    #geom_rect(data = data.frame(xmin=-Inf, xmax=Inf, ymin=-0.0294, ymax=0.0236),
    #          aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='#f21a00', alpha=0.2, inherit.aes = FALSE)+
    #geom_rect(data = data.frame(xmin=-Inf, xmax=Inf, ymin=0.327, ymax=0.359),
    #          aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='#78B7C5', alpha=0.2, inherit.aes = FALSE)+
    #geom_rect(data = data.frame(xmin=-Inf, xmax=Inf, ymin=0.531, ymax=0.567),
    #          aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='#3B9AB2', alpha=0.2, inherit.aes = FALSE)+
    geom_pointrange(aes(qq_comb, m.log.auc, ymin = m.log.auc-se, ymax = m.log.auc+se, col = epi_endo))+
    geom_point(aes(qq_comb, m.log.auc, fill = epi_endo), size = 2.5, shape = 21)+
    #geom_smooth(aes(qq, log2_AUC, fill = epi_endo, col = epi_endo), method = 'lm')+
    #coord_flip()+
    scale_color_manual(values = c('#74a089', '#9a8822'))+
    scale_fill_manual(values = c('#74a089', '#9a8822'))+
    scale_y_continuous(labels=scaleFUN)+
    scale_shape_manual(values = c(24, 23, 21))+
    labs(title = '', y = 'log2(normalised AUC)', x = 'theoretical quantiles')+
    theme_cu()+
    theme(aspect.ratio = .3)+
    ylim(-1, 0.5)+
    xlim(-2.2,2.2)

apr %>% group_by(strain_ID) %>% mutate(se = sd(log.auc)/sqrt(4)) %>%
  filter(epi_endo != 'Control') %>%
  #filter(sum_fluo == 4) %>%
  ggplot()+
  #geom_hline(yintercept = 0, size = 1, color = "grey30")+
  geom_hline(yintercept = 0, size = 1, color = 'grey24')+
  #geom_hline(yintercept = 0.343, size = 1, color = '#78B7C5')+
  #geom_hline(yintercept = 0.549, size = 1, color = '#3B9AB2')+
  #geom_rect(data = data.frame(xmin=-Inf, xmax=Inf, ymin=-0.0294, ymax=0.0236),
  #          aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='#f21a00', alpha=0.2, inherit.aes = FALSE)+
  #geom_rect(data = data.frame(xmin=-Inf, xmax=Inf, ymin=0.327, ymax=0.359),
  #          aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='#78B7C5', alpha=0.2, inherit.aes = FALSE)+
  #geom_rect(data = data.frame(xmin=-Inf, xmax=Inf, ymin=0.531, ymax=0.567),
  #          aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='#3B9AB2', alpha=0.2, inherit.aes = FALSE)+
  #geom_pointrange(aes(qq, mean_log2_AUC, ymin = mean_log2_AUC-se_log2_AUC, ymax = mean_log2_AUC+se_log2_AUC,
  #                    col = epi_endo, shape = sum_fluo), size = 0.5)+
  geom_point(aes(qq, m.log.auc, fill = epi_endo), size = 2.5, shape = 21, alpha =.8)+
  geom_smooth(aes(qq, log.auc, col = epi_endo, fill = epi_endo), method = 'lm')+
  geom_point(data=subset_cand29,aes(qq, m.log.auc), size = 2.6, shape = 8, col='#DA2B27')+
  #coord_flip()+
  scale_color_manual(values = c('#74a089', '#9a8822'))+
  scale_fill_manual(values = c('#74a089', '#9a8822'))+
  scale_y_continuous(labels=scaleFUN)+
  scale_shape_manual(values = c(25, 24, 23, 21))+
  labs(y = 'log2(normalised AUC)', x = 'theoretical quantiles')+
  theme_cu()+
  theme(aspect.ratio=.3)
  #ylim(-1, 0.5)+
  #xlim(-2.2,2.2)

subset_cand29 %>%
  group_by(strain_ID) %>% mutate(se = sd(log.auc)/sqrt(4)) %>%
  ggplot()+
  #geom_hline(yintercept = 0, size = 1, color = "grey30")+
  #geom_hline(yintercept = 0, size = 1, color = 'grey24')+
  #geom_hline(yintercept = 0.343, size = 1, color = '#78B7C5')+
  #geom_hline(yintercept = 0.549, size = 1, color = '#3B9AB2')+
  #geom_rect(data = data.frame(xmin=-Inf, xmax=Inf, ymin=-0.0294, ymax=0.0236),
  #          aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='#f21a00', alpha=0.2, inherit.aes = FALSE)+
  #geom_rect(data = data.frame(xmin=-Inf, xmax=Inf, ymin=0.327, ymax=0.359),
  #          aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='#78B7C5', alpha=0.2, inherit.aes = FALSE)+
  #geom_rect(data = data.frame(xmin=-Inf, xmax=Inf, ymin=0.531, ymax=0.567),
  #          aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='#3B9AB2', alpha=0.2, inherit.aes = FALSE)+
  geom_pointrange(aes(qq_comb, m.log.auc, ymin = m.log.auc-se, ymax = m.log.auc+se, col = plant_ID))+
  geom_point(aes(qq_comb, m.log.auc, fill = plant_ID), size = 2.5, shape = 21)+
  #geom_smooth(aes(qq, log2_AUC, fill = epi_endo, col = epi_endo), method = 'lm')+
  #coord_flip()+
  scale_color_manual(values = c( '#c0affb','#e6a176','#00678a','#984464','#5eccab'))+
  scale_fill_manual(values = c('#c0affb','#e6a176','#00678a','#984464','#5eccab'))+
  scale_y_continuous(labels=scaleFUN)+
  scale_shape_manual(values = c(24, 23, 21))+
  labs(title = '', y = 'log2(normalised AUC)', x = 'theoretical quantiles')+
  theme_cu()+
  theme(aspect.ratio = .3)+
  ylim(-1, 0.5)+
  xlim(-2.2,2.2)
 
setwd('~/My Drive/protoplast_Mo/transpose_clean_forR') 
rna18_strain_ID <- read.csv("strainID_18RNA.csv")
rna18_strain_ID$strain_ID <- as.character(rna18_strain_ID$strain_ID)
subset_rna18 <- semi_join(apr, rna18_strain_ID, by = 'strain_ID')

subset_rna18 %>% filter(epi_endo %in% c('Epiphyte', 'Endophyte')) %>%
  group_by(strain_ID) %>% mutate(se = sd(log.auc)/sqrt(4)) %>%
  ggplot()+
  #geom_hline(yintercept = 0, size = 1, color = "grey30")+
  #geom_hline(yintercept = 0, size = 1, color = 'grey24')+
  #geom_hline(yintercept = 0.343, size = 1, color = '#78B7C5')+
  #geom_hline(yintercept = 0.549, size = 1, color = '#3B9AB2')+
  #geom_rect(data = data.frame(xmin=-Inf, xmax=Inf, ymin=-0.0294, ymax=0.0236),
  #          aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='#f21a00', alpha=0.2, inherit.aes = FALSE)+
  #geom_rect(data = data.frame(xmin=-Inf, xmax=Inf, ymin=0.327, ymax=0.359),
  #          aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='#78B7C5', alpha=0.2, inherit.aes = FALSE)+
  #geom_rect(data = data.frame(xmin=-Inf, xmax=Inf, ymin=0.531, ymax=0.567),
  #          aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='#3B9AB2', alpha=0.2, inherit.aes = FALSE)+
  geom_pointrange(aes(qq_comb, m.log.auc, ymin = m.log.auc-se, ymax = m.log.auc+se, col = epi_endo))+
  geom_point(aes(qq_comb, m.log.auc, fill = epi_endo), size = 2.5, shape = 21)+
  #geom_smooth(aes(qq, log2_AUC, fill = epi_endo, col = epi_endo), method = 'lm')+
  #coord_flip()+
  scale_color_manual(values = c('#74a089', '#9a8822'))+
  scale_fill_manual(values = c('#74a089', '#9a8822'))+
  scale_y_continuous(labels=scaleFUN)+
  scale_shape_manual(values = c(24, 23, 21))+
  labs(title = '', y = 'log2(normalised AUC)', x = 'theoretical quantiles')+
  theme_cu()+
  theme(aspect.ratio = .3)+
  ylim(-1, 0.5)+
  xlim(-2.2,2.2)

sub29 <- subset_cand29 %>% dplyr::select(strain_ID, plant_ID,m.log.auc,sd.log.auc,qq_comb) %>% distinct()
setwd('~/My Drive/protoplast_Mo/protoplastJuly2022')
write_csv(sub29,'sub29.prot.data.csv')

ap29 <- semi_join(ap, cand29_strain_ID, by = 'strain_ID')
apr29 <- semi_join(apr, cand29_strain_ID, by = 'strain_ID')
ap29 %>% filter(replicate_ID == 'replicate4') %>%
  group_by(strain_ID,plant_ID,Time) %>% summarise(spl_min_corr=mean(spl_min_corr),qq_comb) %>%
  ggplot(aes(strain_ID, spl_min_corr, fill=Time,col=Time))+
  geom_point(alpha=.3,size=2,shape=21)+
  scale_colour_viridis_c(name='time (h)')+
  scale_fill_viridis_c(name='time (h)')+
  theme_cu() %+replace% theme(axis.text.x = element_text(angle=0, hjust=1))+
  theme(axis.ticks.x=element_blank())+
  theme(aspect.ratio = 2,
        legend.position = 'none',
        legend.title = element_text(size = 20),
        legend.title.align = 1,
        panel.grid.major = element_blank())+
  labs(x = '', y = '')+
  coord_flip()+
  ylim(0,26000)
apr29 %>% 
  group_by(strain_ID) %>% mutate(se = sd(log.auc)/sqrt(4)) %>%
  ggplot()+
  geom_vline(xintercept=0,col='grey',lwd=.8)+geom_vline(xintercept=-.5,col='grey',lwd=.8)+
  geom_pointrange(aes(m.log.auc,
                  factor(strain_ID,levels=rev(c('229','221','286','369','83','463','854','248',
                                            '103','306','161','690','433','832','315','561',
                                            '449','435','507','754','464','451','759','630',
                                            '732','861','872','761','266'))),
                                    xmin=m.log.auc-sd.log.auc, xmax=m.log.auc+sd.log.auc),shape=21,fill='#e6b254')+
  #geom_point(aes(log.auc,strain_ID))
  theme_void()+
  theme(aspect.ratio = 2,
        legend.position = 'none',
        legend.title = element_text(size = 20),
        legend.title.align = 1,
        #panel.grid.major = element_blank()
        )+
  labs(x = '', y = '')
class(apr29$strain_ID)
apr29 %>% 
  group_by(strain_ID) %>% dplyr::mutate(se = sd(log.auc)/sqrt(n())) %>%
  ggplot()+
  geom_vline(xintercept=0,size=1)+#geom_vline(xintercept=-.5,size=.5)+
  geom_hline(yintercept = 1.5, size = .5)+ geom_hline(yintercept = 2.5, size = .5)+ geom_hline(yintercept = 3.5, size = .5)+ geom_hline(yintercept = 4.5, size = .5)+
  geom_hline(yintercept = 5.5, size = .5)+ geom_hline(yintercept = 6.5, size = .5)+ geom_hline(yintercept = 7.5, size = .5)+ geom_hline(yintercept = 8.5, size = .5)+
  geom_hline(yintercept = 9.5, size = .5)+ geom_hline(yintercept = 10.5, size = .5)+ geom_hline(yintercept = 11.5, size = .5)+ geom_hline(yintercept = 12.5, size = .5)+
  geom_hline(yintercept = 13.5, size = .5)+ geom_hline(yintercept = 14.5, size = .5)+ geom_hline(yintercept = 15.5, size = .5)+ geom_hline(yintercept = 16.5, size = .5)+
  geom_hline(yintercept = 17.5, size = .5)+ geom_hline(yintercept = 18.5, size = .5)+ geom_hline(yintercept = 19.5, size = .5)+ geom_hline(yintercept = 20.5, size = .5)+
  geom_hline(yintercept = 21.5, size = .5)+ geom_hline(yintercept = 22.5, size = .5)+ geom_hline(yintercept = 23.5, size = .5)+ geom_hline(yintercept = 24.5, size = .5)+
  geom_hline(yintercept = 25.5, size = .5)+ geom_hline(yintercept = 26.5, size = .5)+ geom_hline(yintercept = 27.5, size = .5)+ geom_hline(yintercept = 28.5, size = .5)+
  geom_hline(yintercept = 29.5, size = .5)+ #geom_hline(yintercept = 30.5, size = .5)+
  geom_pointrange(aes(m.log.auc,
                      factor(strain_ID,levels=c('464','630','561','248','315','463','872','229','854',
                                                '761','221','451','690','449','306','759','369','832',
                                                '754','266','507','732','861','286','435','433','83',
                                                '161','103')),
                      xmin=m.log.auc-sd.log.auc, xmax=m.log.auc+sd.log.auc),
                  position = position_dodge(width = 0.66), size = 1.3, shape = 21, fill='#FCFC99')+
  #geom_point(aes(log.auc,strain_ID))
  scale_y_discrete(limits=rev)+
  theme_cu()+
  theme(legend.position = 'right',
        aspect.ratio = 3,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.minor.y = element_line(size = 2, colour = 'black'),
        legend.title = element_text(size = 20),
        axis.text = element_text(size = 13),
        #axis.text.y = element_blank(),
        panel.border = element_rect(size = 1, colour = 'black'))+
  labs(x = '', y = '')

aa<-apr29 %>% 
  dplyr::select(m.log.auc,strain_ID) %>% distinct()


factor(data$Treatment, levels=c("Y", "X", "Z"))
  ### heatmap with z-scores to check for reproducibility
  #global heatmap setting
  ht_opt(legend_border = "black",
         heatmap_border = TRUE,
         annotation_border = TRUE)
  
  heat1 <- apr %>% ungroup() %>%
    #filter(plant_ID %in% 'plant1') %>%
    filter(strain_ID != 'mock') %>%
    select(c(strain_ID, replicate_ID, log.auc, qq, epi_endo))
  
  heat_epi_endo1 <- heat1 %>% select(strain_ID, epi_endo) %>% distinct()
  heat_epi_endo1 <- heat_epi_endo1 %>% remove_rownames %>% column_to_rownames((var="strain_ID"))
  heat_epi_endo1 <- as.matrix(heat_epi_endo1)
  heat_epi_endo1 <- heat_epi_endo1[order(as.numeric(row.names(heat_epi_endo1))),]
  
  heat1 <- heat1 %>% select(strain_ID, replicate_ID, log.auc) %>% distinct() %>%
    filter(strain_ID != 'mock') %>% filter(strain_ID != 'SA_0.005') %>% filter(strain_ID != 'SA_0.01')
  heat1 <- spread(heat1, replicate_ID, log.auc)
  heat1 <- heat1 %>% remove_rownames %>% column_to_rownames((var="strain_ID"))
  heat1 <- as.matrix(heat1)
  heat1 <- heat1[order(as.numeric(row.names(heat1))),]
  
  library(circlize)
  col_qq = colorRamp2(c(-2, 0, 2), c("yellow", "white", "purple"))
  
  
  Heatmap(heat1, name = 'log2(normalised AUC))', border = 'black',
          width = ncol(heat1)*unit(10, "mm"), 
          height = nrow(heat1)*unit(.5, "mm"),
          show_row_names = FALSE,
          #cluster_columns = FALSE,
          right_annotation = rowAnnotation(gp = gpar(col = "black"), 
                                           'type' = heat_epi_endo1, 
                                           col = list('type' = c("Endophyte" = "#74a089","Epiphyte" = "#9a8822"))))
  
  
  





