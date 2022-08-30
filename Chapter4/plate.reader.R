### Script to analyse the plate reader output of the protoplast response to SA experiment

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

## data import
setwd("~/My Drive/protoplast_Mo/transpose_clean_forR") #new gdrive

SAconc <- read.csv("2021.03.19_SAconc.csv")

# limit data to first 24 h
# protoplast physiology might worsen beyond the 24 h mark

SAconc24 <- SAconc[1:288,]

##wide to long data conversion
dat <- gather(SAconc24, well, fluorescence_intensity, A1:H12)

# split well into plate_column & plate_row
dat_sample <- as.data.frame(str_split(dat$well,"", simplify = TRUE)) # 'simplify' returns a matrix rather than a list 
names(dat_sample)[1] <- 'plate_row'
#names(dat_sample)[2] <- 'plate_column'
dat_sample <- dat_sample %>% unite('plate_column', c('V2','V3'), sep = '', remove = TRUE)

dat <- bind_cols(dat, dat_sample)

# add plant_type info
dat <- mutate(dat, plant_type = ifelse(plate_row %in% c('A', 'B', 'C', 'D'), 'pPR1:nls-YFP',
                                       ifelse(plate_row %in% c('E', 'F', 'G', 'H'), 'Col-0', 'error')))

# add treatment info
dat <- dat %>% arrange(Time, plate_row) %>% add_column(treatment = rep(c('mock', 'mock', '0.005mM_SA', '0.005mM_SA', '0.01mM_SA', '0.01mM_SA', '0.05mM_SA', '0.05mM_SA',
                                                                         '0.1mM_SA', '0.1mM_SA', '0.5mM_SA', '0.5mM_SA','mock', 'mock', '0.005mM_SA', '0.005mM_SA', '0.01mM_SA', '0.01mM_SA', '0.05mM_SA', '0.05mM_SA',
                                                                         '0.1mM_SA', '0.1mM_SA', '0.5mM_SA', '0.5mM_SA', 1:24), times = nrow(dat)/48))

# add info of origin
dat <- dat %>% group_by(plant_type, Time) %>% add_column(epi_endo = rep(c(rep(c('Control'), times = 24), 
                                                                          rep(c('Endophyte'), times = 24)), times = nrow(dat)/48))

## for this script only interested in SA treated samples
dat <- dat %>% filter(epi_endo %in% 'Control')

## baseline subtraction
dat_baseline <- dat %>% group_by(well) %>% filter(Time < 1) %>% summarise(mean_fluo_int_1h = mean(fluorescence_intensity))

dat <- dat %>% group_by(well) %>% right_join(dat_baseline, dat, by = 'well')

dat <- dat %>% filter(Time >=1)

dat <- dat %>% group_by(well) %>% mutate(bs_fluo_int = fluorescence_intensity - mean_fluo_int_1h)

## spline the baseline subtracted data
dat <- dat %>% group_by(well) %>%
  #filter(well %in% 'A1') %>% 
  mutate(spl = smooth.spline(x = Time, y = bs_fluo_int, nknots = 4)$y)

## the area under curve function either subtracts areas in the negative space or adds them.
## however I'm interested in the area under the whole curve (from baseline and not from 0).
## get lowest value of curve to be 0 by adding the min value of to all values for each sample.

dat <- ungroup(dat)
dat <- dat %>% group_by(plant_type) %>% mutate(pos_min = sqrt(min(spl, na.rm = TRUE)^2))
dat <- dat %>% mutate(spl_min_corr = spl + pos_min)
dat <- dat %>% mutate(bs_fluo_corr = bs_fluo_int + pos_min)

## calculate AUC for each well
dat <- dat %>% group_by(well) %>% mutate(area_under_curve = auc(x = Time, y = spl_min_corr, method = 'linear', absolutearea = TRUE, ties = min)) #'linear' means linear interpolation for the trapezoid method

## normalise AUCs against mock control !Do this for each plant_type

mock_AUC_PR1 <- mean(dat  %>% filter(well %in% c('A1', 'A2', 'B1', 'B2')) %>% pull(var = area_under_curve))

mock_AUC_Col <- mean(dat  %>% filter(well %in% c('E1', 'E2', 'F1', 'F2')) %>% pull(var = area_under_curve))

## normalise only against mock
dat_PR1 <- dat %>% filter(plant_type %in% 'pPR1:nls-YFP') %>% group_by(well) %>% mutate(norm_AUC = area_under_curve / mock_AUC_PR1)
dat_Col <- dat %>% filter(plant_type %in% 'Col-0') %>% group_by(well) %>% mutate(norm_AUC = area_under_curve  / mock_AUC_Col)

dat <- bind_rows(dat_PR1, dat_Col)


dat <- dat %>% mutate(SA_conc = ifelse(treatment %in% 'mock', 0,
                                       ifelse(treatment %in% '0.005mM_SA', 0.005,
                                              ifelse(treatment %in% '0.01mM_SA', 0.01,
                                                     ifelse(treatment %in% '0.05mM_SA', 0.05,
                                                            ifelse(treatment %in% '0.1mM_SA', 0.1,
                                                                   ifelse(treatment %in% '0.5mM_SA', 0.5, 1)))))))

dat <- dat %>% mutate(treat = ifelse(treatment %in% 'mock', 'mock',
                                     ifelse(treatment %in% '0.005mM_SA', '0.005 mM SA',
                                            ifelse(treatment %in% '0.01mM_SA', '0.01 mM SA',
                                                   ifelse(treatment %in% '0.05mM_SA', '0.05 mM SA',
                                                          ifelse(treatment %in% '0.1mM_SA', '0.1 mM SA',
                                                                 ifelse(treatment %in% '0.5mM_SA', '0.5 mM SA', 1)))))))

dat$treat <- factor(dat$treat, levels = c('mock', '0.005 mM SA', '0.01 mM SA', '0.05 mM SA', '0.1 mM SA', '0.5 mM SA'))

## graph data
# custom ggplot theme
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
          legend.title = element_text(size =0),
          legend.text = element_text(colour = "black", family = "sans", size = 15),
          # plot title
          plot.title = element_text(colour = "black", family = "sans", size = 20),
          aspect.ratio = 1
    )
}

theme_cu_animation <- function() {
  theme_bw() +
    theme(panel.border = element_rect(colour = "grey54", linetype = 'solid', size = 1),
          panel.grid.major = element_line(colour = "grey54", linetype = 'solid', size = 0.5),
          panel.grid.minor = element_line(colour = "grey54", linetype = 'solid', size = 0.5),
          axis.text = element_text(colour = "black", family = "sans", size = 25),
          axis.title = element_text(colour = "black", family = "sans", size = 35),
          axis.ticks = element_blank(),
          # legend at the bottom 
          legend.justification = "center",
          legend.position = "bottom",
          legend.title = element_text(size =0),
          legend.text = element_text(colour = "black", family = "sans", size = 25),
          # plot title
          plot.title = element_text(colour = "black", family = "sans", size = 35),
          aspect.ratio = 1
    )
}

## response curve animation
p <- dat  %>% filter(well %in% 'A1') %>% group_by(Time, treat, plant_type) %>% 
  summarise(mean_fluo = mean(spl_min_corr), bs_fluo_corr = bs_fluo_corr) %>%
  ggplot(aes(x = Time, y = mean_fluo))+
  #geom_smooth()+
  geom_line(lwd = 3, col = 'grey24')+
  geom_point(aes(x = Time, y = bs_fluo_corr, group = seq_along(Time)), shape = 21, size = 2, fill = 'grey24' )+
  scale_fill_manual(values = c('grey24', '#78B7C5','#3B9AB2', '#EBCC2A', '#E1AF00', '#F21A00'))+
  scale_color_manual(values = c('grey24', '#78B7C5','#3B9AB2', '#EBCC2A', '#E1AF00', '#F21A00'))+
  theme_cu_animation()+
  labs(title = '', y = 'fluorescence intensity (a.u.)', x = 'time (h)')+
  ylim(0,20000)+
  transition_reveal(Time)

## plots
# show method on example of mock curve
dat  %>% filter(well %in% 'A1') %>% group_by(Time, treat, plant_type) %>% 
  summarise(mean_fluo = mean(spl), fluo_int = fluorescence_intensity) %>%
  ggplot(aes(x = Time, y = mean_fluo))+
  #geom_smooth()+
  geom_hline(yintercept = 0, size = 1.5)+
  geom_point(aes(x = Time, y = fluo_int), shape = 21, size = 2.5, fill = "grey24", alpha = .8)+
  scale_fill_manual(values = c('grey24', '#78B7C5','#3B9AB2', '#EBCC2A', '#E1AF00', '#F21A00'))+
  scale_color_manual(values = c('grey24', '#78B7C5','#3B9AB2', '#EBCC2A', '#E1AF00', '#F21A00'))+
  theme_cu()+
  theme(legend.position = 'right',
        aspect.ratio = 1.5,
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        panel.border = element_rect(size = 2, colour = 'black'))+
  labs(y = 'fluorescence intensity (a.u.)', x = 'time (h)', title = 'raw')
  ylim(-5000,10000)


dat  %>% filter(well %in% 'A1') %>% group_by(Time, treat, plant_type) %>% 
  summarise(mean_fluo = mean(spl), bs_fluo_int = bs_fluo_int) %>%
  ggplot(aes(x = Time, y = mean_fluo))+
  #geom_smooth()+
  geom_hline(yintercept = 0, size = 1.5)+
  geom_point(aes(x = Time, y = bs_fluo_int), shape = 21, size = 2.5, fill = "grey24", alpha = .8)+
  scale_fill_manual(values = c('grey24', '#78B7C5','#3B9AB2', '#EBCC2A', '#E1AF00', '#F21A00'))+
  scale_color_manual(values = c('grey24', '#78B7C5','#3B9AB2', '#EBCC2A', '#E1AF00', '#F21A00'))+
  theme_cu()+
  theme(legend.position = 'right',
        aspect.ratio = 1.5,
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        panel.border = element_rect(size = 2, colour = 'black'))+
  labs(y = 'fluorescence intensity (a.u.)', x = 'time (h)', title = 'background subtracted')+
  ylim(-5000,10000)

dat  %>% filter(well %in% 'A1') %>% group_by(Time, treat, plant_type) %>% 
  summarise(mean_fluo = mean(spl), bs_fluo_int = bs_fluo_int) %>%
  ggplot(aes(x = Time, y = mean_fluo))+
  geom_bar(stat = "identity", size=1, col = "#f4d984")+
  geom_hline(yintercept = 0, size = 1.5)+
  geom_line(lwd = 2)+
  scale_fill_manual(values = c('grey24', '#78B7C5','#3B9AB2', '#EBCC2A', '#E1AF00', '#F21A00'))+
  scale_color_manual(values = c('grey24', '#78B7C5','#3B9AB2', '#EBCC2A', '#E1AF00', '#F21A00'))+
  theme_cu()+
  theme(legend.position = 'right',
        aspect.ratio = 1.5,
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        panel.border = element_rect(size = 2, colour = 'black'))+
  labs(y = 'fluorescence intensity (a.u.)', x = 'time (h)', title = 'spline fitted')+
  ylim(-5000,10000)

dat  %>% filter(well %in% 'A1') %>% group_by(Time, treat, plant_type) %>% 
  summarise(mean_fluo = mean(spl_min_corr), bs_fluo_corr = bs_fluo_corr) %>%
  ggplot(aes(x = Time, y = mean_fluo))+
  geom_bar(stat = "identity", size=1, col = "#f4d984")+
  geom_hline(yintercept = 0, size = 1.5)+
  geom_line(lwd = 2)+
  scale_fill_manual(values = c('grey24', '#78B7C5','#3B9AB2', '#EBCC2A', '#E1AF00', '#F21A00'))+
  scale_color_manual(values = c('grey24', '#78B7C5','#3B9AB2', '#EBCC2A', '#E1AF00', '#F21A00'))+
  theme_cu()+
  theme(legend.position = 'right',
        aspect.ratio = 1.5,
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        panel.border = element_rect(size = 2, colour = 'black'))+
  labs(y = 'fluorescence intensity (a.u.)', x = 'time (h)', title = 'raised')+
  ylim(-5000,10000)

dat  %>% filter(plant_type %in% 'Col-0') %>% group_by(Time, treat, plant_type) %>% summarise(mean_fluo = mean(spl_min_corr)) %>%
  ggplot(aes(x = Time, y = mean_fluo, col = treat, fill = treat))+
  #geom_smooth()+
  geom_line(lwd = 2)+
  scale_fill_manual(values = c('grey24', '#78B7C5','#3B9AB2', '#EBCC2A', '#E1AF00', '#F21A00'))+
  scale_color_manual(values = c('grey24', '#78B7C5','#3B9AB2', '#EBCC2A', '#E1AF00', '#F21A00'))+
  theme_cu()+
  theme(legend.position = 'right',
        aspect.ratio = 1.5,
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 20),
        panel.border = element_rect(size = 2, colour = 'black'))+
  labs(title = 'Col-0', y = 'fluorescence intensity (a.u.)', x = 'time (h)')+
  ylim(0,15000)

dat  %>% filter(plant_type %in% 'pPR1:nls-YFP') %>% group_by(Time, treat, plant_type) %>% summarise(mean_fluo = mean(spl_min_corr)) %>%
  ggplot(aes(x = Time, y = mean_fluo, col = treat, fill = treat))+
  #geom_smooth()+
  geom_line(lwd = 2)+
  scale_fill_manual(values = c('grey24', '#78B7C5','#3B9AB2', '#EBCC2A', '#E1AF00', '#F21A00'))+
  scale_color_manual(values = c('grey24', '#78B7C5','#3B9AB2', '#EBCC2A', '#E1AF00', '#F21A00'))+
  theme_cu()+
  theme(legend.position = 'right',
        aspect.ratio = 1.5,
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 20),
        panel.border = element_rect(size = 2, colour = 'black'))+
  labs(title = 'pPR1:eYFP-nls', y = 'fluorescence intensity (a.u.)', x = 'time (h)')+
  ylim(0,15000)

dat %>% filter(plant_type %in%'Col-0') %>% 
  ggplot(aes(treat, spl_min_corr, fill=Time,col=Time))+
  geom_point(alpha=.3,size=2,shape=21)+
  scale_colour_viridis_c(name='time (h)')+
  scale_fill_viridis_c(name='time (h)')+
  theme_cu() %+replace% theme(axis.text.x = element_text(angle=0, hjust=1))+
  theme(axis.ticks.x=element_blank())+
  theme(aspect.ratio = 2,
        legend.position = 'right',
        legend.title = element_text(size = 20),
        legend.title.align = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(title = 'Col-0', x = 'isolate #', y = 'fluorescence intensity (a.u.)')+
  coord_flip()+
  ylim(0,35000)

dat %>% filter(plant_type %in%'pPR1:nls-YFP') %>% 
  ggplot(aes(treat, spl_min_corr, fill=Time,col=Time))+
  geom_point(alpha=.3,size=2,shape=21)+
  scale_colour_viridis_c(name='time (h)')+
  scale_fill_viridis_c(name='time (h)')+
  theme_cu() %+replace% theme(axis.text.x = element_text(angle=0, hjust=1))+
  theme(axis.ticks.x=element_blank())+
  theme(aspect.ratio = 2,
        legend.position = 'right',
        legend.title = element_text(size = 20),
        legend.title.align = 1,
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  labs(title = 'pPR1:eYFP-nls', x = 'isolate #', y = 'fluorescence intensity (a.u.)')+
  coord_flip()+
  ylim(0,35000)

## plot AUCs
dat %>% 
  group_by(treat, SA_conc, plant_type) %>% summarise(log2AUC = log2(area_under_curve)) %>% mutate(meanlog = mean(log2AUC), sd = sd(log2AUC)) %>%
  distinct() %>%
  ggplot()+
  #geom_hline(yintercept = 0, size = 1, color = '#F21A00')+
  geom_pointrange(aes(treat, meanlog, shape = plant_type, ymin = meanlog-sd, ymax = meanlog+sd, group = plant_type, fill = treat),
                  size = 1.8, position = position_dodge(width = 0.5))+
  geom_point(aes(treat, log2AUC, group = plant_type, shape = plant_type), position = position_dodge(width = 0.5), size = 3, fill = 'grey24', alpha =.6)+
  scale_fill_manual(values = c('grey24', '#78B7C5','#3B9AB2', '#EBCC2A', '#E1AF00', '#F21A00'))+
  scale_shape_manual(values = c(21, 22))+
  theme_cu()+
  theme(legend.position = 'right',
        aspect.ratio = 1.5,
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 20),
        panel.border = element_rect(size = 2, colour = 'black'),
        axis.text.x = element_blank()
        )+
  labs(x = '', y = 'log2(AUC/AUC[mock])', color = 'plant type')+
  guides(fill = 'none')+
  theme(legend.key.width=unit(1,"cm"))

## reference plot for protoplast+bacteria experiment
dat %>% filter(treatment %in% '0.01mM_SA') %>%
  filter(plant_type %in% 'pPR1:nls-YFP') %>%
  group_by(treat, SA_conc, plant_type) %>%  summarise(log2AUC = log2(norm_AUC)) %>% 
  distinct() %>%
  mutate(meanlog = mean(log2AUC), sd = sd(log2AUC)) %>%
  mutate(meanlogpse = meanlog+(sd/4), meanlogmse = meanlog-(sd/4))


dat %>% filter(treatment %in% c('mock', '0.005mM_SA', '0.01mM_SA')) %>%
  filter(plant_type %in% 'pPR1:nls-YFP') %>%
  group_by(treat, SA_conc, plant_type) %>%  summarise(log2AUC = log2(norm_AUC)) %>% 
  distinct() %>%
  mutate(meanlog = mean(log2AUC), sd = sd(log2AUC)) %>%
  ggplot()+
  geom_hline(yintercept = 0, size = 1, color = '#F21A00')+
  geom_hline(yintercept = 0.343, size = 1, color = '#78B7C5')+
  geom_hline(yintercept = 0.549, size = 1, color = '#3B9AB2')+
  geom_rect(data = data.frame(xmin=-Inf, xmax=Inf, ymin=-0.0294, ymax=0.0236),
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='#f21a00', alpha=0.2, inherit.aes = FALSE)+
  geom_rect(data = data.frame(xmin=-Inf, xmax=Inf, ymin=0.327, ymax=0.359),
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='#78B7C5', alpha=0.2, inherit.aes = FALSE)+
  geom_rect(data = data.frame(xmin=-Inf, xmax=Inf, ymin=0.531, ymax=0.567),
            aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill='#3B9AB2', alpha=0.2, inherit.aes = FALSE)+
  geom_pointrange(aes(treat, meanlog, fill = treat, ymin = meanlog-sd, ymax = meanlog+sd, group = plant_type),
                  size = 1.8, shape = 22, position = position_dodge(width = 0.35))+
  geom_point(aes(treat, log2AUC, group = plant_type), shape = 22, position = position_dodge(width = 0.35), size = 3, alpha=.6, fill = 'grey24')+
  scale_fill_manual(values = c('grey24', '#78B7C5','#3B9AB2', '#EBCC2A', '#E1AF00', '#F21A00'))+
  theme_cu()+
  theme(legend.position = 'right',
        aspect.ratio = 2.2,
        #panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text = element_text(size = 20),
        panel.border = element_rect(size = 2, colour = 'black'),
        axis.text.x = element_blank())+
  ylim(-0.5, 0.75)+
  labs(x = '', y = 'log2(AUC/AUC[mock])', color = 'plant type')
scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  theme(legend.key.width=unit(1,"cm"))

## statistics
# Kruskal Wallis test is non-parametric test. 
# It thus tests for significance in difference between medians.
dat.pr1 <- dat %>% filter(plant_type %in% 'pPR1:nls-YFP') %>%
  group_by(treat, SA_conc, plant_type) %>% summarise(log2AUC = log2(area_under_curve)) %>% distinct()
kruskal.test(log2AUC ~ treat, data = dat.pr1)
dunn.test(dat.pr1$log2AUC, dat.pr1$treat, method = 'sidak')

dat.wt <- dat %>% filter(plant_type %in% 'Col-0') %>%
  group_by(treat, SA_conc, plant_type) %>% summarise(log2AUC = log2(area_under_curve)) %>% distinct()
kruskal.test(log2AUC ~ treat, data = dat.wt)
dunn.test(dat.wt$log2AUC, dat.wt$treat, method = 'sidak')

# One-way ANOVA instead of Kruskal Wallis test.
summary(aov(log2AUC ~ treat, data = dat.wt))
# In one-way ANOVA test, a significant p-value indicates that some of the group means are different.
# Tukey test for multiple pairwise-comparison.
TukeyHSD(aov(log2AUC ~ treat, data = dat.wt))
#
summary(aov(log2AUC ~ treat, data = dat.pr1))
TukeyHSD(aov(log2AUC ~ treat, data = dat.pr1))


dat.pr1.sub <- dat %>% filter(treat %in% c('mock', '0.005 mM SA', '0.01 mM SA')) %>% 
  group_by(treat, SA_conc, plant_type) %>% summarise(log2AUC = log2(area_under_curve)) %>% distinct()
kruskal.test(log2AUC ~ treat, data = dat.pr1.sub)
dunn.test(dat.pr1$log2AUC, dat.pr1$treat, method = 'sidak')

dat.m <- dat %>% filter(treat %in% 'mock') %>% group_by(treat, SA_conc, plant_type) %>%
  summarise(log2AUC = log2(area_under_curve)) %>% distinct()
pairwise.wilcox.test(dat.m$log2AUC, dat.m$plant_type, p.adjust.method = 'none')
dat.1 <- dat %>% filter(treat %in% '0.005 mM SA') %>% group_by(treat, SA_conc, plant_type) %>%
  summarise(log2AUC = log2(area_under_curve)) %>% distinct()
pairwise.wilcox.test(dat.1$log2AUC, dat.1$plant_type, p.adjust.method = 'none')
dat.2 <- dat %>% filter(treat %in% '0.01 mM SA') %>% group_by(treat, SA_conc, plant_type) %>%
  summarise(log2AUC = log2(area_under_curve)) %>% distinct()
pairwise.wilcox.test(dat.2$log2AUC, dat.2$plant_type, p.adjust.method = 'none')
dat.3 <- dat %>% filter(treat %in% '0.05 mM SA') %>% group_by(treat, SA_conc, plant_type) %>%
  summarise(log2AUC = log2(area_under_curve)) %>% distinct()
pairwise.wilcox.test(dat.3$log2AUC, dat.3$plant_type, p.adjust.method = 'none')
dat.4 <- dat %>% filter(treat %in% '0.1 mM SA') %>% group_by(treat, SA_conc, plant_type) %>%
  summarise(log2AUC = log2(area_under_curve)) %>% distinct()
pairwise.wilcox.test(dat.4$log2AUC, dat.4$plant_type, p.adjust.method = 'none')
dat.5 <- dat %>% filter(treat %in% '0.5 mM SA') %>% group_by(treat, SA_conc, plant_type) %>%
  summarise(log2AUC = log2(area_under_curve)) %>% distinct()
pairwise.wilcox.test(dat.5$log2AUC, dat.5$plant_type, p.adjust.method = 'none')


################################################################################
dat %>% summarise(well = well, treatment = treatment) %>% distinct() %>% filter(treatment %in% 'mock')

export_dat <- dat %>% 
  group_by(treat, SA_conc, plant_type) %>% summarise(log2normAUC = log2(norm_AUC), area_under_curve = area_under_curve) %>%
  mutate(meanlog2normAUC = mean(log2normAUC), selog2normAUC = sd(log2normAUC)/sqrt(4), log2AUC = log2(area_under_curve)) %>%
  distinct()

setwd("/Volumes/GoogleDrive/My Drive/protoplast_Mo/SA_response")
write.csv(export_dat, 'protoplast_assay_SA_summary.csv')

## all curves
dat  %>% filter(plant_type %in% 'pPR1:nls-YFP') %>% group_by(Time, treat, plant_type, well) %>% 
  summarise(spl_min_corr = spl_min_corr) %>%
  ggplot()+
  #geom_smooth()+
  geom_line(aes(x = Time, y = spl_min_corr, col = treat, group = well), lwd = 1)+
  scale_fill_manual(values = c('grey24','#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00'))+
  scale_color_manual(values = c('grey24','#3B9AB2', '#78B7C5', '#EBCC2A', '#E1AF00', '#F21A00'))+
  theme_cu()+
  labs(y = 'fluorescence intensity (a.u.)', x = 'time (h)')+
  ylim(0,16000)


