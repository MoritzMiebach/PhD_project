## import packages
library(tidyverse)
library(iNEXT)
library(ggplot2)


setwd("~/My Drive/16S/16S.2/trimmed.truseq")

da <- read_csv('dada2.16S2.output.csv')

da2 <- da %>% filter(Family != 'Mitochondria') %>% filter(Order != "Chloroplast") %>%
  filter(!is.na('Phylum')) %>% 
  group_by(ASV) %>% 
  mutate(sum.counts = endo1+endo2+endo3+endo4+endo5+epi1+epi2+epi3+epi4+epi5)
da2 %>% filter(sum.counts == 1)

da.rare <- as.data.frame(da2 %>% ungroup() %>% dplyr::select(endo1, endo2, endo3, endo4, endo5,
                                                             epi1, epi2, epi3, epi4, epi5))

da.r <- iNEXT(da.rare)
ggiNEXT(da.r,facet.var="order")

da.endo <- as.data.frame(da2 %>% ungroup() %>% dplyr::select(endo1, endo2, endo3, endo4, endo5))
da.epi <- as.data.frame(da2 %>% ungroup() %>% dplyr::select(epi1, epi2, epi3, epi4, epi5))


da.ep <- iNEXT(da.epi, q=1, datatype="abundance", knots=80)
ggiNEXT(da.ep, type=1, facet.var="order")+
  #scale_x_continuous(labels = function(x) format(x, scientific = TRUE))+
  #scale_x_continuous(guide = guide_axis(n.dodge=2))+
  scale_x_continuous(breaks=seq(0, 100000, 100000))+
  scale_shape_manual(values=c(16,16,16,16,16))+
  scale_color_manual(values=c('#4053d3','#ddb310','#b51d14','#00beff','#00b25d'))+
  scale_fill_manual(values=c('#4053d3','#ddb310','#b51d14','#00beff','#00b25d'))+
  theme_bw()+
  theme(legend.position = 'right',
        aspect.ratio = 1.5,
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.grid.minor.y = element_blank(),
        title = element_text(size = 15),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        panel.border = element_rect(size = 1, colour = 'black'),
        axis.title = element_text(colour = "black", size = 15),
        axis.ticks = element_blank(),
        # legend at the bottom 
        legend.text = element_text(colour = "black", size = 18))+
  labs(x='sequencing depth', y='species diversity', title='Shannon diversity')#+
  ylim(0,620)

da.ep$AsyEst


da.en <- iNEXT(da.endo, datatype="abundance",q=2, knots=80)
ggiNEXT(da.en, type=1, facet.var="order")+
  #scale_x_continuous(labels = function(x) format(x, scientific = TRUE))+
  #scale_x_continuous(guide = guide_axis(n.dodge=2))+
  scale_x_continuous(breaks=seq(0, 10000, 2000))+
  scale_shape_manual(values=c(16,16,16,16,16))+
  scale_color_manual(values=c('#4053d3','#ddb310','#b51d14','#00beff','#00b25d'))+
  scale_fill_manual(values=c('#4053d3','#ddb310','#b51d14','#00beff','#00b25d'))+
  theme_bw()+
  theme(legend.position = 'right',
        aspect.ratio = 1.5,
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.grid.minor.y = element_blank(),
        title = element_text(size = 15),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 15),
        panel.border = element_rect(size = 1, colour = 'black'),
        axis.title = element_text(colour = "black", size = 15),
        axis.ticks = element_blank(),
        # legend at the bottom 
        legend.text = element_text(colour = "black", size = 18))+
  labs(x='sequencing depth', y='species diversity', title='Simpson diversity')#+
  ylim(0,60)

da.en$AsyEst

class(da.en$AsyEst)

div <- union(da.ep$AsyEst, da.en$AsyEst)

write_csv(div,'iNEXT.diversity.csv')
