## Overlap of epiphytic and endophytic ASVs

library(tidyverse)
library(phyloseq); packageVersion("phyloseq")
library(Biostrings); packageVersion("Biostrings")
library(vegan)
library(ComplexHeatmap)
library(viridis)
library(circlize)
library(UpSetR)
library("ggVennDiagram")

setwd("~/My Drive/16S/16S.2/trimmed.truseq")

da <- read_csv('dada2.16S2.output.csv')

da2 <- da %>% filter(Family != 'Mitochondria') %>% filter(Order != "Chloroplast") %>%
  filter(!is.na('Phylum'))

## converet to presence/absence data
dapa <- da2 %>% 
  pivot_longer(9:18,names_to = 'community',values_to = 'relab') %>%
  mutate(pre=if_else(relab>0,1,0)) %>% dplyr::select(-relab) %>%
  pivot_wider(names_from = community, values_from = pre) %>% mutate(asv=1:1447)
 

## Venn diagram based on ASVs
epi1 <- dapa %>% filter(epi1 %in% 1) %>% summarise(epi1=asv)
endo1 <- dapa %>% filter(endo1 %in% 1) %>% summarise(endo1=asv)
plant1 <- list(epi1$epi1, endo1$endo1)
ggVennDiagram(plant1, label_alpha = 0, col = "grey24", label = 'count', edge_size=2)+
  scale_fill_gradient(low="white",high = "white")+
  scale_color_manual(values = c('#9a8822','#74a089'))

epi2 <- dapa %>% filter(epi2 %in% 1) %>% summarise(epi2=asv)
endo2 <- dapa %>% filter(endo2 %in% 1) %>% summarise(endo2=asv)
plant2 <- list(epi2$epi2, endo2$endo2)
ggVennDiagram(plant2, label_alpha = 0, col = "grey24", label = 'count', edge_size=2)+
  scale_fill_gradient(low="white",high = "white")+
  scale_color_manual(values = c('#9a8822','#74a089'))

epi3 <- dapa %>% filter(epi3 %in% 1) %>% summarise(epi3=asv)
endo3 <- dapa %>% filter(endo3 %in% 1) %>% summarise(endo3=asv)
plant3 <- list(epi3$epi3, endo3$endo3)
ggVennDiagram(plant3, label_alpha = 0, col = "grey24", label = 'count', edge_size=2)+
  scale_fill_gradient(low="white",high = "white")+
  scale_color_manual(values = c('#9a8822','#74a089'))

epi4 <- dapa %>% filter(epi4 %in% 1) %>% summarise(epi4=asv)
endo4 <- dapa %>% filter(endo4 %in% 1) %>% summarise(endo4=asv)
plant4 <- list(epi4$epi4, endo4$endo4)
ggVennDiagram(plant4, label_alpha = 0, col = "grey24", label = 'count', edge_size=2)+
  scale_fill_gradient(low="white",high = "white")+
  scale_color_manual(values = c('#9a8822','#74a089'))

epi5 <- dapa %>% filter(epi5 %in% 1) %>% summarise(epi5=asv)
endo5 <- dapa %>% filter(endo5 %in% 1) %>% summarise(endo5=asv)
plant5 <- list(epi5$epi5, endo5$endo5)
ggVennDiagram(plant5, label_alpha = 0, col = "grey24", edge_size=2,label = 'count')+
  scale_fill_gradient(low="white",high = "white")+
  scale_color_manual(values = c('#9a8822','#74a089'))

## Venn diagram based on all samples combined
dapa2 <- dapa %>% pivot_longer(cols=endo1:epi5, names_to='sample', values_to='abundance') %>%
  mutate(ep.en = if_else(sample %in% c('endo1','endo2','endo3','endo4','endo5'),'endo','epi'))

epi <- dapa2 %>% filter(ep.en %in% 'epi') %>% filter(abundance==1) %>% summarise(epi=asv) %>% distinct()
endo <- dapa2 %>% filter(ep.en %in% 'endo') %>% filter(abundance==1) %>% summarise(endo=asv) %>% distinct()
plant <- list(epi$epi, endo$endo)
ggVennDiagram(plant, label_alpha = 0, col = "grey24", edge_size=2, label = 'count')+
  scale_fill_gradient(low="white",high = "white")+
  scale_color_manual(values = c('#9a8822','#74a089'))

## reduce ASV level data to Genus level data
dag <- da2 %>% pivot_longer(9:18,names_to = 'community',values_to = 'relab') %>%
  group_by(Genus,community) %>% summarise(Kingdom,Phylum,Class,Order,Family,
                                        relab=sum(relab)) %>% distinct() %>%
  filter(!is.na(Genus)) %>%
  mutate(pre=if_else(relab>0,1,0)) %>% dplyr::select(-relab) %>%
  pivot_wider(names_from = community, values_from = pre) %>% ungroup() %>%
  mutate(asv=1:307)


## Venn diagram based on genera
dag2 <- dag %>% pivot_longer(cols=endo1:epi5, names_to='sample', values_to='abundance') %>%
  mutate(ep.en = if_else(sample %in% c('endo1','endo2','endo3','endo4','endo5'),'endo','epi'))

epi <- dag2 %>% filter(ep.en %in% 'epi') %>% filter(abundance==1) %>% summarise(epi=asv) %>% distinct()
endo <- dag2 %>% filter(ep.en %in% 'endo') %>% filter(abundance==1) %>% summarise(endo=asv) %>% distinct()
plant <- list(epi$epi, endo$endo)
ggVennDiagram(plant, label_alpha = 0, col = "grey24", edge_size=2, label = 'count')+
  scale_fill_gradient(low="white",high = "white")+
  scale_color_manual(values = c('#9a8822','#74a089'))

