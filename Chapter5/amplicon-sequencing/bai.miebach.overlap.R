# Overlap in culture-independent data from Bai2015 and Miebach2022 enpiphyte data 
# based on OTU/ASVs with rel abundance >0.1%

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

da <- read_csv('Bai.Miebach.comparison.csv')

bai <- da %>% filter(study == 'Bai') %>% distinct() %>% mutate(v=1)
mie <- da %>% filter(study == 'Miebach') %>% distinct() %>% mutate(v=1)

da2 <- bind_rows(bai,mie) %>% pivot_wider(names_from = study, values_from = v)

comp <- list(bai$genus, mie$genus)

ggVennDiagram(comp, label_alpha = 0, col = "grey24", label = 'count', edge_size=2)+
  scale_fill_gradient(low="white",high = "white")+
  scale_color_manual(values = c('#9a8822','#74a089'))



da <- read_csv('Bai.Miebach.family.csv')

bai <- da %>% filter(study == 'Bai') %>% distinct() %>% mutate(v=1)
mie <- da %>% filter(study == 'Miebach') %>% distinct() %>% mutate(v=1)

da2 <- bind_rows(bai,mie) %>% pivot_wider(names_from = study, values_from = v)

comp <- list(bai$family, mie$family)

ggVennDiagram(comp, label_alpha = 0, col = "grey24", label = 'count', edge_size=2)+
  scale_fill_gradient(low="white",high = "white")+
  scale_color_manual(values = c('#9a8822','#74a089'))


da <- read_csv('Bai.Miebach.order.csv')

bai <- da %>% filter(study == 'Bai') %>% distinct() %>% mutate(v=1)
mie <- da %>% filter(study == 'Miebach') %>% distinct() %>% mutate(v=1)

da2 <- bind_rows(bai,mie) %>% pivot_wider(names_from = study, values_from = v)

comp <- list(bai$order, mie$order)

ggVennDiagram(comp, label_alpha = 0, col = "grey24", label = 'count', edge_size=2)+
  scale_fill_gradient(low="white",high = "white")+
  scale_color_manual(values = c('#9a8822','#74a089'))

