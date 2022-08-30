## Work with DADA2 output

if(!requireNamespace("BiocManager")){
  install.packages("BiocManager")
}
BiocManager::install("phyloseq")

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

write_csv(da2,'dada2.16S.mt+chl.rm.csv') 

da2 <- da2 %>% pivot_longer(cols=endo1:epi5, names_to='sample', values_to='abundance') %>%
  group_by (sample) %>% mutate(rel.abundance = abundance/sum(abundance)) %>%
  mutate(environment = ifelse(sample %in% c('endo1','endo2','endo3','endo4','endo5'), 'endophytic','epiphytic')) %>%
  mutate(plant = ifelse(sample %in% c('endo1', 'epi1'), 'sample 1',
                                ifelse(sample %in% c('endo2', 'epi2'), 'sample 2',
                                               ifelse(sample %in% c('endo3', 'epi3'), 'sample 3',
                                                              ifelse(sample %in% c('endo4', 'epi4'), 'sample 4',
                                                                             'sample 5')))))


### species richness
da2 %>% filter(sample %in% 'endo4') %>% filter(abundance > 0)
### 16S read counts
da2 %>% filter(sample %in% 'endo1') %>% filter(abundance > 0) %>% summarise(reads=sum(abundance))

### reduce ASV level data to Genus level data

dag <- da2 %>% #filter(environment %in% 'epiphytic') %>%
  group_by(Genus, sample) %>% summarise(Kingdom,Phylum,Class,Order,Family,
                                              abundance=sum(abundance)) %>% distinct() %>%
  filter(!is.na(Genus)) %>%
  group_by (sample) %>% mutate(rel.abundance = abundance/sum(abundance)) %>% 
  dplyr::select(-abundance) %>% pivot_wider(names_from = sample,values_from = rel.abundance)

#write_csv(dag,'dada2.16S.rel.abundance.genus.level.csv') 

sub29 <- c(rep('Pseudomonas',times=6),rep('Stenotrophomonas',times=2),'Variovorax','Pelomonas',
           'Massilia',rep('Methylobacterium-Methylorubrum',times=4),'Belnapia',rep('Sphingomonas',times=2),
           'Brevundimonas','Chryseobacterium','Pedobacter','Rathayibacter','Frigoribacterium',
           'Kineococcus','Klenkia','Rhodococcus',rep('Staphylococcus',times=2), 'Bacillus')

sub29 <- as_tibble(sub29) %>% dplyr::rename(Genus=value)
dag29 <- semi_join(dag,sub29,by='Genus')

#write_csv(dag29,'dag29.csv')
dag29s <- read_csv('dag29.sorted.csv')

dag29s %>% filter(Genus %in% c('Pelomonas','Belnapia','Rathayibacter', 
                               'Frigoribacterium','Kinecoccus','Klenkia','Rhodococcus','Bacillus')) %>%
  pivot_longer(7:16,names_to = 'community',values_to = 'relab') %>%
  summarise(max(relab))
  
## heatmap of rel. abundances of the genera that were sequenced
#global heatmap setting
ht_opt(legend_border = "black",
       heatmap_border = TRUE,
       annotation_border = TRUE)

df29 <- as.data.frame(dag29s[,7:16])
col_fun = colorRamp2(c( 0, 0.4), c("white", "black"))
Heatmap(df29, cluster_rows = FALSE, cluster_columns = FALSE,col=col_fun, name='relative abundance')

bai <- read_csv('bai.phyl.dat.csv')
bai %>% group_by(Genus)

## heatmap of rel. abundances of the genera that were sequenced
dag %>% 
  pivot_longer(7:16,names_to = 'community',values_to = 'relab') %>%
  summarise(max(relab))

daf <- as.data.frame(dag[,7:16])
col_fun = colorRamp2(c( 0, 0.4), c("white", "black"))
Heatmap(daf, cluster_rows = FALSE, cluster_columns = FALSE,col=col_fun, name='relative abundance')

## heatmap of presence/absence of the genera that were sequenced
dag2 <- dag %>% 
  pivot_longer(7:16,names_to = 'community',values_to = 'relab') %>%
  mutate(pre=if_else(relab>0,1,0)) %>% dplyr::select(-relab) %>%
  pivot_wider(names_from = community, values_from = pre)

daf2 <- as.data.frame(dag2[,c(3,7:16)])
daf2p <- as.data.frame(dag2[,c(3,12:16)])
col_fun = colorRamp2(c( 0, 1), c("white", "black"))
Heatmap(daf2, cluster_rows = FALSE, cluster_columns = FALSE,col=col_fun, name='relative abundance')

upset(daf2p, nsets = 10,order.by = "freq", nintersects = 70,
      text.scale = 1.5, point.size = 6, line.size = 1.5, mb.ratio = c(0.6, 0.4),
      mainbar.y.label = "ASVs", sets.x.label = "ASVs",
      queries = list(
        list(query = elements, 
             params = list("Phylum", c("Verrucomicrobiota","Sumerlaeota", "Proteobacteria", 'Planctomycetota','Nitrospirota', 'Myxococcota','Gemmatimonadota','Firmicutes','Fibrobacterota','Desulfobacterota','Deinococcota','Cyanobacteria','Chloroflexi','Bdellovibrionota','Bacteroidota','Aquificota','Actinobacteriota','Acidobacteriota','Abditibacteriota')), color = "#fffac8", active = T),
        list(query = elements, 
             params = list("Phylum", c("Verrucomicrobiota","Sumerlaeota", "Proteobacteria", 'Planctomycetota','Nitrospirota', 'Myxococcota','Gemmatimonadota','Firmicutes','Fibrobacterota','Desulfobacterota','Deinococcota','Cyanobacteria','Chloroflexi','Bdellovibrionota','Bacteroidota','Aquificota','Actinobacteriota','Acidobacteriota')), color = "#ffd8b1", active = T),
        list(query = elements, 
             params = list("Phylum", c("Verrucomicrobiota","Sumerlaeota", "Proteobacteria", 'Planctomycetota','Nitrospirota', 'Myxococcota','Gemmatimonadota','Firmicutes','Fibrobacterota','Desulfobacterota','Deinococcota','Cyanobacteria','Chloroflexi','Bdellovibrionota','Bacteroidota','Aquificota','Actinobacteriota')), color = "#fabed4", active = T),
        list(query = elements, 
             params = list("Phylum", c("Verrucomicrobiota","Sumerlaeota", "Proteobacteria", 'Planctomycetota','Nitrospirota', 'Myxococcota','Gemmatimonadota','Firmicutes','Fibrobacterota','Desulfobacterota','Deinococcota','Cyanobacteria','Chloroflexi','Bdellovibrionota','Bacteroidota','Aquificota')), color = "#a9a9a9", active = T),
        list(query = elements, 
             params = list("Phylum", c("Verrucomicrobiota","Sumerlaeota", "Proteobacteria", 'Planctomycetota','Nitrospirota', 'Myxococcota','Gemmatimonadota','Firmicutes','Fibrobacterota','Desulfobacterota','Deinococcota','Cyanobacteria','Chloroflexi','Bdellovibrionota','Bacteroidota')), color = "#f032e6", active = T),
        list(query = elements, 
             params = list("Phylum", c("Verrucomicrobiota","Sumerlaeota", "Proteobacteria", 'Planctomycetota','Nitrospirota', 'Myxococcota','Gemmatimonadota','Firmicutes','Fibrobacterota','Desulfobacterota','Deinococcota','Cyanobacteria','Chloroflexi','Bdellovibrionota')), color = "#911eb4", active = T),
        list(query = elements, 
             params = list("Phylum", c("Verrucomicrobiota","Sumerlaeota", "Proteobacteria", 'Planctomycetota','Nitrospirota', 'Myxococcota','Gemmatimonadota','Firmicutes','Fibrobacterota','Desulfobacterota','Deinococcota','Cyanobacteria','Chloroflexi')), color = "#4363d8", active = T),
        list(query = elements, 
             params = list("Phylum", c("Verrucomicrobiota","Sumerlaeota", "Proteobacteria", 'Planctomycetota','Nitrospirota', 'Myxococcota','Gemmatimonadota','Firmicutes','Fibrobacterota','Desulfobacterota','Deinococcota','Cyanobacteria')), color = "#42d4f4", active = T),
        list(query = elements, 
             params = list("Phylum", c("Verrucomicrobiota","Sumerlaeota", "Proteobacteria", 'Planctomycetota','Nitrospirota', 'Myxococcota','Gemmatimonadota','Firmicutes','Fibrobacterota','Desulfobacterota','Deinococcota')), color = "#3cb44b", active = T),
        list(query = elements, 
             params = list("Phylum", c("Verrucomicrobiota","Sumerlaeota", "Proteobacteria", 'Planctomycetota','Nitrospirota', 'Myxococcota','Gemmatimonadota','Firmicutes','Fibrobacterota','Desulfobacterota')), color = "#bfef45", active = T),
        list(query = elements, 
             params = list("Phylum", c("Verrucomicrobiota","Sumerlaeota", "Proteobacteria", 'Planctomycetota','Nitrospirota', 'Myxococcota','Gemmatimonadota','Firmicutes','Fibrobacterota')), color = "#ffe119", active = T),
        list(query = elements, 
             params = list("Phylum", c("Verrucomicrobiota","Sumerlaeota", "Proteobacteria", 'Planctomycetota','Nitrospirota', 'Myxococcota','Gemmatimonadota','Firmicutes')), color = "#f58231", active = T),
        list(query = elements, 
             params = list("Phylum", c("Verrucomicrobiota","Sumerlaeota", "Proteobacteria", 'Planctomycetota','Nitrospirota', 'Myxococcota','Gemmatimonadota')), color = "#e6194B", active = T),
        list(query = elements, 
             params = list("Phylum", c("Verrucomicrobiota","Sumerlaeota", "Proteobacteria", 'Planctomycetota','Nitrospirota', 'Myxococcota')), color = "#000000", active = T),
        list(query = elements, 
             params = list("Phylum", c("Verrucomicrobiota","Sumerlaeota", "Proteobacteria", 'Planctomycetota','Nitrospirota')), color = "#000075", active = T),
        list(query = elements, 
             params = list("Phylum", c("Verrucomicrobiota","Sumerlaeota", "Proteobacteria", 'Planctomycetota')), color = "#469990", active = T),
        list(query = elements, 
             params = list("Phylum", c("Verrucomicrobiota","Sumerlaeota", "Proteobacteria")), color = "#808000", active = T),
        list(query = elements, 
             params = list("Phylum", c("Verrucomicrobiota","Sumerlaeota")), color = "#9A6324", active = T),
        list(query = elements, 
             params = list("Phylum", "Verrucomicrobiota"), color = '#800000', active = T)))

## Venndiagram comparing epi vs endo per plant
dag2 <- dag2 %>% mutate(asv=1:307)
#
epi1 <- dag2 %>% filter(epi1==1) %>% summarise(epi1=asv)
endo1 <- dag2 %>% filter(endo1==1) %>% summarise(endo1=asv)
plant1 <- list(epi1$epi1, endo1$endo1)
ggVennDiagram(plant1, label_alpha = 0, col = "grey24", label = 'count', edge_size=1)+
  scale_fill_gradient(low="white",high = "grey54")+
  scale_color_manual(values = c('#9a8822','#74a089'))
#
epi2 <- dag2 %>% filter(epi2==1) %>% summarise(epi2=asv)
endo2 <- dag2 %>% filter(endo2==1) %>% summarise(endo2=asv)
plant2 <- list(epi2$epi2, endo2$endo2)
ggVennDiagram(plant2, label_alpha = 0, col = "grey24", label = 'count', edge_size=1)+
  scale_fill_gradient(low="white",high = "grey54")+
  scale_color_manual(values = c('#9a8822','#74a089'))
#
epi3 <- dag2 %>% filter(epi3==1) %>% summarise(epi3=asv)
endo3 <- dag2 %>% filter(endo3==1) %>% summarise(endo3=asv)
plant3 <- list(epi3$epi3, endo3$endo3)
ggVennDiagram(plant3, label_alpha = 0, col = "grey24", label = 'count', edge_size=1)+
  scale_fill_gradient(low="white",high = "grey54")+
  scale_color_manual(values = c('#9a8822','#74a089'))
#
epi4 <- dag2 %>% filter(epi4==1) %>% summarise(epi4=asv)
endo4 <- dag2 %>% filter(endo4==1) %>% summarise(endo4=asv)
plant4 <- list(epi4$epi4, endo4$endo4)
ggVennDiagram(plant4, label_alpha = 0, col = "grey24", label = 'count', edge_size=1)+
  scale_fill_gradient(low="white",high = "grey54")+
  scale_color_manual(values = c('#9a8822','#74a089'))
#
epi5 <- dag2 %>% filter(epi5==1) %>% summarise(epi5=asv)
endo5 <- dag2 %>% filter(endo5==1) %>% summarise(endo5=asv)
plant5 <- list(epi5$epi5, endo5$endo5)
ggVennDiagram(plant5, label_alpha = 0, col = "grey24", label = 'count', edge_size=1)+
  scale_fill_gradient(low="white",high = "grey54")+
  scale_color_manual(values = c('#9a8822','#74a089'))
### graph relative abundance data
da2 %>% 
  group_by(Phylum,environment,plant) %>% #summarise(x=sum(rel.abundance)) %>% filter(x > 0.01)
  filter(sum(rel.abundance)/10 > 0.01) %>% # filters for phyla that represent at least 1% of all samples rel.abundance
  group_by(plant,environment) %>% mutate(total.ab = sum(abundance)) %>% 
  group_by(Phylum,environment,plant) %>% 
  summarise(rel.abundance = sum(abundance)/total.ab,environment=environment,plant=plant) %>% distinct() %>%
   ggplot()+ 
  geom_bar(aes(x=plant, y = rel.abundance, fill = Phylum), stat = 'identity', size=.1, col='black')+ 
  #scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
  #                           '#fdbf6f','#ff7f00','#cab2d6'))
  facet_wrap(~environment, scales="free_x")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = 'right',
        aspect.ratio = 2,
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        panel.border = element_rect(size = 1, colour = 'black'),
        axis.title = element_text(colour = "black", size = 18),
        axis.ticks = element_blank(),
        # legend at the bottom 
        legend.text = element_text(colour = "black", size = 12),
        # facet wrap properties
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black", face = 'italic'),
        strip.text.y = element_text(size = 12, color = "black", face = 'italic'))+
  labs(x='')

t <- da2 %>% 
    group_by(Phylum,environment,plant) %>% summarise(x=sum(rel.abundance)) %>% filter(x > 0.01) %>%
mutate(perc = x*100)
## 9 Phyla with a relative abundance above 1% in at least one of the ten samples
da2 %>% filter(Phylum %in% c('Acidobacteriota','Actinobacteriota','Bacteroidota',
                             'Cyanobacteria','Desulfobacterota','Firmicutes',
                             'Nitrospirota','Proteobacteria','Verrucomicrobiota')) %>%
  group_by(plant,environment) %>% mutate(total.ab = sum(abundance)) %>% 
  group_by(Phylum,environment,plant) %>% 
  summarise(rel.abundance = sum(abundance)/total.ab,environment=environment,plant=plant) %>% distinct() %>%
  ggplot()+ 
  geom_bar(aes(x=plant, y = rel.abundance, fill = Phylum), stat = 'identity', size=.2, col='black')+ 
  scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
                             '#fdbf6f','#cab2d6','#ff7f00'))+
  facet_wrap(~environment, scales="free_x")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = 'right',
        aspect.ratio = 3,
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.border = element_rect(size = 1, colour = 'black'),
        axis.title = element_text(colour = "black", size = 15),
        axis.ticks = element_blank(),
        # legend at the bottom 
        legend.text = element_text(colour = "black", size = 12),
        # facet wrap properties
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black", face = 'italic'),
        strip.text.y = element_text(size = 12, color = "black", face = 'italic'))+
  labs(x='')
## Vorholt colours
da2 %>% filter(Phylum %in% c('Acidobacteriota','Actinobacteriota','Bacteroidota',
                             'Cyanobacteria','Desulfobacterota','Firmicutes',
                             'Nitrospirota','Proteobacteria','Verrucomicrobiota')) %>%
  group_by(plant,environment) %>% mutate(total.ab = sum(abundance)) %>% 
  group_by(Phylum,environment,plant) %>% 
  summarise(rel.abundance = sum(abundance)/total.ab,environment=environment,plant=plant) %>% distinct() %>%
  ggplot()+ 
  geom_bar(aes(x=plant, y = rel.abundance, fill = Phylum), stat = 'identity', size=.2, col='black')+ 
  scale_fill_manual(values=c('#1f78b4','#fdc3bf','#97d8ff','#33a02c','#ff6f6e','#ffe998',
                             '#cab2d6','#d8eebe','#ff7f00'))+
  facet_wrap(~environment, scales="free_x")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = 'right',
        aspect.ratio = 3,
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.border = element_rect(size = 1, colour = 'black'),
        axis.title = element_text(colour = "black", size = 15),
        axis.ticks = element_blank(),
        # legend at the bottom 
        legend.text = element_text(colour = "black", size = 12),
        # facet wrap properties
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black", face = 'italic'),
        strip.text.y = element_text(size = 12, color = "black", face = 'italic'))+
  labs(x='')

## relative abundance of whole dataset combined
da2 %>% filter(Phylum %in% c('Acidobacteriota','Actinobacteriota','Bacteroidota',
                             'Cyanobacteria','Desulfobacterota','Firmicutes',
                             'Nitrospirota','Proteobacteria','Verrucomicrobiota')) %>%
  group_by(sample) %>%
  mutate(total.ab = sum(abundance)) %>% 
  group_by(Phylum, sample) %>% 
  summarise(rel.abundance = sum(abundance)/total.ab*100) %>% distinct()

phyl29 <- da2 %>% filter(Phylum %in% c('Acidobacteriota','Actinobacteriota','Bacteroidota',
                                       'Cyanobacteria','Desulfobacterota','Firmicutes',
                                       'Nitrospirota','Proteobacteria','Verrucomicrobiota')) %>%
  group_by(sample) %>%
  mutate(total.ab = sum(abundance)) %>% 
  group_by(Phylum, sample) %>% 
  summarise(rel.abundance = sum(abundance)/total.ab*100) %>% distinct() %>%
  pivot_wider(names_from = sample,values_from = rel.abundance)

## graph 29 isolate data
d29 <- data.frame()
d29$Phylum <- c('Actinobacteriota','Bacteroidota',
                'Cyanobacteria','Desulfobacterota','Bacillota',
                'Nitrospirota','Proteobacteria','Verrucomicrobiota')
## higher resolution of individual families
da2 %>% filter(Phylum %in% 'Proteobacteria') %>%
  group_by(plant,environment) %>% mutate(total.ab = sum(abundance)) %>% 
  group_by(Class,environment,plant) %>% 
  summarise(rel.abundance = sum(abundance)/total.ab,environment=environment,plant=plant) %>% distinct() %>%
  ggplot()+
  geom_bar(aes(x=plant, y = rel.abundance, fill = Class), stat = 'identity', size=.2, col='black')+ 
  facet_wrap(~environment, scales="free_x")+
  scale_fill_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
                             '#fdbf6f','#cab2d6','#ff7f00'))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = 'right',
        aspect.ratio = 2,
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.border = element_rect(size = 1, colour = 'black'),
        axis.title = element_text(colour = "black", size = 15),
        axis.ticks = element_blank(),
        # legend at the bottom 
        legend.text = element_text(colour = "black", size = 12),
        # facet wrap properties
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black", face = 'italic'),
        strip.text.y = element_text(size = 12, color = "black", face = 'italic'))+
  labs(x='', title = 'Proteobacteria')

da2 %>% filter(Phylum %in% 'Proteobacteria') %>%
  group_by(plant,environment) %>% mutate(total.ab = sum(abundance)) %>% 
  group_by(Class,environment,plant) %>% 
  summarise(rel.abundance = sum(abundance)/total.ab,environment=environment,plant=plant) %>% distinct() %>%
  ggplot()+
  geom_bar(aes(x=plant, y = rel.abundance, fill = Class), stat = 'identity', size=.2, col='black')+ 
  facet_wrap(~environment, scales="free_x")+
  scale_fill_manual(values=c('#98e4bb','#bed1b0'))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = 'right',
        aspect.ratio = 3,
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.border = element_rect(size = 1, colour = 'black'),
        axis.title = element_text(colour = "black", size = 15),
        axis.ticks = element_blank(),
        # legend at the bottom 
        legend.text = element_text(colour = "black", size = 12),
        # facet wrap properties
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black", face = 'italic'),
        strip.text.y = element_text(size = 12, color = "black", face = 'italic'))+
  labs(x='')
## relative abundance of whole dataset combined
da2 %>% filter(Phylum %in% 'Proteobacteria') %>%
  ungroup() %>%
  mutate(total.ab = sum(abundance)) %>% 
  group_by(Class) %>% 
  summarise(rel.abundance = sum(abundance)/total.ab*100) %>% distinct()

prot29 <- da2 %>% filter(Phylum %in% 'Proteobacteria') %>%
  group_by(sample) %>%
  mutate(total.ab = sum(abundance)) %>% 
  group_by(Class,sample) %>% 
  summarise(rel.abundance = sum(abundance)/total.ab*100) %>% distinct() %>%
  pivot_wider(names_from = sample,values_from = rel.abundance)

relab29 <- bind_rows(phyl29,prot29)
write_csv(relab29,'relab29.csv')

### Shannon diversity index (species diversity)
da2w <- da2 %>% dplyr::select(-rel.abundance) %>% dplyr::select(-environment) %>% dplyr::select(-plant) %>%
  pivot_wider(names_from=sample,values_from=abundance) %>% distinct()

shadi <- function(x){ #x is a vector containing the abundances! (not rel.abundances)
  N <- sum(x)
  p <- x/sum(x) # calculates relative abundance
  H <- -sum(p*log(p)) #Shannon diversity index
  return(H)
}

shadi(da2w %>% filter(epi1>0) %>% summarise(epi1))
shadi(da2w %>% filter(epi2>0) %>% summarise(epi2))
shadi(da2w %>% filter(epi3>0) %>% summarise(epi3))
shadi(da2w %>% filter(epi4>0) %>% summarise(epi4))
shadi(da2w %>% filter(epi5>0) %>% summarise(epi5))

shadi(da2w %>% filter(endo1>0) %>% summarise(endo1))
shadi(da2w %>% filter(endo2>0) %>% summarise(endo2))
shadi(da2w %>% filter(endo3>0) %>% summarise(endo3))
shadi(da2w %>% filter(endo4>0) %>% summarise(endo4))
shadi(da2w %>% filter(endo5>0) %>% summarise(endo5))

### Shannon equitability (species evenness) -> ranges from 0-1 where 1 indicates complete eveness

shadi(da2w %>% filter(epi1>0) %>% summarise(epi1))/log(da2w %>% filter(epi1>0) %>% summarise(epi1) %>% summarise(n()))
shadi(da2w %>% filter(epi2>0) %>% summarise(epi2))/log(da2w %>% filter(epi2>0) %>% summarise(epi2) %>% summarise(n()))
shadi(da2w %>% filter(epi3>0) %>% summarise(epi3))/log(da2w %>% filter(epi3>0) %>% summarise(epi3) %>% summarise(n()))
shadi(da2w %>% filter(epi4>0) %>% summarise(epi4))/log(da2w %>% filter(epi4>0) %>% summarise(epi4) %>% summarise(n()))
shadi(da2w %>% filter(epi5>0) %>% summarise(epi5))/log(da2w %>% filter(epi5>0) %>% summarise(epi5) %>% summarise(n()))

shadi(da2w %>% filter(endo1>0) %>% summarise(endo1))/log(da2w %>% filter(endo1>0) %>% summarise(endo1) %>% summarise(n()))
shadi(da2w %>% filter(endo2>0) %>% summarise(endo2))/log(da2w %>% filter(endo2>0) %>% summarise(endo2) %>% summarise(n()))
shadi(da2w %>% filter(endo3>0) %>% summarise(endo3))/log(da2w %>% filter(endo3>0) %>% summarise(endo3) %>% summarise(n()))
shadi(da2w %>% filter(endo4>0) %>% summarise(endo4))/log(da2w %>% filter(endo4>0) %>% summarise(endo4) %>% summarise(n()))
shadi(da2w %>% filter(endo5>0) %>% summarise(endo5))/log(da2w %>% filter(endo5>0) %>% summarise(endo5) %>% summarise(n()))

da2w %>% filter(epi1>0) %>% summarise(epi1) %>% summarise(n())

table.asv <- da2 %>% select(ASV,sample, abundance) %>% pivot_wider(names_from=sample, values_from=abundance) %>%
  tibble::column_to_rownames("ASV") 

table.taxa <- da2 %>% ungroup() %>% select(ASV, Kingdom, Phylum, Class, Order, Family, Genus, Species) %>% 
  distinct() %>% tibble::column_to_rownames("ASV") 

table.samples <- da2 %>% ungroup() %>% select(sample, environment, plant) %>%
  distinct() %>% tibble::column_to_rownames("sample") 

matrix.asv <- as.matrix(table.asv)
matrix.taxa <- as.matrix(table.taxa)

## Transform to phyloseq object

OTU <- otu_table(matrix.asv, taxa_are_rows = TRUE)
TAX <- tax_table(matrix.taxa)
samples <- sample_data(table.samples)

ps <- phyloseq(OTU, TAX, samples)
ps

sample_names(ps)
rank_names(ps)
sample_variables(ps)

rarecurve(t(otu_table(ps)), step=50, cex=0.5)

## normalise number of reads in each sample using median sequencing depth
total <- median(sample_sums(ps))
standf <- function(x, t=total) round(t * (x / sum(x)))
ps <- transform_sample_counts(ps, standf)

ps.phylum <- tax_glom(ps, taxrank="Phylum", NArm=FALSE)
# keep only the 5 most abundant taxa
phylum.sum <- tapply(taxa_sums(ps.phylum), tax_table(ps.phylum)[, "Phylum"], sum, na.rm=TRUE)
top5phyla <- names(sort(phylum.sum, TRUE))[1:5]
ps.phylum <- prune_taxa((tax_table(ps.phylum)[, "Phylum"] %in% top5phyla), ps.phylum)

plot_bar(ps.phylum, fill = "Phylum") + facet_wrap(~environment, scales= "free_x", nrow=1)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = 'right',
        aspect.ratio = 2,
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        panel.border = element_rect(size = 1, colour = 'black'),
        axis.title = element_text(colour = "black", size = 18),
        axis.ticks = element_blank(),
        # legend at the bottom 
        legend.text = element_text(colour = "black", size = 12),
        # facet wrap properties
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black", face = 'italic'),
        strip.text.y = element_text(size = 12, color = "black", face = 'italic'))+
  labs(x='')

  
ps_epi <- subset_samples(ps, environment %in% "epiphytic")
plot_heatmap(ps_epi, method = "NMDS", distance = "bray")
# use only most abundant ASVs for heatmap
ps_abund <- filter_taxa(ps, function(x) sum(x > total*0.20) > 0, TRUE)
ps_abund

plot_heatmap(ps_abund, method = "NMDS", distance = "bray", taxa.label = "Family")
# different distance method and ASVs labeled by Class
plot_heatmap(ps_abund, method = "MDS", distance = "(A+B-2*J)/(A+B-J)", 
             taxa.label = "Class", taxa.order = "Class", 
             trans=NULL, low="beige", high="red", na.value="beige")

## alpha diversity measures

plot_richness(ps, measures=c("Shannon", "Simpson"))
plot_richness(ps, measures=c("Shannon", "Simpson"), x='sample')+
  geom_point(aes(fill=environment),size=5, shape=21)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = 'right',
        aspect.ratio = 2,
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        panel.border = element_rect(size = 1, colour = 'black'),
        axis.title = element_text(colour = "black", size = 18),
        axis.ticks = element_blank(),
        # legend at the bottom 
        legend.text = element_text(colour = "black", size = 12),
        # facet wrap properties
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black", face = 'italic'),
        strip.text.y = element_text(size = 12, color = "black", face = 'italic'))+
  labs(x='')

## Ordination

ps.ord <- ordinate(ps, "NMDS", "bray")

plot_ordination(ps, ps.ord, type="taxa", color="Class", 
                title="ASVs")+
  theme_bw()+
  theme(legend.position = 'right',
        aspect.ratio = 1,
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        panel.border = element_rect(size = 1, colour = 'black'),
        axis.title = element_text(colour = "black", size = 18),
        axis.ticks = element_blank(),
        # legend at the bottom 
        legend.text = element_text(colour = "black", size = 12),
        # facet wrap properties
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black", face = 'italic'),
        strip.text.y = element_text(size = 12, color = "black", face = 'italic'))

plot_ordination(ps, ps.ord, type="taxa", color="Class", 
                title="ASVs") + 
  facet_wrap(~Phylum, 4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        legend.position = 'right',
        aspect.ratio = 1,
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        panel.border = element_rect(size = 1, colour = 'black'),
        axis.title = element_text(colour = "black", size = 18),
        axis.ticks = element_blank(),
        # legend at the bottom 
        legend.text = element_text(colour = "black", size = 12),
        # facet wrap properties
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black", face = 'italic'),
        strip.text.y = element_text(size = 12, color = "black", face = 'italic'))

plot_ordination(ps, ps.ord, type="samples", col="plant", shape="environment", title="Samples")+ 
  geom_point(aes(fill= plant),size=5, col='black')+
  scale_shape_manual(values = c(21,22))+
  theme_bw()+
  theme(legend.position = 'right',
        aspect.ratio = 1.5,
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        panel.border = element_rect(size = 1, colour = 'black'),
        axis.title = element_text(colour = "black", size = 18),
        axis.ticks = element_blank(),
        # legend at the bottom 
        legend.text = element_text(colour = "black", size = 12),
        # facet wrap properties
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black", face = 'italic'),
        strip.text.y = element_text(size = 12, color = "black", face = 'italic'))

## Network analysis
# 'binary' Sørenson distance / 'minimum' Bray-Curtis distance
plot_net(ps, distance = "(A+B-2*J)/(A+B)", type = "taxa", 
         maxdist = 0.8, color="Family") 




