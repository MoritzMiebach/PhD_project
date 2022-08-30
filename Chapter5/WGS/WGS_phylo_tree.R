### Phylogeny with 'ggtree'

install("ggtree")
BiocManager::install("remotes")
BiocManager::install("YuLab-SMU/treedataverse")

library(BiocManager)
library(ggtree)
library(treedataverse)


## Setting working directory
setwd("~/My Drive/WGS") 

## Data import
tr <- read.newick('AConlyTREE.txt')
md <- read.csv("WGS_metadata.csv")
#acn <- read_csv('ACnames.csv')

tr <- as_tibble(tr)

md$label <- md$isolate
md <- md %>% select(-c(isolate, plant_response_num))

tr <- full_join(tr,md, by= 'label')
head(tr)

tr <- as.treedata(tr)

ggtree(tr, layout="daylight")+
  geom_hilight(node=38, fill="#ffe998", alpha=.6)+
  geom_hilight(node=59, fill="#bed1b0", alpha=.6)+
  geom_hilight(node=47, fill="#98e4bb", alpha=.6)+
  geom_hilight(node=32, fill="#fdc3bf", alpha=.6)+
  geom_hilight(node=44, fill="#97d8ff", alpha=.6)+
  geom_tippoint(aes(fill = plant_response), shape = 21, size = 5, alpha = .5)+
  scale_fill_manual(values = c('#f21a00','#e1af00','#3b9ab2'))+
 
ggtree(tr, layout="daylight")+
  geom_hilight(node=38, fill="#ffe998", alpha=.6)+
  geom_hilight(node=59, fill="#bed1b0", alpha=.6)+
  geom_hilight(node=47, fill="#98e4bb", alpha=.6)+
  geom_hilight(node=32, fill="#fdc3bf", alpha=.6)+
  geom_hilight(node=44, fill="#97d8ff", alpha=.6)+
  geom_tippoint(aes(fill = family), shape = 21, size = 5, alpha = .5) 

ggtree(tr
       ,branch.length='none'
       )+
  geom_hilight(node=38, fill="#ffe998", alpha=.5)+
  geom_hilight(node=59, fill="#bed1b0", alpha=.5)+
  geom_hilight(node=47, fill="#98e4bb", alpha=.5)+
  geom_hilight(node=32, fill="#fdc3bf", alpha=.5)+
  geom_hilight(node=44, fill="#97d8ff", alpha=.5)+
  geom_tippoint(aes(fill = family), shape = 21, size = 5, alpha = .8)+
  #geom_treescale(width=0.1)+
  #geom_tiplab(size = 5, family = 'mono', offset = .1, fontface = 1)+
  theme(aspect.ratio = 2.7,
        legend.position='none')
  geom_tiplab(aes(color = plant_response), size = 5, family = 'mono', offset = .1, fontface = 1)+
  #scale_color_manual(values = c('#f21a00','#e1af00','#3b9ab2'))  

ggtree(tr,
       layout="circular")+
  geom_hilight(node=38, fill="#ffe998", alpha=.5)+
  geom_hilight(node=59, fill="#bed1b0", alpha=.5)+
  geom_hilight(node=47, fill="#98e4bb", alpha=.5)+
  geom_hilight(node=32, fill="#fdc3bf", alpha=.5)+
  geom_hilight(node=44, fill="#97d8ff", alpha=.5)+
  geom_tippoint(aes(fill = family), shape = 21, size = 5, alpha = .8)+
  geom_tiplab(size = 5, family = 'mono', offset = .1, fontface = 1)
  theme(aspect.ratio = 2,
        legend.position='bottom')
geom_tiplab(aes(color = plant_response), size = 5, family = 'mono', offset = .1, fontface = 1)+
#scale_color_manual(values = c('#f21a00','#e1af00','#3b9ab2'))  

gheatmap(ggtree(tr, layout="circular"), md$plant_response,offset=.8, width=.2,
         colnames_angle=95, colnames_offset_y = .25)+
  geom_hilight(node=38, fill="#ffe998", alpha=.5)+
  geom_hilight(node=59, fill="#bed1b0", alpha=.5)+
  geom_hilight(node=47, fill="#98e4bb", alpha=.5)+
  geom_hilight(node=32, fill="#fdc3bf", alpha=.5)+
  geom_hilight(node=44, fill="#97d8ff", alpha=.5)+
  geom_tippoint(aes(fill = family), shape = 21, size = 5, alpha = .8)+
  geom_tiplab(aes(color = plant_response), size = 5, family = 'mono', offset = .1, fontface = 1)+
  scale_color_manual(values = c('#f21a00','#e1af00','#3b9ab2'))
  
ggtree(tr, layout="circular")+
  geom_hilight(node=38, fill="#ffe998", alpha=.6)+
  geom_hilight(node=59, fill="#bed1b0", alpha=.6)+
  geom_hilight(node=47, fill="#98e4bb", alpha=.6)+
  geom_hilight(node=32, fill="#fdc3bf", alpha=.6)+
  geom_hilight(node=44, fill="#97d8ff", alpha=.6)+
  geom_text2(aes(subset=!isTip, label=node), hjust=-.3)+ 
  geom_tiplab(aes(color = family))

ggtree(tr, layout = 'daylight')+
  geom_tiplab(aes(color =  family), size = 3, family = 'mono')
                           
