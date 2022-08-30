### This script explores functional profiles (SEED, COG, PFAM, TIGRFAM domains) of 30 bacterial genomes.
### The input data was generated in kbase with the "functional domains" app.
### Genomes were assembled with RAST and annotated with the "annotation" app.

library(ComplexHeatmap)
library(tidyverse)
library(UpSetR)
library(RColorBrewer)
library(circlize)
library(viridis)

## Setting working directory
setwd("~/My Drive/WGS") 

## Data import
sr <- read_csv("SEED_raw_count.csv")
co <- read_csv("COG_raw_count.csv")
pf <- read_csv("PFAM_raw_count.csv")
tr <- read_csv("TIGRFAM_raw_count.csv")
md <- read.csv("WGS_metadata.ac.csv")

## Transform data into 'Heatmap' readable matrix
TH <- function(x){
  domains <- as.data.frame(x[,1])
  x <- as.data.frame(x[,2:30])
  row.names(x) <- domains$domain
  x <- as.matrix(x)
  return(x)
}

srm <- TH(sr)
com <- TH(co)
pfm <- TH(pf)
trm <- TH(tr)

## Transform data into presence/absence data
PA <- function(x){
  x[x > 0] <- 1
  return(x)
}

srp <- PA(srm)
cop <- PA(com)
pfp <- PA(pfm)
trp <- PA(trm)

## Transform data into z-scores of raw counts
ZS <- function(x){
  x <- na.omit(t(scale(t(x))))
  return(x)
}

srz <- ZS(srm)
coz <- ZS(com)
pfz <- ZS(pfm)
trz <- ZS(trm)

## Heatmap annotation
prot_col <- colorRamp2(c(-0.8,0), c("#FFFFF0", "#FFBF00"))
prot_col <- colorRamp2(c(-0.8,0), c("blue", "white"))
ha = HeatmapAnnotation(classification = md$classification,
                       protoplast.response = md$m.log.auc,
                       family = anno_text(md$symbol,location = 0.5, just = "center",
                                          gp = gpar(border = "grey24")),
                       col = list(classification = c("Gammaproteobacteria" = "#bed1b0", "Bacillota" = "#ffe998", 
                                                     "Alphaproteobacteria" = "#98e4bb", "Bacteroidota" = "#97d8ff",
                                                     "Actinomycetota" = "#fdc3bf"),
                                  protoplast.response=prot_col),
                       gp = gpar(col = "grey24"))

## Heatmap - raw count data
Heatmap(srm, col = rev(magma(100)), bottom_annotation = ha, name = 'raw_count', show_row_names = FALSE)
Heatmap(com, col = rev(rainbow(10)), bottom_annotation = ha, name = 'raw_count', show_row_names = FALSE)
Heatmap(pfm, col = rev(rainbow(10)), bottom_annotation = ha, name = 'raw_count', show_row_names = FALSE)
Heatmap(trm, col = rev(rainbow(10)), bottom_annotation = ha, name = 'raw_count', show_row_names = FALSE)

## Heatmap - z score
Heatmap(srz,bottom_annotation=ha,name='z-score',border=TRUE,show_row_names=FALSE)
Heatmap(coz,bottom_annotation=ha,name='z-score',border=TRUE,show_row_names=FALSE)
Heatmap(pfz,bottom_annotation=ha,name='z-score',border=TRUE,show_row_names=FALSE)
Heatmap(trz,bottom_annotation=ha,name='z-score',border=TRUE,show_row_names=FALSE)

## Heatmap presence/absence
Heatmap(srp, bottom_annotation = ha, name = 'presence/absence', border = TRUE,
        show_row_names = FALSE,
        col = c('grey89', 'grey24'))

## Determine functional domains that might explain plant response
## Correlate each functional domain against plant response data (Spearman correlation)
CF <- function(x){
  df.cor <- data.frame()
  counter <- nrow(x)
  for(i in 1:counter){
    df.cor[i,1] <- cor.test(md$m.log.auc, x[i,])$estimate
    # positive values indicate positive correlation, negative values indicate negative correlation
    df.cor[i,2] <- cor.test(md$m.log.auc, x[i,])$p.value
    # a significant p-value for presence of domain y means that it will give us a reliable guess of the plant response.
  }
  df.cor <- as_tibble(df.cor)
  df.cor <- df.cor %>% dplyr::rename(rho = V1, p.value = V2)
  df.cor$domain <- row.names(x)
  return(df.cor)
}

srz.cor <- CF(srz)
coz.cor <- CF(coz)
pfz.cor <- CF(pfz)
trz.cor <- CF(trz)

srm.cor <- CF(srm)

srp.cor <- CF(srp)
cor.test(md$m.log.auc,srz[262,])$estimate

plot(md$m.log.auc,srz[262,])

## plot significant correlation data
srp.cor %>% 
  filter(p.value < 0.05) %>%
  mutate(nlogp=-1*log10(p.value)) %>%
  ggplot(aes(y = reorder(domain, rho)))+
  geom_bar(aes(x = rho, fill=nlogp),stat = "identity", col = 'grey24', width=.01)+
  geom_point(aes(x = rho, fill=nlogp),stat = "identity", col = 'grey24', shape=21, size=4)+
 scale_fill_viridis_c(option='A', direction=-1)+
  theme_bw()+
  theme(panel.border = element_rect(size = 2, colour = 'black'),
        panel.grid.minor = element_blank(),
        panel.grid.minor.y = element_line(size = 2, colour = 'black'),
        axis.text = element_text(colour = "black", size = 15),
        axis.title = element_text(colour = "black", size = 18),
        legend.justification = "center",
        legend.position = "right",
        legend.title = element_text(size =15),
        legend.text = element_text(colour = "black", size = 15),
        plot.title = element_text(colour = "black", size = 18),
        aspect.ratio = 2.8)+
  labs(y='', x='rho', title='')+
  xlim(-0.8,0.8)

## Heatmap of only the domains that seem to have an influence on the plant response
md <- md %>% arrange(m.log.auc) %>% mutate(ord=1:29)
HS <- function(x,y){
  x <- x %>% 
    filter(p.value < 0.05) %>%
    mutate(min.rho = min(rho), max.rho = max(rho))
  
  dom <- rownames(y) 
  y <- as_tibble(y)
  y$domain <- dom
  
  z <- semi_join(y,x, by = 'domain')
  z.dom <- z$domain
  z <- as.matrix(z %>% select(-domain))
  row.names(z) <- z.dom
  
  ra = rowAnnotation(rho = x$rho,
                     gp = gpar(col = 'grey24'),
                     col = list(
                       rho = colorRamp2(c(-0.5, 0, 0.5), c("blue", "white", "red")),
                       show_legend = FALSE
                       ))
                       
  
  x <- x %>% arrange(rho)
  Heatmap(z, bottom_annotation = ha, name = 'z-score', border = TRUE, row_order = rev(x$domain),
          right_annotation = ra, show_heatmap_legend = TRUE, column_order = md$isolate,
          col = rev(magma(100)))
}

HS(srz.cor, srz)
HS(pfz.cor, pfz)
HS(trz.cor, trz)

HS(srp.cor, srp)

HS(srm.cor, srm)

