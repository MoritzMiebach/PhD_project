## This script analyses RNAseq data from featureCount output to functional enrichment analysis
## DEG analysis is performed via edgeR

# load the edgeR library
library(edgeR)
library(statmod)
library(tidyverse)
library(ComplexHeatmap) # ???
library(RColorBrewer)
library(circlize)
library(gplots)
library(S4Vectors)
library(org.At.tair.db)
library(org.Mm.eg.db)
library(UpSetR)
library(ggpubr)
library(magick)
library("ggVennDiagram")
library(ggdendro)
library(dendextend)

#global heatmap setting
ht_opt(legend_border = "black",
       heatmap_border = TRUE,
       annotation_border = TRUE)

# Data import
setwd("~/My Drive/Transcriptome/RNAseq_AC/v2")
df <- read.delim('counts_a_AC.txt', row.names ='Geneid')
head(df)

## removing genes that were consistently expressed at a given timepoint in the same direction in all three samples
ko.at <- read_csv('ko.at.overlap>1.csv')
dft <- as_tibble(df)
dft$AGI <- rownames(df)

df <- as.data.frame(anti_join(dft,ko.at,by='AGI'))
rownames(df) <- df[,64]
df <- df[,1:63]
head(df)


# Define sample groups
group <- factor(rep(c('MOCKa', 'AC229', 'AC266', 'AC464', 'AC754', 'AC759', 'AC761',
                      'MOCKb', 'AC103', 'AC248', 'AC433', 'AC435', 'AC561', 'AC854',
                      'MOCKc', 'AC221', 'AC286', 'AC463', 'AC507', 'AC832', 'AC861'), each = 3))
table(group)

y <- DGEList(group = group, counts = df)
y$samples
head(y$counts)

# Filtering to remove low counts
# keep genes that have CPM values above 0.5 in at least 3 libraries (each group contains 3 replicates)
keep <- rowSums(cpm(y) > 0.5) >= 3
table(keep)

# subset the DGEList to only retain the non-filtered genes
y <- y[keep, ,keep.lib.sizes=FALSE]

# Normalisation for composition bias
y <- calcNormFactors(y)
y$samples

# data exploration
plotMDS(y, col=rep(1:7, each=3))

plotMDS(y, bg=rep(c("grey24", "white", "white", "white", "white", "white", "white"), each = 3),
        pch = rep(c(21,21,21,22,22,22,23,23,23), each=7), cex=2, bty="n")
box(lwd=2.5)
legend('bottomright', bg="transparent", legend=c('mock','ACisolates'), 
       pt.bg=c("grey24", "white"),
       pch = 21, box.lty = 0, cex=1.1) 

plotMDS(y, bg=alpha(rep(c("grey24", "#bed1b0", "#ffe998", "#97d8ff", "#97d8ff", "#fdc3bf", "#ffe998",
                          "grey24", "#bed1b0", "#bed1b0", "#98e4bb", "#98e4bb", "#98e4bb", "#bed1b0",
                          "grey24", "#bed1b0", "#bed1b0", "#bed1b0", "#98e4bb", "#98e4bb", "#fdc3bf"), each=3),0.7),
        pch = rep(c(21,21,21,22,22,22,24,24,24), each=7), cex=2, bty="n", alpha=.5)
box(lwd=2.5)
legend('bottomright', bg="transparent", legend=c('mock','Alphaproteobacteria', 'Gammaproteobacteria', 'Firmicutes', 'Actinobacteria', 'Bacteriodetes'), 
       pt.bg=c("grey24", "#98e4bb", "#bed1b0", "#ffe998", "#fdc3bf", "#97d8ff"),
       pch = 21, box.lty = 0, cex=1.1) 

plotMDS(y, bg=alpha(rep(c("#7F7F7F", "#F5C242", "#F8D374", "#FFFFF1", "#F8D376", "#F8D880", "#FAE4A4",
                          "#7F7F7F", "#F5C242", "#FAE8B0", "#F7CB5D", "#F7CF69", "#FBE9B1", "#FAE4A5",
                          "#7F7F7F", "#FADF98", "#F7CF69", "#F8D880", "#F8D373", "#F8D579", "#F7D06C"), each=3),1),
        pch = rep(c(21,21,21,22,22,22,24,24,24), each=7), cex=2, bty="n")
box(lwd=2.5)
#exported as 6x10 inch pdf to reflect 15% vs 25% (4x1.5 & 4x2.5)

# Expression profiles of individual samples (Mean-Difference plots)
par(mfcol=c(3,4))
plotMD(y, column=6, ylim=c(-4,6)) #change column from 1 to 12
box(lwd=2.5)
abline(h=0, col="red", lwd=2)

dev.off()
# design matrix
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
design

# dispersion estimation
y <- estimateDisp(y, design, robust=TRUE)
plotBCV(y)

fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)

plotQLDisp(fit)

summary(fit$df.prior)

# Differential expression relative to a fold-change threshold
# Determine optimal fold-change threshold by rarefaction curve
Mv229 <- makeContrasts(AC229 - MOCKa, levels=design)
Mv266 <- makeContrasts(AC266 - MOCKa, levels=design)
Mv464 <- makeContrasts(AC464 - MOCKa, levels=design)
Mv754 <- makeContrasts(AC754 - MOCKa, levels=design)
Mv759 <- makeContrasts(AC759 - MOCKa, levels=design)
Mv761 <- makeContrasts(AC761 - MOCKa, levels=design)

Mv103 <- makeContrasts(AC103 - MOCKb, levels=design)
Mv248 <- makeContrasts(AC248 - MOCKb, levels=design)
Mv433 <- makeContrasts(AC433 - MOCKb, levels=design)
Mv435 <- makeContrasts(AC435 - MOCKb, levels=design)
Mv561 <- makeContrasts(AC561 - MOCKb, levels=design)
Mv854 <- makeContrasts(AC854 - MOCKb, levels=design)

Mv221 <- makeContrasts(AC221 - MOCKc, levels=design)
Mv286 <- makeContrasts(AC286 - MOCKc, levels=design)
Mv463 <- makeContrasts(AC463 - MOCKc, levels=design)
Mv507 <- makeContrasts(AC507 - MOCKc, levels=design)
Mv832 <- makeContrasts(AC832 - MOCKc, levels=design)
Mv861 <- makeContrasts(AC861 - MOCKc, levels=design)

## Differential expression relative to a fold-change threshold
DE.FC <- function(x){
  #y <- glmTreat(fit, contrast=x, lfc=log2(1))
  y <- glmQLFTest(fit, contrast=x)
  topTags(y)
  is.de <- decideTestsDGE(y)
  return(is.de)
  plotMD(y, status=is.de, values=c(1,-1), col=c("red","blue"), legend='none',ylim=c(-8,6))
  box(lwd=2.5)
  return(is.de)
}

is.de103 <- DE.FC(Mv103)
is.de221 <- DE.FC(Mv221)
is.de229 <- DE.FC(Mv229)
is.de248 <- DE.FC(Mv248)
is.de266 <- DE.FC(Mv266)
is.de286 <- DE.FC(Mv286)
is.de433 <- DE.FC(Mv433)
is.de435 <- DE.FC(Mv435)
is.de463 <- DE.FC(Mv463)
is.de464 <- DE.FC(Mv464)
is.de507 <- DE.FC(Mv507)
is.de561 <- DE.FC(Mv561)
is.de754 <- DE.FC(Mv754)
is.de759 <- DE.FC(Mv759)
is.de761 <- DE.FC(Mv761)
is.de832 <- DE.FC(Mv832)
is.de854 <- DE.FC(Mv854)
is.de861 <- DE.FC(Mv861)

## UpSet plot
DE <- function(x,y){
  z <- as.data.frame(x)
  z$ID <- rownames(z)
  colnames(z) <- c(y,'ID')
  z <- as_tibble(z)
  return(z)
}
de103 <- DE(is.de103, "AC103.ud")
de221 <- DE(is.de221, "AC221.ud")
de229 <- DE(is.de229, "AC229.ud")
de248 <- DE(is.de248, "AC248.ud")
de266 <- DE(is.de266, "AC266.ud")
de286 <- DE(is.de286, "AC286.ud")
de433 <- DE(is.de433, "AC433.ud")
de435 <- DE(is.de435, "AC435.ud")
de463 <- DE(is.de463, "AC463.ud")
de464 <- DE(is.de464, "AC464.ud")
de507 <- DE(is.de507, "AC507.ud")
de561 <- DE(is.de561, "AC561.ud")
de754 <- DE(is.de754, "AC754.ud")
de759 <- DE(is.de759, "AC759.ud")
de761 <- DE(is.de761, "AC761.ud")
de832 <- DE(is.de832, "AC832.ud")
de854 <- DE(is.de854, "AC854.ud")
de861 <- DE(is.de861, "AC861.ud")

de <- merge(merge(de854,de861, by = 'ID'),merge(merge(merge(merge(de103,de221, by = 'ID'),merge(de229,de248, by = 'ID'), by = 'ID'),
                                                      merge(merge(de266,de286, by = 'ID'),merge(de433,de435, by = 'ID'), by = 'ID'), by = 'ID'),
                                                merge(merge(merge(de463,de464, by = 'ID'),merge(de507,de561, by = 'ID'), by = 'ID'),
                                                      merge(merge(de754,de759, by = 'ID'),merge(de761,de832, by = 'ID'), by = 'ID'), by = 'ID'), by = 'ID'), by = 'ID')


head(de)
de <- de %>% group_by(ID) %>% mutate(Vario103 = AC103.ud^2,Pseudo221 = AC221.ud^2,Pseudo229 = AC229.ud^2,Steno248 = AC248.ud^2,Bacil266 = AC266.ud^2,Pseudo286 = AC286.ud^2,
                                     Meth433 = AC433.ud^2,Sphingo435 = AC435.ud^2,Pseudo463 = AC463.ud^2,Pedo464 = AC464.ud^2,Brevu507 = AC507.ud^2,Belna561 = AC561.ud^2,
                                     Chrys754 = AC754.ud^2,Frigo759 = AC759.ud^2,Staph761 = AC761.ud^2,Meth832 = AC832.ud^2,Steno854 = AC854.ud^2,Rhodo861 = AC861.ud^2)
summary(de)
de <- de %>% group_by(ID) %>% mutate(sum.de = sum(c(Vario103,Pseudo221,Pseudo229,Steno248,Bacil266,Pseudo286,Meth433,Sphingo435,Pseudo463,
                                                    Pedo464,Brevu507,Belna561,Chrys754,Frigo759,Staph761,Meth832,Steno854,Rhodo861)),
                                     sum.ud = sum(c(AC103.ud,AC221.ud,AC229.ud,AC248.ud,AC266.ud,AC286.ud,AC433.ud,AC435.ud,AC463.ud,
                                                    AC464.ud,AC507.ud,AC561.ud,AC754.ud,AC759.ud,AC761.ud,AC832.ud,AC854.ud,AC861.ud))) %>%
  filter(sum.de > 0)
de <- de %>% group_by(ID) %>% mutate(updown = ifelse(sum.de == sum.ud, 'up', 
                                                     ifelse(sum.de^2 ==sum.ud^2, 'down', 'up_down')))
## !!!! get up and down data. below is only for upregulated as I filtered negatives away
de_upset <- as.data.frame(dplyr::select(de, -c(sum.ud, sum.de, AC103.ud,AC221.ud,AC229.ud,AC248.ud,AC266.ud,AC286.ud,AC433.ud,AC435.ud,AC463.ud,
                                               AC464.ud,AC507.ud,AC561.ud,AC754.ud,AC759.ud,AC761.ud,AC832.ud,AC854.ud,AC861.ud)))
head(de_upset)
summary(de_upset)
upset(de_upset, order.by = "freq", nsets = 18, #nintersects = 200,
      text.scale = 1.5, point.size = 6, line.size = 1.5, mb.ratio = c(0.6, 0.4),
      mainbar.y.label = "DEGenes", sets.x.label = "DEGenes",
      sets.bar.color=rev(c("#bed1b0","#bed1b0","#bed1b0","#ffe998","#bed1b0","#bed1b0",
                           "#bed1b0","#97d8ff","#98e4bb","#98e4bb","#97d8ff","#fdc3bf",
                           "#98e4bb","#ffe998","#fdc3bf","#98e4bb","#98e4bb","#bed1b0"
      )),
      #group.by = 'sets',cutoff = 10,
      queries = list(
        list(
          query = elements, 
          params = list("updown", 'up'),
          color = "grey54", 
          active = T
        )
      )
)
# sorted by phylogeny
upset(de_upset, order.by = "freq", nsets = 18, #nintersects = 200,
      text.scale = 1.5, point.size = 6, line.size = 1.5, mb.ratio = c(0.6, 0.4),
      mainbar.y.label = "DEGenes", sets.x.label = "DEGenes",
      sets = c('Pseudo221','Vario103','Pseudo229','Pseudo286','Steno248','Pseudo463','Steno854',
               'Brevu507','Meth433','Meth832','Sphingo435','Belna561',
               'Pedo464','Chrys754','Bacil266','Staph761','Frigo759','Rhodo861'),
      keep.order = TRUE,
      sets.bar.color=rev(c("#fdc3bf","#fdc3bf","#ffe998","#ffe998","#97d8ff","#97d8ff",
                           "#98e4bb","#98e4bb","#98e4bb","#98e4bb","#98e4bb",
                           "#bed1b0","#bed1b0","#bed1b0","#bed1b0","#bed1b0","#bed1b0","#bed1b0"
      )),
      # group.by = 'sets',cutoff = 2,
      queries = list(
        list(
          query = elements, 
          params = list("updown", 'up'),
          color = "grey54", 
          active = T
        )
      )
)
# sorted by experiment
upset(de_upset, order.by = "freq", nsets = 18, #nintersects = 200,
      text.scale = 1.5, point.size = 6, line.size = 1.5, mb.ratio = c(0.6, 0.4),
      mainbar.y.label = "DEGenes", sets.x.label = "DEGenes",
      sets = c('Pseudo229','Bacil266','Pedo464','Chrys754','Frigo759','Staph761',
               'Vario103','Steno248','Meth433','Sphingo435','Belna561','Steno854',
               'Pseudo221','Pseudo286','Pseudo463','Brevu507','Meth832','Rhodo861'),
      keep.order = TRUE,
      sets.bar.color=rev(c("#214d65","#214d65","#214d65","#214d65","#214d65","#214d65",
                           "#287DAB","#287DAB","#287DAB","#287DAB","#287DAB","#287DAB",
                           "#E5BF86","#E5BF86","#E5BF86","#E5BF86","#E5BF86","#E5BF86"
      )),
      #group.by = 'sets',cutoff = 2,
      queries = list(
        list(
          query = elements, 
          params = list("updown", 'up'),
          color = "grey54", 
          active = T
        )
      )
)
# exported as 10x12 inch pdf


## Heatmap clustering

logCPM <- cpm(y, prior.count=2, log=TRUE)  
#prior.count = average count to be added to each observation to avoid taking log of zero
#prior.count reduces the variability of logCPM values of low count genes.

#logCPM <- cpm(y, log=TRUE)

logCPMgeneid <- rownames(logCPM)
logCPM <- as_tibble(logCPM)
logCPM$ID <- logCPMgeneid
logCPM <- right_join(logCPM, de, by='ID')
logCPM <- logCPM %>% dplyr::select('ID', MOCKa_1,MOCKa_2,MOCKa_3,MOCKb_1,MOCKb_2,MOCKb_3,MOCKc_1,MOCKc_2,MOCKc_3,
                                   AC103_1,AC103_2,AC103_3,AC221_1,AC221_2,AC221_3,AC229_1,AC229_2,AC229_3,
                                   AC248_1,AC248_2,AC248_3,AC266_1,AC266_2,AC266_3,AC286_1,AC286_2,AC286_3,
                                   AC433_1,AC433_2,AC433_3,AC435_1,AC435_2,AC435_3,AC463_1,AC463_2,AC463_3,
                                   AC464_1,AC464_2,AC464_3,AC507_1,AC507_2,AC507_3,AC561_1,AC561_2,AC561_3,
                                   AC754_1,AC754_2,AC754_3,AC759_1,AC759_2,AC759_3,AC761_1,AC761_2,AC761_3,
                                   AC832_1,AC832_2,AC832_3,AC854_1,AC854_2,AC854_3,AC861_1,AC861_2,AC861_3)
logCPM <- as.data.frame(logCPM)
logCPMgeneid <- logCPM[,1]
logCPM <- logCPM[,2:64]
row.names(logCPM) <- logCPMgeneid
#scale each row/gene to have mean zero and standard deviation one
#logCPM <- t(scale(t(logCPM)))
#logCPM <- as_tibble(logCPM)
#logCPM <- logCPM %>% rowwise() %>% summarise(Mock1 = MOCK_1, Mock2 = MOCK_2, Mock3 = MOCK_3,
# Pst1 = PST_1, Pst2 = PST_2, Pst3 = PST_3,
# Willi1 = WILLI_1, Willi2 = WILLI_2, Willi3 = WILLI_3,
# Micro1 = MICRO_1, Micro2 = MICRO_2, Micro3 = MICRO_3)
#or
#global z-score
#class(logCPM)
#m <- data.frame()
#sde <- data.frame()
#for(i in 1:12){
#m[i,1] <- mean(logCPM[,i], na.rm = TRUE)
#sde[i,1] <- sd(logCPM[,i], na.rm = TRUE)}
#m <- mean(m[,1])
#sde <- mean(sde[,1])
#logCPM <- (logCPM-m)/sde
#logCPM <- as_tibble(logCPM)
#logCPM <- logCPM %>% rowwise() %>% summarise(Mock1 = MOCK_1, Mock2 = MOCK_2, Mock3 = MOCK_3,
# Pst1 = PST_1, Pst2 = PST_2, Pst3 = PST_3,
# Willi1 = WILLI_1, Willi2 = WILLI_2, Willi3 = WILLI_3,
# Micro1 = MICRO_1, Micro2 = MICRO_2, Micro3 = MICRO_3)
#or
#scale the rows by subtracting the row mean or median (see Vogel2016)

logCPM <- as_tibble(logCPM)
logCPM <- logCPM %>% rowwise() %>% mutate(m = mean(c(MOCKa_1,MOCKa_2,MOCKa_3,MOCKb_1,MOCKb_2,MOCKb_3,MOCKc_1,MOCKc_2,MOCKc_3,
                                                     AC103_1,AC103_2,AC103_3,AC221_1,AC221_2,AC221_3,AC229_1,AC229_2,AC229_3,
                                                     AC248_1,AC248_2,AC248_3,AC266_1,AC266_2,AC266_3,AC286_1,AC286_2,AC286_3,
                                                     AC433_1,AC433_2,AC433_3,AC435_1,AC435_2,AC435_3,AC463_1,AC463_2,AC463_3,
                                                     AC464_1,AC464_2,AC464_3,AC507_1,AC507_2,AC507_3,AC561_1,AC561_2,AC561_3,
                                                     AC754_1,AC754_2,AC754_3,AC759_1,AC759_2,AC759_3,AC761_1,AC761_2,AC761_3,
                                                     AC832_1,AC832_2,AC832_3,AC854_1,AC854_2,AC854_3,AC861_1,AC861_2,AC861_3)))
logCPM <- logCPM %>% rowwise() %>% summarise(MOCKa_1 = MOCKa_1 - m, MOCKa_2 = MOCKa_2 - m, MOCKa_3 = MOCKa_3 -m,
                                             MOCKb_1 = MOCKb_1 - m, MOCKb_2 = MOCKb_2 - m, MOCKb_3 = MOCKb_3 -m,
                                             MOCKc_1 = MOCKc_1 - m, MOCKc_2 = MOCKc_2 - m, MOCKc_3 = MOCKc_3 -m,
                                             AC103_1 = AC103_1 - m, AC103_2 = AC103_2 - m, AC103_3 = AC103_3 -m,
                                             AC221_1 = AC221_1 - m, AC221_2 = AC221_2 - m, AC221_3 = AC221_3 -m,
                                             AC229_1 = AC229_1 - m, AC229_2 = AC229_2 - m, AC229_3 = AC229_3 -m,
                                             AC248_1 = AC248_1 - m, AC248_2 = AC248_2 - m, AC248_3 = AC248_3 -m,
                                             AC266_1 = AC266_1 - m, AC266_2 = AC266_2 - m, AC266_3 = AC266_3 -m,
                                             AC286_1 = AC286_1 - m, AC286_2 = AC286_2 - m, AC286_3 = AC286_3 -m,
                                             AC433_1 = AC433_1 - m, AC433_2 = AC433_2 - m, AC433_3 = AC433_3 -m,
                                             AC435_1 = AC435_1 - m, AC435_2 = AC435_2 - m, AC435_3 = AC435_3 -m,
                                             AC463_1 = AC463_1 - m, AC463_2 = AC463_2 - m, AC463_3 = AC463_3 -m,
                                             AC464_1 = AC464_1 - m, AC464_2 = AC464_2 - m, AC464_3 = AC464_3 -m,
                                             AC507_1 = AC507_1 - m, AC507_2 = AC507_2 - m, AC507_3 = AC507_3 -m,
                                             AC561_1 = AC561_1 - m, AC561_2 = AC561_2 - m, AC561_3 = AC561_3 -m,
                                             AC754_1 = AC754_1 - m, AC754_2 = AC754_2 - m, AC754_3 = AC754_3 -m,
                                             AC759_1 = AC759_1 - m, AC759_2 = AC759_2 - m, AC759_3 = AC759_3 -m,
                                             AC761_1 = AC761_1 - m, AC761_2 = AC761_2 - m, AC761_3 = AC761_3 -m,
                                             AC832_1 = AC832_1 - m, AC832_2 = AC832_2 - m, AC832_3 = AC832_3 -m,
                                             AC854_1 = AC854_1 - m, AC854_2 = AC854_2 - m, AC854_3 = AC854_3 -m,
                                             AC861_1 = AC861_1 - m, AC861_2 = AC861_2 - m, AC861_3 = AC861_3 -m)


## since data from 3 experiments was combined the gene expression changes need to be normalised to the mock treatment 
e1 <- as_tibble(y$counts) %>% dplyr::select(c('MOCKa_1','MOCKa_2','MOCKa_3','AC229_1','AC229_2','AC229_3',
                                              'AC266_1','AC266_2','AC266_3','AC464_1','AC464_2','AC464_3',
                                              'AC754_1','AC754_2','AC754_3','AC759_1','AC759_2','AC759_3',
                                              'AC761_1','AC761_2','AC761_3'))
e1 <- e1+2 # add a priot count. otherwise if a gene has 0 counts for mock every FC will be Inf.
e1$ID <- rownames(y$count)
e1 <- semi_join(e1, de, by='ID')
e1 <- e1 %>% pivot_longer(1:21,names_to='sample',values_to='n') 
e1m <- e1 %>% filter(sample %in% c('MOCKa_1','MOCKa_2','MOCKa_3')) %>% group_by(ID) %>%
  summarise(mean.mock=mean(n))
e1 <- left_join(e1,e1m,by='ID') %>% mutate(norm.n=log2(n/mean.mock)) %>% dplyr::select(ID,sample,norm.n) %>% 
  pivot_wider(names_from='sample',values_from='norm.n') 

e2 <- as_tibble(y$counts) %>% dplyr::select(c('MOCKb_1','MOCKb_2','MOCKb_3','AC103_1','AC103_2','AC103_3',
                                              'AC248_1','AC248_2','AC248_3','AC433_1','AC433_2','AC433_3',
                                              'AC435_1','AC435_2','AC435_3','AC561_1','AC561_2','AC561_3',
                                              'AC854_1','AC854_2','AC854_3'))
e2 <- e2+2 # add a priot count. otherwise if a gene has 0 counts for mock every FC will be Inf.
e2$ID <- rownames(y$count)
e2 <- semi_join(e2, de, by='ID')
e2 <- e2 %>% pivot_longer(1:21,names_to='sample',values_to='n') 
e2m <- e2 %>% filter(sample %in% c('MOCKb_1','MOCKb_2','MOCKb_3')) %>% group_by(ID) %>%
  summarise(mean.mock=mean(n))
e2 <- left_join(e2,e2m,by='ID') %>% mutate(norm.n=log2(n/mean.mock)) %>% dplyr::select(ID,sample,norm.n) %>% 
  pivot_wider(names_from='sample',values_from='norm.n') 

e3 <- as_tibble(y$counts) %>% dplyr::select(c('MOCKc_1','MOCKc_2','MOCKc_3','AC221_1','AC221_2','AC221_3',
                                              'AC286_1','AC286_2','AC286_3','AC463_1','AC463_2','AC463_3',
                                              'AC507_1','AC507_2','AC507_3','AC832_1','AC832_2','AC832_3',
                                              'AC861_1','AC861_2','AC861_3'))
e3 <- e3+2 # add a priot count. otherwise if a gene has 0 counts for mock every FC will be Inf.
e3$ID <- rownames(y$count)
e3 <- semi_join(e3, de, by='ID')
e3 <- e3 %>% pivot_longer(1:21,names_to='sample',values_to='n') 
e3m <- e3 %>% filter(sample %in% c('MOCKc_1','MOCKc_2','MOCKc_3')) %>% group_by(ID) %>%
  summarise(mean.mock=mean(n))
e3 <- left_join(e3,e3m,by='ID') %>% mutate(norm.n=log2(n/mean.mock)) %>% dplyr::select(ID,sample,norm.n) %>% 
  pivot_wider(names_from='sample',values_from='norm.n')  

ea <- left_join(left_join(e1,e2,by='ID'),e3,by='ID')


####
an <- HeatmapAnnotation(taxa = ea$ID,
                        prot.resp = ea$ID,
                        show_legend = FALSE,
                        col = list(taxa = c('MOCKa_1'='#7F7F7F','MOCKa_2'='#7F7F7F','MOCKa_3'='#7F7F7F',
                                            'AC229_1'='#bed1b0','AC229_2'='#bed1b0','AC229_3'='#bed1b0',
                                            'AC266_1'='#ffe998','AC266_2'='#ffe998','AC266_3'='#ffe998',
                                            'AC464_1'='#97d8ff','AC464_2'='#97d8ff','AC464_3'='#97d8ff',
                                            'AC754_1'='#97d8ff','AC754_2'='#97d8ff','AC754_3'='#97d8ff',
                                            'AC759_1'='#fdc3bf','AC759_2'='#fdc3bf','AC759_3'='#fdc3bf',
                                            'AC761_1'='#ffe998','AC761_2'='#ffe998','AC761_3'='#ffe998',
                                            'MOCKb_1'='#7F7F7F','MOCKb_2'='#7F7F7F','MOCKb_3'='#7F7F7F',
                                            'AC103_1'='#bed1b0','AC103_2'='#bed1b0','AC103_3'='#bed1b0',
                                            'AC248_1'='#bed1b0','AC248_2'='#bed1b0','AC248_3'='#bed1b0',
                                            'AC433_1'='#98e4bb','AC433_2'='#98e4bb','AC433_3'='#98e4bb',
                                            'AC435_1'='#98e4bb','AC435_2'='#98e4bb','AC435_3'='#98e4bb',
                                            'AC561_1'='#98e4bb','AC561_2'='#98e4bb','AC561_3'='#98e4bb',
                                            'AC854_1'='#bed1b0','AC854_2'='#bed1b0','AC854_3'='#bed1b0',
                                            'MOCKc_1'='#7F7F7F','MOCKc_2'='#7F7F7F','MOCKc_3'='#7F7F7F',
                                            'AC221_1'='#bed1b0','AC221_2'='#bed1b0','AC221_3'='#bed1b0',
                                            'AC286_1'='#bed1b0','AC286_2'='#bed1b0','AC286_3'='#bed1b0',
                                            'AC463_1'='#bed1b0','AC463_2'='#bed1b0','AC463_3'='#bed1b0',
                                            'AC507_1'='#98e4bb','AC507_2'='#98e4bb','AC507_3'='#98e4bb',
                                            'AC832_1'='#98e4bb','AC832_2'='#98e4bb','AC832_3'='#98e4bb',
                                            'AC861_1'='#fdc3bf','AC861_2'='#fdc3bf','AC861_3'='#fdc3bf'),
                                   prot.resp = c('MOCKa_1'="#7F7F7F",'MOCKa_2'="#7F7F7F",'MOCKa_3'="#7F7F7F",
                                                 'AC229_1'="#F5C242",'AC229_2'="#F5C242",'AC229_3'="#F5C242",
                                                 'AC266_1'="#F8D374",'AC266_2'="#F8D374",'AC266_3'="#F8D374",
                                                 'AC464_1'="#FFFFF1",'AC464_2'="#FFFFF1",'AC464_3'="#FFFFF1",
                                                 'AC754_1'="#F8D376",'AC754_2'="#F8D376",'AC754_3'="#F8D376",
                                                 'AC759_1'="#F8D880",'AC759_2'="#F8D880",'AC759_3'="#F8D880",
                                                 'AC761_1'="#FAE4A4",'AC761_2'="#FAE4A4",'AC761_3'="#FAE4A4",
                                                 'MOCKb_1'="#7F7F7F",'MOCKb_2'="#7F7F7F",'MOCKb_3'="#7F7F7F",
                                                 'AC103_1'="#F5C242",'AC103_2'="#F5C242",'AC103_3'="#F5C242",
                                                 'AC248_1'="#FAE8B0",'AC248_2'="#FAE8B0",'AC248_3'="#FAE8B0",
                                                 'AC433_1'="#F7CB5D",'AC433_2'="#F7CB5D",'AC433_3'="#F7CB5D",
                                                 'AC435_1'="#F7CF69",'AC435_2'="#F7CF69",'AC435_3'="#F7CF69",
                                                 'AC561_1'="#FBE9B1",'AC561_2'="#FBE9B1",'AC561_3'="#FBE9B1",
                                                 'AC854_1'="#FAE4A5",'AC854_2'="#FAE4A5",'AC854_3'="#FAE4A5",
                                                 'MOCKc_1'="#7F7F7F",'MOCKc_2'="#7F7F7F",'MOCKc_3'="#7F7F7F",
                                                 'AC221_1'="#FADF98",'AC221_2'="#FADF98",'AC221_3'="#FADF98",
                                                 'AC286_1'="#F7CF69",'AC286_2'="#F7CF69",'AC286_3'="#F7CF69",
                                                 'AC463_1'="#F8D880",'AC463_2'="#F8D880",'AC463_3'="#F8D880",
                                                 'AC507_1'="#F8D373",'AC507_2'="#F8D373",'AC507_3'="#F8D373",
                                                 'AC832_1'="#F8D579",'AC832_2'="#F8D579",'AC832_3'="#F8D579",
                                                 'AC861_1'="#F7D06C",'AC861_2'="#F7D06C",'AC861_3'="#F7D06C")))

# kmeans clustering outside of complexHeatmap function
# determine optimal number of gene clusters via the elbow method
ea.id <- ea$ID
ea <- as.data.frame(ea[,2:64])
rownames(ea) <- ea.id

kmean_withinss <- function(k) {
  cluster <- kmeans(ea, k, nstart=50, iter.max=100)
  return (cluster$tot.withinss)
}

wss <- sapply(2:40, kmean_withinss)
elbow <-data.frame(2:40, wss)
ggplot(elbow, aes(x = X2.40, y = wss))+
  geom_line()+
  geom_point(shape=21, fill = 'white', size = 2)+
  #scale_x_continuous(breaks = seq(1, 20, by = 1))+
  theme_bw()+
  theme(legend.position = 'right',
        aspect.ratio = 1.25,
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        panel.border = element_rect(size = 2, colour = 'black'),
        axis.title = element_text(colour = "black", size = 15),
        axis.ticks = element_blank(),
        # legend at the bottom 
        legend.text = element_text(colour = "black", size = 15),
        # plot title
        plot.title = element_text(colour = "black", size = 18),
        #axis.text.x = element_text(angle = 90)
  )+
  labs(x='# of clusters', y='total within sum of squares')

kclus.row <- kmeans(ea, 10, nstart=50, iter.max=100)
kclus.row$cluster
split.row <- factor(paste0("\n", kclus.row$cluster), 
                    levels=c("\n2","\n1","\n4","\n10","\n3",
                             "\n8","\n5","\n7","\n9", "\n6"))

ht2 <- Heatmap(ea, #cluster_column_slices = FALSE, column_split=split.col,
               row_split=split.row, cluster_row_slices = FALSE, show_row_names=FALSE,
               name = "logcpm", show_column_names=FALSE#, top_annotation = an
)

#Heatmap(ea, show_row_names=FALSE,rect_gp = gpar(type = "none"))

ht2 <- draw(ht2)
# exported as 5x7 inch pdf
r.dend <- row_dend(ht2)  #Extract row dendrogram
rcl.list <- row_order(ht2)  #Extract clusters (output is a list)
lapply(rcl.list, function(x) length(x))  #check/confirm size clusters
rcl.list[2]
# loop to extract genes for each cluster.
for (i in 1:length(row_order(ht2))){
  if (i == 1) {
    clu <- t(t(row.names(ea[row_order(ht2)[[i]],])))
    out <- cbind(clu, paste("cluster", i, sep=""))
    colnames(out) <- c("GeneID", "Cluster")
  } else {
    clu <- t(t(row.names(ea[row_order(ht2)[[i]],])))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    out <- rbind(out, clu)
  }
}
#check
out 

#export
#write.table(out, file= "ko.o2.cl10.txt", sep="\t", quote=F, row.names=FALSE)

#load out
setwd("~/My Drive/Transcriptome/RNAseq_AC/v2")
out <- read.delim('ko.o2.cl10.txt')

out <- as_tibble(out)

lcpm.rn <- rownames(ea) 
lcpm <- as_tibble(ea)
lcpm$GeneID <- lcpm.rn
lcpm <- lcpm %>% pivot_longer(1:63,names_to = 'treat', values_to = 'logFC') %>%
  mutate(treatment = if_else(treat %in% c("MOCKa_1", "MOCKa_2", "MOCKa_3"), 'MOCKa',
                             if_else(treat %in% c("MOCKb_1", "MOCKb_2", "MOCKb_3"), 'MOCKb',
                                     if_else(treat %in% c("MOCKc_1", "MOCKc_2", "MOCKc_3"), 'MOCKc',
                                             if_else(treat %in% c("AC103_1", "AC103_2", "AC103_3"), 'AC103',
                                                     if_else(treat %in% c("AC221_1", "AC221_2", "AC221_3"), 'AC221',
                                                             if_else(treat %in% c("AC229_1", "AC229_2", "AC229_3"), 'AC229',
                                                                     if_else(treat %in% c("AC248_1", "AC248_2", "AC248_3"), 'AC248',
                                                                             if_else(treat %in% c("AC266_1", "AC266_2", "AC266_3"), 'AC266',
                                                                                     if_else(treat %in% c("AC286_1", "AC286_2", "AC286_3"), 'AC286',
                                                                                             if_else(treat %in% c("AC433_1", "AC433_2", "AC433_3"), 'AC433',
                                                                                                     if_else(treat %in% c("AC435_1", "AC435_2", "AC435_3"), 'AC435',
                                                                                                             if_else(treat %in% c("AC463_1", "AC463_2", "AC463_3"), 'AC463',
                                                                                                                     if_else(treat %in% c("AC464_1", "AC464_2", "AC464_3"), 'AC464',
                                                                                                                             if_else(treat %in% c("AC507_1", "AC507_2", "AC507_3"), 'AC507',
                                                                                                                                     if_else(treat %in% c("AC561_1", "AC561_2", "AC561_3"), 'AC561',
                                                                                                                                             if_else(treat %in% c("AC754_1", "AC754_2", "AC754_3"), 'AC754',
                                                                                                                                                     if_else(treat %in% c("AC759_1", "AC759_2", "AC759_3"), 'AC759',
                                                                                                                                                             if_else(treat %in% c("AC761_1", "AC761_2", "AC761_3"), 'AC761',
                                                                                                                                                                     if_else(treat %in% c("AC832_1", "AC832_2", "AC832_3"), 'AC832',
                                                                                                                                                                             if_else(treat %in% c("AC854_1", "AC854_2", "AC854_3"), 'AC854',
                                                                                                                                                                                     if_else(treat %in% c("AC861_1", "AC861_2", "AC861_3"), 'AC861',
                                                                                                                                                                                             "x"))))))))))))))))))))),
         replicate = rep(1:3, times = n()/3))  %>%
  group_by(treatment, GeneID) %>% mutate(mean.logFC = mean(logFC)) 

lcpm <- left_join(lcpm, out , by = "GeneID")

wgs <- read_csv('WGS_metadata2.csv')

scaleFUN <- function(x) sprintf("%.1f", x)

##
left_join(lcpm,wgs,by='treatment') %>% dplyr::select(-c(treat,logFC,replicate)) %>% distinct() %>%
  filter(Cluster %in% 'cluster1') %>% filter(treatment != 'MOCKa') %>%
  filter(treatment != 'MOCKb') %>% filter(treatment != 'MOCKc') %>%
  group_by(treatment) %>% mutate(m = mean(mean.logFC,na.rm=TRUE), sd = sd(mean.logFC),n=n()) %>%
  ggplot()+
  labs(y='',x='', title='cluster1')+
  #geom_line(aes(factor(treatment,levels = c('Mock','Micro347', 'Willi354', 'Pst')), mean.logCPM, group = GeneID), alpha = .3)+
  #geom_line(aes(factor(treatment,levels = c('Mock','Micro347', 'Willi354', 'Pst')), m, group = GeneID), col = 'red', size = 1)+
  #geom_point(aes(factor(treatment,levels = c('Mock','Micro347', 'Willi354', 'Pst')), mean.logCPM), shape=21, size=3, alpha=.3, fill = "grey24")+
  geom_violin(aes(reorder(treatment,m.log.auc), mean.logFC),
              size=.5, alpha=.5, scale='width',fill='#FFFDD0')+
  geom_pointrange(aes(reorder(treatment,m.log.auc) ,m,ymax =m+sd,ymin=m-sd), size=.5)+
  scale_fill_manual(values=c('#fdc3bf','#98e4bb','#ffe998','#97d8ff','#bed1b0','grey24'))+
  #scale_color_manual(values = c("#e1af00","grey24","#f21a00","#3b9ab2"))+
  theme_bw()+
  theme(legend.position = 'none',
        aspect.ratio = .5,
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        #panel.grid.minor.y = element_line(size = 2, colour = 'black'),
        legend.title = element_text(size = 15),
        axis.text = element_text(color = 'black',size = 15),
        panel.border = element_rect(size = 2, colour = 'black'),
        axis.title = element_text(colour = "black", size = 18),
        axis.ticks = element_blank(),
        # legend at the bottom 
        legend.text = element_text(colour = "black", size = 15),
        axis.text.x = element_text(angle = 90),
        plot.title = element_text(colour = "black", size = 25)
  )+
  #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  scale_y_continuous(labels=scaleFUN)
#labs(y = expression(centered~log[2](CPM)), x = '') 


cor <- left_join(lcpm,wgs,by='treatment') %>% filter(!is.na(Cluster)) %>% filter(!is.na(m.log.auc)) %>%
  dplyr::select(treatment,mean.logFC,m.log.auc,se.log.auc,Cluster,classification,experiment,epi.endo) %>% distinct() %>%
  group_by(treatment,Cluster) %>% 
  mutate(m.logFC=mean(mean.logFC),se.logFC=sd(mean.logFC)/sqrt(n()),n=n())

cor %>% filter(Cluster=='cluster1') %>% dplyr::select(-c(GeneID,mean.logFC)) %>% distinct() %>%
  ggplot()+
  labs(x='log2(fold change)',y='log2(normalised AUC)', title='cluster10')+
  geom_errorbarh(aes(y=m.log.auc,xmin=m.logFC-se.logFC,xmax=m.logFC+se.logFC),size=.5)+
  geom_pointrange(aes(m.logFC,m.log.auc,ymin=m.log.auc-se.log.auc,ymax=m.log.auc+se.log.auc),
                  shape=21,fill='white', size=.5)+
  geom_smooth(aes(x=m.logFC,y=m.log.auc), method='lm', col='black')+
  geom_point(aes(x=m.logFC,y=m.log.auc,fill=experiment),size=4,shape=21)+
  #scale_fill_manual(values=c('#fdc3bf','#98e4bb','#ffe998','#97d8ff','#bed1b0','grey24'))+
  scale_fill_manual(values=c('#4053d3','#ddb310','#b51d14'))+
  theme_bw()+
  theme(legend.position = 'right',
        aspect.ratio = 1,
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.border = element_rect(size = 1, colour = 'black'),
        axis.title = element_text(colour = "black", size = 15),
        plot.title = element_text(colour = "black", size = 15),
        axis.ticks = element_blank(),
        # legend at the bottom 
        legend.text = element_text(colour = "black", size = 12),
        # facet wrap properties
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black", face = 'italic'),
        strip.text.y = element_text(size = 12, color = "black", face = 'italic'))

cor.cl1 <- cor %>% filter(Cluster=='cluster1') %>% dplyr::select(-c(GeneID,mean.logFC)) %>% distinct()
summary(lm(m.logFC~m.log.auc,data=cor.cl1))
cor.cl2 <- cor %>% filter(Cluster=='cluster2') %>% dplyr::select(-c(GeneID,mean.logFC)) %>% distinct()
summary(lm(m.logFC~m.log.auc,data=cor.cl2))
cor.cl3 <- cor %>% filter(Cluster=='cluster3') %>% dplyr::select(-c(GeneID,mean.logFC)) %>% distinct()
summary(lm(m.logFC~m.log.auc,data=cor.cl3))
cor.cl4 <- cor %>% filter(Cluster=='cluster4') %>% dplyr::select(-c(GeneID,mean.logFC)) %>% distinct()
summary(lm(m.logFC~m.log.auc,data=cor.cl4))
cor.cl5 <- cor %>% filter(Cluster=='cluster5') %>% dplyr::select(-c(GeneID,mean.logFC)) %>% distinct()
summary(lm(m.logFC~m.log.auc,data=cor.cl5))
cor.cl6 <- cor %>% filter(Cluster=='cluster6') %>% dplyr::select(-c(GeneID,mean.logFC)) %>% distinct()
summary(lm(m.logFC~m.log.auc,data=cor.cl6))
cor.cl7 <- cor %>% filter(Cluster=='cluster7') %>% dplyr::select(-c(GeneID,mean.logFC)) %>% distinct()
summary(lm(m.logFC~m.log.auc,data=cor.cl7))
cor.cl8 <- cor %>% filter(Cluster=='cluster8') %>% dplyr::select(-c(GeneID,mean.logFC)) %>% distinct()
summary(lm(m.logFC~m.log.auc,data=cor.cl8))
cor.cl9 <- cor %>% filter(Cluster=='cluster9') %>% dplyr::select(-c(GeneID,mean.logFC)) %>% distinct()
summary(lm(m.logFC~m.log.auc,data=cor.cl9))
cor.cl10 <- cor %>% filter(Cluster=='cluster10') %>% dplyr::select(-c(GeneID,mean.logFC)) %>% distinct()
summary(lm(m.logFC~m.log.auc,data=cor.cl10))


any(is.na(ea))
tt <- t(as.data.frame(ea))
hc.ea <- hclust(dist(tt), method='ward.D2')
plot(hc.ea)

ggdendrogram(hc.ea, size = 2)

ea.m <- left_join(lcpm,wgs,by='treatment') 
ea.m <- ea.m %>% ungroup() %>% dplyr::select(GeneID,symbol,mean.logFC) %>%
  distinct() %>%
  pivot_wider(names_from=symbol,values_from=mean.logFC)

ea.m.rn <- ea.m$GeneID
ea.m <- as.data.frame(ea.m[,2:22])
ea.m <- t(ea.m)

hc.ea.m <- hclust(dist(ea.m), method='ward.D2')
plot(hc.ea.m)
hc.ea.m <- as.dendrogram(hc.ea.m)

par(mar=c(10,1,1,1))
hc.ea.m %>% 
  #set("branches_k_color", value = c("#b80058", "#5954d6", "#878500"), k = 3) %>%
  set("branches_lwd", 2) %>%
  plot(axes=FALSE,edgePar = list(lwd = 4))
colored_bars(colors = c('#cacaca','#878500','#878500','#878500','#878500','#878500','#878500',
                        '#cacaca','#ddb310','#ddb310','#ddb310','#ddb310','#ddb310','#ddb310',
                        '#cacaca','#b80058','#b80058','#b80058','#b80058','#b80058','#b80058'),
             dend = hc.ea.m, rowLabels = "experiment")
colored_bars(colors = c('black','#4053d3','#4053d3','#4053d3','#4053d3','#4053d3','#4053d3',
                        '#cacaca','#ddb310','#ddb310','#ddb310','#ddb310','#ddb310','#ddb310',
                        '#cacaca','#b51d14','#b51d14','#b51d14','#b51d14','#b51d14','#b51d14'),
             dend = hc.ea.m, rowLabels = "log2(normalised AUC)")


## hclust annotation in ggplot (colored_bars looks shit)
wgs %>% mutate(c=1) %>%
  ggplot()+
  geom_point(aes(factor(symbol,levels=c('Pseudo221','Pseudo286','Rhodo861','Pseudo463',
                                        'Brevu507','Meth832','Meth433','Sphingo435','Belna561','Steno854','Frigo759','Staph761',
                                        'Pedo464','Chrys754','Pseudo229','Bacil266','Mock2','Mock1','Mock3','Vario103','Steno248')),
                 c,fill=m.log.auc),shape=21,size=7)+
  scale_fill_gradient(name='protoplast response',low='#FFFFF0',high='#FFBF00')+
  theme_void()+
  theme(legend.position='none')
wgs %>% mutate(c=1) %>%
  ggplot()+
  geom_point(aes(factor(symbol,levels=c('Pseudo221','Pseudo286','Rhodo861','Pseudo463',
                                        'Brevu507','Meth832','Meth433','Sphingo435','Belna561','Steno854','Frigo759','Staph761',
                                        'Pedo464','Chrys754','Pseudo229','Bacil266','Mock2','Mock1','Mock3','Vario103','Steno248')),
                 c,fill=classification),shape=21,size=7)+
  scale_fill_manual(values=c('#fdc3bf','#98e4bb','#97d8ff','#ffe998','#bed1b0'), name='classification')+
  theme_void()+
  theme(legend.position='right')
wgs %>% mutate(c=1) %>%
  ggplot()+
  geom_point(aes(factor(symbol,levels=c('Pseudo221','Pseudo286','Rhodo861','Pseudo463',
                                        'Brevu507','Meth832','Meth433','Sphingo435','Belna561','Steno854','Frigo759','Staph761',
                                        'Pedo464','Chrys754','Pseudo229','Bacil266','Mock2','Mock1','Mock3','Vario103','Steno248')),
                 c,fill=experiment),shape=21,size=7,alpha=.8)+
  scale_fill_manual(values=c('#4053d3','#ddb310','#b51d14'), name='experiment')+
  theme_void()+
  theme(legend.position='right')

wgs %>% mutate(c=1) %>%
  ggplot()+
  geom_point(aes(factor(symbol,levels=c('Pseudo221','Pseudo286','Rhodo861','Pseudo463',
                                        'Brevu507','Meth832','Meth433','Sphingo435','Belna561','Steno854','Frigo759','Staph761',
                                        'Pedo464','Chrys754','Pseudo229','Bacil266','Mock2','Mock1','Mock3','Vario103','Steno248')),
                 c,fill=epi.endo),shape=21,size=7,alpha=.8)+
  scale_fill_manual(values=c('#74a089', '#9a8822'), name='experiment')+
  theme_void()+
  theme(legend.position='right')

### subset of genes belonging to GO terms related to defense, immunity, resistance, bacterium,
### salicylic acid, ethylene and jasmonic acid

def <- read_csv('def.go.csv') %>% dplyr::rename(GeneID=AGI)

dm <- semi_join(lcpm,def,by='GeneID')
dm <- left_join(dm,wgs,by='treatment') 
dm <- dm %>% ungroup() %>% dplyr::select(GeneID,symbol,mean.logFC) %>%
  distinct() %>%
  pivot_wider(names_from=symbol,values_from=mean.logFC)

dm.rn <- dm$GeneID
dm <- as.data.frame(dm[,2:22])
dm <- t(dm)

hc.dm <- hclust(dist(dm), method='ward.D2')
plot(hc.dm)
hc.dm <- as.dendrogram(hc.dm)

par(mar=c(10,1,1,1))
hc.dm %>% 
  #set("branches_k_color", value = c("#b80058", "#5954d6", "#878500"), k = 3) %>%
  set("branches_lwd", 2) %>%
  plot(axes=FALSE,edgePar = list(lwd = 4))
## hclust annotation in ggplot (colored_bars looks shit)
wgs %>% mutate(c=1) %>%
  ggplot()+
  geom_point(aes(factor(symbol,levels=c('Rhodo861','Pseudo463','Brevu507','Meth832','Pseudo229','Pseudo221','Pseudo286',
                                        'Steno854','Belna561','Meth433','Sphingo435','Bacil266','Frigo759','Staph761',
                                        'Pedo464','Chrys754','Mock2','Mock1','Mock3','Vario103','Steno248')),
                 c,fill=m.log.auc),shape=21,size=7)+
  scale_fill_gradient(name='protoplast response',low='#FFFFF0',high='#FFBF00')+
  theme_void()+
  theme(legend.position='none')
wgs %>% mutate(c=1) %>%
  ggplot()+
  geom_point(aes(factor(symbol,levels=c('Rhodo861','Pseudo463','Brevu507','Meth832','Pseudo229','Pseudo221','Pseudo286',
                                        'Steno854','Belna561','Meth433','Sphingo435','Bacil266','Frigo759','Staph761',
                                        'Pedo464','Chrys754','Mock2','Mock1','Mock3','Vario103','Steno248')),
                 c,fill=classification),shape=21,size=7)+
  scale_fill_manual(values=c('#fdc3bf','#98e4bb','#97d8ff','#ffe998','#bed1b0'), name='classification')+
  theme_void()+
  theme(legend.position='none')
wgs %>% mutate(c=1) %>%
  ggplot()+
  geom_point(aes(factor(symbol,levels=c('Rhodo861','Pseudo463','Brevu507','Meth832','Pseudo229','Pseudo221','Pseudo286',
                                        'Steno854','Belna561','Meth433','Sphingo435','Bacil266','Frigo759','Staph761',
                                        'Pedo464','Chrys754','Mock2','Mock1','Mock3','Vario103','Steno248')),
                 c,fill=experiment),shape=21,size=7,alpha=.8)+
  scale_fill_manual(values=c('#4053d3','#ddb310','#b51d14'), name='experiment')+
  theme_void()+
  theme(legend.position='none')

wgs %>% mutate(c=1) %>%
  ggplot()+
  geom_point(aes(factor(symbol,levels=c('Rhodo861','Pseudo463','Brevu507','Meth832','Pseudo229','Pseudo221','Pseudo286',
                                        'Steno854','Belna561','Meth433','Sphingo435','Bacil266','Frigo759','Staph761',
                                        'Pedo464','Chrys754','Mock2','Mock1','Mock3','Vario103','Steno248')),
                 c,fill=epi.endo),shape=21,size=7,alpha=.8)+
  scale_fill_manual(values=c('#74a089', '#9a8822'), name='experiment')+
  theme_void()+
  theme(legend.position='right')
### do SA related PR gene expression correlate with protoplast response?

sag <- read_csv('SA.PR.genes.csv') %>% dplyr::rename(GeneID=AGI)

sa <- semi_join(cor,sag,by='GeneID') %>% dplyr::select(-c(m.logFC,se.logFC))

sa %>% filter(GeneID %in% 'AT2G14610') %>%
  filter(experiment %in% 'exp3') %>%
  group_by(treatment) %>% 
  mutate(m.logFC=mean(mean.logFC),se.logFC=sd(mean.logFC)/sqrt(n()),n=n()) %>%
  dplyr::select(-c(GeneID,mean.logFC)) %>%
  distinct() %>%
  ggplot()+
  labs(x='log2(fold change)',y='log2(normalised AUC)')+
  geom_errorbarh(aes(y=m.log.auc,xmin=m.logFC-se.logFC,xmax=m.logFC+se.logFC),size=.5)+
  geom_pointrange(aes(m.logFC,m.log.auc,ymin=m.log.auc-se.log.auc,ymax=m.log.auc+se.log.auc),
                  shape=21,fill='white', size=.5)+
  geom_smooth(aes(x=m.logFC,y=m.log.auc), method='lm', col='black')+
  geom_point(aes(x=m.logFC,y=m.log.auc,fill=epi.endo),size=4,shape=21)+
  scale_fill_manual(values=c('#fdc3bf','#97d8ff','#ffe998','#bed1b0'))+ #exp1
  #scale_fill_manual(values=c('#98e4bb','#bed1b0'))+ #exp2
  #scale_fill_manual(values=c('#fdc3bf','#98e4bb','#bed1b0'))+ #exp3
  #scale_fill_manual(values=c('#4053d3','#ddb310','#b51d14'))+
  theme_bw()+
  theme(legend.position = 'right',
        aspect.ratio = 1,
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size = 14),
        axis.text = element_text(size = 12),
        panel.border = element_rect(size = 1, colour = 'black'),
        axis.title = element_text(colour = "black", size = 15),
        plot.title = element_text(colour = "black", size = 15),
        axis.ticks = element_blank(),
        # legend at the bottom 
        legend.text = element_text(colour = "black", size = 12),
        # facet wrap properties
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid"),
        strip.text.x = element_text(size = 12, color = "black", face = 'italic'),
        strip.text.y = element_text(size = 12, color = "black", face = 'italic'))+
  ylim(-1,.25)

sae1 <- sa %>% filter(GeneID %in% 'AT2G14610') %>%
  filter(experiment %in% 'exp1') %>%
  group_by(treatment) %>% 
  mutate(m.logFC=mean(mean.logFC),se.logFC=sd(mean.logFC)/sqrt(n()),n=n())
sae2 <- sa %>% filter(GeneID %in% 'AT2G14610') %>%
  filter(experiment %in% 'exp2') %>%
  group_by(treatment) %>% 
  mutate(m.logFC=mean(mean.logFC),se.logFC=sd(mean.logFC)/sqrt(n()),n=n())
sae3 <- sa %>% filter(GeneID %in% 'AT2G14610') %>%
  filter(experiment %in% 'exp3') %>%
  group_by(treatment) %>% 
  mutate(m.logFC=mean(mean.logFC),se.logFC=sd(mean.logFC)/sqrt(n()),n=n())
summary(lm(m.logFC~m.log.auc,data=sae1))
summary(lm(m.logFC~m.log.auc,data=sae2))
summary(lm(m.logFC~m.log.auc,data=sae3))


## number of DEGs
de_upset %>% filter(Vario103==1) %>% summarise(n=n())
de_upset %>% filter(Pseudo221==1) %>% summarise(n=n())
de_upset %>% filter(Pseudo229==1) %>% summarise(n=n())
de_upset %>% filter(Steno248==1) %>% summarise(n=n())
de_upset %>% filter(Bacil266==1) %>% summarise(n=n())
de_upset %>% filter(Pseudo286==1) %>% summarise(n=n())
de_upset %>% filter(Meth433==1) %>% summarise(n=n())
de_upset %>% filter(Sphingo435==1) %>% summarise(n=n())
de_upset %>% filter(Pseudo463==1) %>% summarise(n=n())
de_upset %>% filter(Pedo464==1) %>% summarise(n=n())
de_upset %>% filter(Brevu507==1) %>% summarise(n=n())
de_upset %>% filter(Belna561==1) %>% summarise(n=n())
de_upset %>% filter(Chrys754==1) %>% summarise(n=n())
de_upset %>% filter(Frigo759==1) %>% summarise(n=n())
de_upset %>% filter(Staph761==1) %>% summarise(n=n())
de_upset %>% filter(Meth832==1) %>% summarise(n=n())
de_upset %>% filter(Steno854==1) %>% summarise(n=n())
de_upset %>% filter(Rhodo861==1) %>% summarise(n=n())

