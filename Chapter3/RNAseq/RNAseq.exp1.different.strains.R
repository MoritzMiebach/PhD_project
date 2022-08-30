## This script analyses RNAseq data from featureCount output to functional enrichment analysis
## DEG analysis is performed via edgeR

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("S4Vectors")
BiocManager::install("ComplexHeatmap")
BiocManager::install("GO.db")
BiocManager::install("org.At.eg.db")
BiocManager::install("edgeR")
BiocManager::install("org.At.tair.db")
BiocManager::install("org.Mm.eg.db")


# load the edgeR library
library(edgeR)
library(statmod)
library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(gplots)
library(S4Vectors)
library(org.At.tair.db)
library(org.Mm.eg.db)
library(UpSetR)
library(ggpubr)
library("ggVennDiagram")

#global heatmap setting
ht_opt(legend_border = "black",
       heatmap_border = TRUE,
       annotation_border = TRUE)

# Data import
setwd("/Volumes/GoogleDrive/My Drive/Transcriptome/RNAseq_Bai_clean")
df <- read.delim('counts_a_Rinput.txt', row.names ='Geneid')
head(df)

# Define sample groups
group <- factor(rep(c('MOCK', 'PST', 'WILLI', 'MICRO'), each = 3))
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
plotMDS(y, col=rep(1:5, each=3))

plotMDS(y, bg=rep(c("grey24", "#f21a00", "#3b9ab2", "#e1af00"), each = 3), pch = 21, cex=2, bty="n")
box(lwd=2.5)
legend('topleft', bg="transparent", legend=c('mock','Pst','Willi354','Micro347'), 
       pt.bg=c("grey24", "#f21a00", "#3b9ab2", "#e1af00"),
       pch = 21, box.lty = 0, cex=1.1) 
#exported as 5x5 inch pdf

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
PvsM <- makeContrasts(PST - MOCK, levels=design)
MvsM <- makeContrasts(MICRO - MOCK, levels=design)
WvsM <- makeContrasts(WILLI - MOCK, levels=design)

oFC <- function(x, y){
  df <- data.frame()
  for(i in 1:200){
    j <- i/100+1
    gT <- glmTreat(fit, contrast=x, lfc=log2(j))
    dT <- decideTestsDGE(gT)
    s <- summary(dT)
    df[i,1] <- s[1]
    df[i,2] <- s[3]
    df[i,3] <- j
    df[i,4] <- y
    colnames(df) <- c('Down', 'Up', 'FC_cutoff', 'sample')
  }
  return(df)
}

PR <- as_tibble(oFC(PvsM, 'Pst'))
MR <- as_tibble(oFC(MvsM, 'Micro347'))
WR <- as_tibble(oFC(WvsM, 'Willi354'))
R <- bind_rows(PR, MR, WR)
R <- gather(R, UpDown, DEG, Down:Up)

R %>% group_by(sample, FC_cutoff) %>%
  summarise(DEG = sum(DEG)) %>%
  ggplot(aes(FC_cutoff, DEG))+
  geom_line(aes(col = factor(sample, level =c('Micro347', 'Willi354', 'Pst'))), size = 1)+
  theme_bw()+
  theme(legend.position = 'none',
        legend.title= element_blank(),
        aspect.ratio = 1.25,
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text = element_text(size = 12),
        panel.border = element_rect(size = 2, colour = 'black'),
        axis.title = element_text(colour = "black", size = 15),
        axis.ticks = element_blank(),
        # legend at the bottom 
        legend.text = element_text(colour = "black", size = 15),
        # plot title
        plot.title = element_text(colour = "black", size = 18))+
  scale_color_manual(values = c('#e1af00', '#3b9ab2', '#f21a00'))+
  labs(x='fold change threshold for TREAT', y='# of DEG')
# !! Testing for |logFC| > 1 by TREAT is not the same as selecting genes with |logFC| > 1
# Genes will need to exceed this threshold by some margin, depending on the data, before being declared statistically significant.
# -> ‘the fold-change below which we are definitely not interested in the gene’

## Differential expression relative to a fold-change threshold
DE.FC <- function(x){
  y <- glmTreat(fit, contrast=x, lfc=log2(1.3))
  topTags(y)
  is.de <- decideTestsDGE(y)
  return(is.de)
  plotMD(y, status=is.de, values=c(1,-1), col=c("red","blue"), legend='none',ylim=c(-8,6))
  box(lwd=2.5)
  return(is.de)
}
par(mfrow= c(1,3))
is.deM <- DE.FC(MvsM)
box(lwd=2.5)
is.deW <- DE.FC(WvsM)
box(lwd=2.5)
is.deP <- DE.FC(PvsM)
box(lwd=2.5)

## UpSet plot
DE <- function(x,y){
  z <- as.data.frame(x)
  z$ID <- rownames(z)
  colnames(z) <- c(y,'ID')
  z <- as_tibble(z)
  return(z)
}
deM <- DE(is.deM, "Micro.ud")
deW <- DE(is.deW, "Willi.ud")
deP <- DE(is.deP, "Pst.ud")

dePM <- merge(deP,deM,by ='ID')
de <- merge(dePM,deW,by='ID')
head(de)
de <- de %>% group_by(ID) %>% mutate(Pst = Pst.ud^2,
                                     Micro347 = Micro.ud^2,
                                     Willi354 = Willi.ud^2)
summary(de)
de <- de %>% group_by(ID) %>% mutate(sum.de = sum(c(Pst, Micro347,Willi354)),
                                     sum.ud = sum(c(Pst.ud, Micro.ud,Willi.ud))) %>%
  filter(sum.de > 0)
de <- de %>% group_by(ID) %>% mutate(updown = ifelse(sum.de == sum.ud, 'up', 
                                                     ifelse(sum.de^2 ==sum.ud^2, 'down', 'up_down')))
## !!!! get up and down data. below is only for upregulated as I filtered negatives away
de_upset <- as.data.frame(dplyr::select(de, -c(sum.ud, sum.de, Pst.ud, Micro.ud, Willi.ud)))
head(de_upset)
summary(de_upset)
upset(de_upset, order.by = "freq",
      text.scale = 1.5, point.size = 6, line.size = 1.5, mb.ratio = c(0.6, 0.4),
      mainbar.y.label = "DEGenes", sets.x.label = "DEGenes",
      sets.bar.color=c("#f21a00","#3b9ab2","#e1af00"),
      queries = list(
        list(
          query = elements, 
          params = list("updown", 'up'),
          color = "grey54", 
          active = T
        )
      )
)
# exported as 6x6 inch pdf
de_upset %>% filter(updown=='down')

Pv <- de_upset %>% filter(Pst %in% 1) %>% summarise(Pst = ID)
Wv <- de_upset %>% filter(Willi354 %in% 1) %>% summarise(Willi354 = ID)
Mv <- de_upset %>% filter(Micro347 %in% 1) %>% summarise(Micro347 = ID)
Vn <- list(Pst = Pv$Pst, Willi354 = Wv$Willi, Micro347 = Mv$Micro)

ggVennDiagram(Vn, label_alpha = 0, col = "grey24", label = 'count', edge_size=1)+
  scale_fill_gradient(low="white",high = "grey24")+
  scale_color_manual(values = c("#f21a00","#3b9ab2","#e1af00"))

#export
#write.table(de_upset, file= "DEGs_different_strains.txt", sep="\t", quote=F, row.names=FALSE)

## Heatmap clustering

logCPM <- cpm(y, prior.count=2, log=TRUE)  
#prior.count = average count to be added to each observation to avoid taking log of zero
#prior.count reduces the variability of logCPM values of low count genes.

#logCPM <- cpm(y, log=TRUE)

logCPMgeneid <- rownames(logCPM)
logCPM <- as_tibble(logCPM)
logCPM$ID <- logCPMgeneid
logCPM <- right_join(logCPM, de, by='ID')
logCPM <- logCPM %>% dplyr::select('ID', 'MOCK_1','MOCK_2','MOCK_3',
                            'MICRO_1', 'MICRO_2', 'MICRO_3',
                            'PST_1', 'PST_2', 'PST_3',
                            'WILLI_1', 'WILLI_2', 'WILLI_3')
logCPM <- as.data.frame(logCPM)
logCPMgeneid <- logCPM[,1]
logCPM <- logCPM[,2:13]
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
logCPM <- logCPM %>% rowwise() %>% mutate(m = mean(c(MOCK_1,MOCK_2,MOCK_3,PST_1,PST_2,PST_3,
                                                       WILLI_1,WILLI_2,WILLI_3,MICRO_1,MICRO_2,MICRO_3)))
logCPM <- logCPM %>% rowwise() %>% summarise(Mock1 = MOCK_1 - m, Mock2 = MOCK_2 - m, Mock3 = MOCK_3 -m,
                                             Pst1 = PST_1 - m, Pst2 = PST_2 - m, Pst3 = PST_3 -m,
                                             Willi1 = WILLI_1 - m, Willi2 = WILLI_2 - m, Willi3 = WILLI_3 -m,
                                             Micro1 = MICRO_1 - m, Micro2 = MICRO_2 - m, Micro3 = MICRO_3 -m,)

######
logCPM <- as.matrix(logCPM)
row.names(logCPM) <- logCPMgeneid
coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)
an <- HeatmapAnnotation(treat = colnames(logCPM),
                        show_legend = FALSE,
                        col = list(treat = c("Mock1" = "grey24","Mock2" = "grey24","Mock3" = "grey24",
                                             "Pst1" =  "#f21a00","Pst2" =  "#f21a00","Pst3" =  "#f21a00",
                                             "Willi1" ="#3b9ab2","Willi2" ="#3b9ab2","Willi3" ="#3b9ab2", 
                                             "Micro1" ="#e1af00","Micro2" ="#e1af00","Micro3" ="#e1af00")))
# kmeans clustering outside of complexHeatmap function
# determine optimal number of gene clusters via the elbow method
kmean_withinss <- function(k) {
  cluster <- kmeans(logCPM, k, nstart=50, iter.max=100)
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

kclus.row <- kmeans(logCPM, 10, nstart=50, iter.max=100)
kclus.row$cluster
split.row <- factor(paste0("\n", kclus.row$cluster), 
                    levels=c("\n2","\n1","\n4","\n10","\n3",
                             "\n8","\n5","\n7","\n9", "\n6"))

kclus.col <- kmeans(t(logCPM), 3, nstart=50, iter.max=100)
kclus.col$cluster
split.col <- factor(paste0("\n", kclus.col$cluster), 
                    levels=c("\n3","\n2","\n1"))

ht2 <- Heatmap(logCPM, column_split=split.col, cluster_column_slices = FALSE,
        row_split=split.row, cluster_row_slices = FALSE, show_row_names=FALSE,
        name = "logcpm", top_annotation = an, show_column_names=FALSE)



ht2 <- draw(ht2)
# exported as 5x7 inch pdf
r.dend <- row_dend(ht2)  #Extract row dendrogram
rcl.list <- row_order(ht2)  #Extract clusters (output is a list)
lapply(rcl.list, function(x) length(x))  #check/confirm size clusters
rcl.list[2]
# loop to extract genes for each cluster.
for (i in 1:length(row_order(ht2))){
  if (i == 1) {
    clu <- t(t(row.names(logCPM[row_order(ht2)[[i]],])))
    out <- cbind(clu, paste("cluster", i, sep=""))
    colnames(out) <- c("GeneID", "Cluster")
  } else {
    clu <- t(t(row.names(logCPM[row_order(ht2)[[i]],])))
    clu <- cbind(clu, paste("cluster", i, sep=""))
    out <- rbind(out, clu)
  }
}
#check
out 

#export
#write.table(out, file= "gene_cluster10.june.txt", sep="\t", quote=F, row.names=FALSE)

#load out
#setwd("/Volumes/GoogleDrive/My Drive/Transcriptome/RNAseq_Bai_clean")
#out <- read_csv('gene_clusters.csv')

out <- as_tibble(out)
lcpm.rn <- rownames(logCPM) 
lcpm <- as_tibble(logCPM)
lcpm$GeneID <- lcpm.rn
lcpm <- lcpm %>% pivot_longer(1:12,names_to = 'treat', values_to = 'logCPM') %>%
  mutate(treatment = if_else(treat %in% c("Mock1", "Mock2", "Mock3"), 'Mock',
                             if_else(treat %in% c("Pst1", "Pst2", "Pst3"), 'Pst',
                                     if_else(treat %in% c("Willi1", "Willi2", "Willi3"), 'Willi354',
                                             if_else(treat %in% c("Micro1", "Micro2", "Micro3"), 'Micro347',"x")))),
         replicate = rep(1:3, times = n()/3))  %>%
  group_by(treatment, GeneID) %>% mutate(mean.logCPM = mean(logCPM))

lcpm <- left_join(lcpm, out , by = "GeneID")

scaleFUN <- function(x) sprintf("%.1f", x)

lcpm %>% dplyr::select(-c(treat,logCPM,replicate)) %>% distinct() %>%
  filter(Cluster %in% "cluster10") %>% group_by(treatment) %>% mutate(m = mean(mean.logCPM), sd = sd(mean.logCPM)) %>%
  ggplot()+
  #geom_line(aes(factor(treatment,levels = c('Mock','Micro347', 'Willi354', 'Pst')), mean.logCPM, group = GeneID), alpha = .3)+
  #geom_line(aes(factor(treatment,levels = c('Mock','Micro347', 'Willi354', 'Pst')), m, group = GeneID), col = 'red', size = 1)+
  #geom_point(aes(factor(treatment,levels = c('Mock','Micro347', 'Willi354', 'Pst')), mean.logCPM), shape=21, size=3, alpha=.3, fill = "grey24")+
  geom_violin(aes(factor(treatment,levels = c('Mock','Micro347', 'Willi354', 'Pst')), mean.logCPM, fill = treatment), shape=21, size=1, alpha=.8)+
  geom_pointrange(aes(factor(treatment,levels = c('Mock','Micro347', 'Willi354', 'Pst')),m,ymax =m+sd,ymin=m-sd), size=1)+
  scale_fill_manual(values = c("#e1af00","grey24","#f21a00","#3b9ab2"))+
  scale_color_manual(values = c("#e1af00","grey24","#f21a00","#3b9ab2"))+
  theme_bw()+
  theme(legend.position = 'none',
        aspect.ratio = 1.3,
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        #panel.grid.minor.y = element_line(size = 2, colour = 'black'),
        legend.title = element_text(size = 15),
        axis.text = element_text(color = 'black',size = 25),
        panel.border = element_rect(size = 2, colour = 'black'),
        axis.title = element_text(colour = "black", size = 18),
        axis.ticks = element_blank(),
        # legend at the bottom 
        legend.text = element_text(colour = "black", size = 15),
        # plot title
        plot.title = element_text(colour = "black", size = 25),
        axis.text.x = element_blank()
        )+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  scale_y_continuous(labels=scaleFUN)+
  #labs(y = expression(centered~log[2](CPM)), x = '') 
  labs(y='',x='', title='cluster10') 



tuk.test <- function(x){
lcpm <- lcpm %>% dplyr::select(-c(treat,logCPM,replicate)) %>% distinct() %>%
  filter(Cluster %in% x) 
re <- TukeyHSD(aov(mean.logCPM ~ treatment,lcpm))
return(re)
}

tuk.test('cluster10')

lcpm %>% 
  filter(Cluster %in% "cluster1") %>% group_by(treat) %>% mutate(m = mean(logCPM), sd = sd(logCPM)) %>%
  ggplot()+
  #geom_line(aes(factor(treatment,levels = c('Mock','Micro347', 'Willi354', 'Pst')), mean.logCPM, group = GeneID), alpha = .3)+
  #geom_line(aes(factor(treatment,levels = c('Mock','Micro347', 'Willi354', 'Pst')), m, group = GeneID), col = 'red', size = 1)+
  #geom_point(aes(factor(treatment,levels = c('Mock','Micro347', 'Willi354', 'Pst')), mean.logCPM), shape=21, size=3, alpha=.3, fill = "grey24")+
  geom_violin(aes(factor(treat,levels = c('Mock1','Mock2','Mock3',
                                          'Micro1','Micro2','Micro3',
                                          'Willi1','Willi2','Willi3',
                                          'Pst1','Pst2','Pst3')), logCPM, fill = treatment), shape=21, alpha=.5)+
  geom_pointrange(aes(factor(treat,levels = c('Mock1','Mock2','Mock3',
                                                  'Micro1','Micro2','Micro3',
                                                  'Willi1','Willi2','Willi3',
                                                  'Pst1','Pst2','Pst3')),m,ymax =m+sd,ymin=m-sd))+
  scale_fill_manual(values = c("#e1af00","grey24","#f21a00","#3b9ab2"))+
  theme_bw()+
  theme(legend.position = 'right',
        aspect.ratio = 1,
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        #panel.grid.minor.y = element_line(size = 2, colour = 'black'),
        legend.title = element_text(size = 15),
        axis.text = element_text(color = 'black',size = 18),
        panel.border = element_rect(size = 2, colour = 'black'),
        axis.title = element_text(colour = "black", size = 18),
        axis.ticks = element_blank(),
        # legend at the bottom 
        legend.text = element_text(colour = "black", size = 15),
        # plot title
        plot.title = element_text(colour = "black", size = 15),
        axis.text.x = element_blank()
  )+
  scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  labs(y = expression(log[2](CPM)), x = 'inoculum density (log(cfu/ml))')+
  ylim(-2,2) 

tuk.test.all <- function(x){
  lcpm <- lcpm %>% 
    filter(Cluster %in% x) 
  re <- TukeyHSD(aov(logCPM ~ treat,lcpm))
  return(re)
}

tuk.test.all('cluster6')

#################################################################################################
### Functional enrichment (GO or KEGG) based on heatmap clusters
## KEGG
kegF <- function(x,y,i,j){
  z <- kegga(as.list(x), species="At")
  #z <- topKEGG(z, n=15)
  z$FC <- z$ID/(z$N/(27562+1883)*nrow(x)) #27562 protein genes and 1883 RNA genes in KEGG reference list
  z <- as_tibble(z[,c(1,i,j)]) %>% mutate(sample = rep(y, nrow(z)))
  return(z)
}

cl1 <- as_tibble(out) %>% filter(Cluster %in% 'cluster1') %>% summarise(ID = GeneID)
cl2 <- as_tibble(out) %>% filter(Cluster %in% 'cluster2') %>% summarise(ID = GeneID)
cl3 <- as_tibble(out) %>% filter(Cluster %in% 'cluster3') %>% summarise(ID = GeneID)
cl4 <- as_tibble(out) %>% filter(Cluster %in% 'cluster4') %>% summarise(ID = GeneID)
cl5 <- as_tibble(out) %>% filter(Cluster %in% 'cluster5') %>% summarise(ID = GeneID)
cl6 <- as_tibble(out) %>% filter(Cluster %in% 'cluster6') %>% summarise(ID = GeneID)
cl7 <- as_tibble(out) %>% filter(Cluster %in% 'cluster7') %>% summarise(ID = GeneID)
cl8 <- as_tibble(out) %>% filter(Cluster %in% 'cluster8') %>% summarise(ID = GeneID)
cl9 <- as_tibble(out) %>% filter(Cluster %in% 'cluster9') %>% summarise(ID = GeneID)
cl10 <- as_tibble(out) %>% filter(Cluster %in% 'cluster10') %>% summarise(ID = GeneID)

de.df <- de %>% dplyr::select(sum.de) %>% as.data.frame()
head(de.df)

cl1 <- merge(de.df, cl1, by = 'ID')
cl2 <- merge(de.df, cl2, by = 'ID')
cl3 <- merge(de.df, cl3, by = 'ID')
cl4 <- merge(de.df, cl4, by = 'ID')
cl5 <- merge(de.df, cl5, by = 'ID')
cl6 <- merge(de.df, cl6, by = 'ID')
cl7 <- merge(de.df, cl7, by = 'ID')
cl8 <- merge(de.df, cl8, by = 'ID')
cl9 <- merge(de.df, cl9, by = 'ID')
cl10 <- merge(de.df, cl10, by = 'ID')

kcl1 <- kegF(cl1, "cluster1", 5, 7)
kcl2 <- kegF(cl2, "cluster2", 5, 7)
kcl3 <- kegF(cl3, "cluster3", 5, 7)
kcl4 <- kegF(cl4, "cluster4", 5, 7)
kcl5 <- kegF(cl5, "cluster5", 5, 7)
kcl6 <- kegF(cl6, "cluster6", 5, 7)
kcl7 <- kegF(cl7, "cluster7", 5, 7)
kcl8 <- kegF(cl8, "cluster8", 5, 7)
kcl9 <- kegF(cl9, "cluster9", 5, 7)
kcl10 <- kegF(cl10, "cluster10", 5, 7)

kA <- bind_rows(kcl1, kcl2, kcl3, kcl4, kcl5, kcl6, kcl7, kcl8, kcl9, kcl10) %>%
  mutate(nlogPvalue = -1*log10(P.ID))

kA %>% filter(P.ID < .05) %>% #filter(FC > 5) %>%
  ggplot(aes(factor(sample, level = c(
    'cluster10', 'cluster9','cluster8', 'cluster7', 'cluster6',
    'cluster5', 'cluster4', 'cluster3', 'cluster2', 'cluster1')),Pathway))+
  geom_point(aes(fill=nlogPvalue, size=FC), shape=21)+
  scale_fill_gradient2(name = '-log10(p-value)', low = "white", mid = 'orange', high = "purple",
                       limits=c(0,max(kA$nlogPvalue)), midpoint = max(kA$nlogPvalue)/2,
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", lwd=3))+
  scale_size(name = 'fold enrichment',range = c(4,11))+
  theme_bw()+
  theme(legend.position = 'right',
        aspect.ratio = .35,
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.grid.minor.y = element_line(size = 2, colour = 'black'),
        legend.title = element_text(size = 15),
        axis.text = element_text(color = 'black',size = 16),
        panel.border = element_rect(size = 2, colour = 'black'),
        axis.title = element_text(colour = "black", size = 15),
        axis.ticks = element_blank(),
        # legend at the bottom 
        legend.text = element_text(colour = "black", size = 15),
        # plot title
        plot.title = element_text(colour = "black", size = 15),
        axis.text.x = element_text(angle = 66, vjust=1, hjust = 1))+
  labs(x='',y='')+
  coord_flip()
#exported as 12x8 inch pdf (p < 0.1)

## GO based on heatmap clusters
setwd("/Volumes/GoogleDrive/My Drive/Transcriptome/RNAseq_Bai_clean/bai.plots.june22")

go <- read_csv('go.cl10.june22.csv') %>% dplyr::select(c('GOn','GObp', 'clusterID', 'count', 'expected', 'fold_Enrichment', 'Pvalue'))
rev <- read_csv('revigo.out.csv')
go <- semi_join(go,rev,by='GOn')
go <- go %>% mutate(nlogPvalue = ifelse(-1*log10(Pvalue) < Inf, -1*log10(Pvalue), 0))

go %>% filter(Pvalue < .05) %>%
  filter(fold_Enrichment > 4) %>%
  ggplot(aes(factor(clusterID, level = c(
    'cluster10', 'cluster9','cluster8', 'cluster7', 'cluster6',
    'cluster5', 'cluster4', 'cluster3', 'cluster2', 'cluster1')),GObp))+
  geom_point(aes(fill=nlogPvalue, size=fold_Enrichment), shape=21)+
  scale_fill_gradient2(name = '-log10(p-value)', low = "white", mid = 'orange', high = "purple",
                       limits=c(0,max(go$nlogPvalue)), midpoint = max(go$nlogPvalue)/2,
                       guide = guide_colorbar(frame.colour = "black", ticks.colour = "black", lwd=3))+
  scale_size(name = 'fold enrichment',range = c(4,11))+
  theme_bw()+
  theme(legend.position = 'right',
        aspect.ratio = .35,
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.grid.minor.y = element_line(size = 2, colour = 'black'),
        legend.title = element_text(size = 15),
        axis.text = element_text(color = 'black',size = 16),
        panel.border = element_rect(size = 2, colour = 'black'),
        axis.title = element_text(colour = "black", size = 15),
        axis.ticks = element_blank(),
        # legend at the bottom 
        legend.text = element_text(colour = "black", size = 15),
        # plot title
        plot.title = element_text(colour = "black", size = 15),
        axis.text.x = element_text(angle = 66, vjust=1, hjust = 1))+
  labs(x='',y='')+
  coord_flip()
#exported as 10x5 inch pdf (p < 0.1)

go.h <- go %>% dplyr::select(c('GObp', 'Cluster', 'fold_Enrichment')) %>%
  pivot_wider(names_from = Cluster, values_from = fold_Enrichment)

go.h <- go.h %>% dplyr::select(-'GObp') %>% as.matrix()

Heatmap(go.h) # too many NAs

kclus.go <- kmeans(go.h, 10, nstart=50, iter.max=100)
