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
library(ComplexHeatmap) # ???
library(RColorBrewer)
library(circlize)
library(gplots)
library(S4Vectors)
library(org.At.tair.db)
library(org.Mm.eg.db)
library(UpSetR)
library("ggVennDiagram")

#global heatmap setting
ht_opt(legend_border = "black",
       heatmap_border = TRUE,
       annotation_border = TRUE)

# Data import
setwd("/Volumes/GoogleDrive/My Drive/Transcriptome/RNAseq_Willi_clean")
df <- read.delim('counts_willi_aR.txt', row.names ='Geneid')
head(df)

# define sample groups
group <- factor(rep(c('MOCK', 'W5', 'W6', 'W7', 'W8'), each = 3))
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

plotMDS(y, bg=rep(c("#8C510A", "#80CDC1", "#35978F", "#01665E", "#003C30"), each = 3), pch = 21, cex=2, bty="n")
box(lwd=2.5)
legend('topleft', bg="transparent", legend=c('Mock',expression(Willi354~10^5),expression(Willi354~10^6),expression(Willi354~10^7),expression(Willi354~10^8)), 
       pt.bg=c("#8C510A", "#80CDC1", "#35978F", "#01665E", "#003C30"),
       pch = 21, box.lty = 0, cex=1.2)

# Expression profiles of individual samples (Mean-Difference plots)
par(mfcol=c(3,5))
plotMD(y, column=15, ylim=c(-4,6)) #change column from 1 to 12
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
W8vsM <- makeContrasts(W8 - MOCK, levels=design)
W7vsM <- makeContrasts(W7 - MOCK, levels=design)
W6vsM <- makeContrasts(W6 - MOCK, levels=design)
W5vsM <- makeContrasts(W5 - MOCK, levels=design)

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
R8 <- as_tibble(oFC(W8vsM, 'Willi354 10^8'))
R7 <- as_tibble(oFC(W7vsM, 'Willi354 10^7'))
R6 <- as_tibble(oFC(W6vsM, 'Willi354 10^6'))
R5 <- as_tibble(oFC(W5vsM, 'Willi354 10^5'))
R <- bind_rows(R8, R7, R6, R5)
R <- gather(R, UpDown, DEG, Down:Up)


R %>% group_by(sample, FC_cutoff) %>%
  summarise(DEG = sum(DEG)) %>%
  ggplot(aes(FC_cutoff, DEG))+
  geom_line(aes(col = sample), size = 1)+
  theme_bw()+
  theme(legend.position = 'none',
        aspect.ratio = 1.25,
        #panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        panel.border = element_rect(size = 2, colour = 'black'),
        axis.title = element_text(colour = "black", size = 18),
        axis.ticks = element_blank(),
        # legend at the bottom 
        legend.text = element_text(colour = "black", size = 15),
        # plot title
        plot.title = element_text(colour = "black", size = 18)
  )+
  scale_color_manual(values = c("#80CDC1", "#35978F", "#01665E", "#003C30"))+
  labs(x='', y='# of DEG')
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
# genes that are significantly differential expressed are highlighted
par(mfrow= c(2,2))
is.de5 <- DE.FC(W5vsM)
box(lwd=2.5)
is.de6 <- DE.FC(W6vsM)
box(lwd=2.5)
is.de7 <- DE.FC(W7vsM)
box(lwd=2.5)
is.de8 <- DE.FC(W8vsM)
box(lwd=2.5)

################ Volcano plot
t <- glmTreat(fit, contrast=W8vsM)
t <- topTags(t, n=100000)$table
ggplot(data=t) + geom_point(aes(x=logFC,y=-log(FDR),color=logCPM)) +
  scale_colour_gradientn(colours=c("#000000" ,"#FF0000" ))+
  theme_bw()
is.de <- decideTestsDGE(t)
################

## UpSet plot
DE <- function(x,y){
  z <- as.data.frame(x)
  z$ID <- rownames(z)
  colnames(z) <- c(y,'ID')
  z <- as_tibble(z)
  return(z)
}
de5 <- DE(is.de5, "W5")
de6 <- DE(is.de6, "W6")
de7 <- DE(is.de7, "W7")
de8 <- DE(is.de8, "W8")

de56 <- merge(de5,de6,by ='ID')
de78 <- merge(de7,de8,by='ID')
de <- merge(de56,de78,by='ID')
head(de)
colnames(de) <- c('ID', 'W5updown', 'W6updown', 'W7updown', 'W8updown')
de <- de %>% group_by(ID) %>% mutate(W5 = W5updown^2,
                                     W6 = W6updown^2,
                                     W7 = W7updown^2,
                                     W8 = W8updown^2)
summary(de)
de <- de %>% group_by(ID) %>% mutate(sum.de = sum(c(W5, W6,W7,W8)),
                                     sum.ud = sum(c(W5updown, W6updown,W7updown,W8updown))) %>%
  filter(sum.de > 0)
de <- de %>% group_by(ID) %>% mutate(updown = ifelse(sum.de == sum.ud, 'up', 
                                                     ifelse(sum.de^2 ==sum.ud^2, 'down', 'up_down')))

de_upset <- as.data.frame(dplyr::select(de, -c(sum.ud, sum.de, W5updown, W6updown,W7updown,W8updown)))
head(de_upset)
summary(de_upset)
upset(de_upset, order.by = "freq",
      text.scale = 2, point.size = 4, line.size = 1.5, mb.ratio = c(0.65, 0.35),
      mainbar.y.label = "DEGenes", sets.x.label = "DEGenes",
      sets.bar.color=c("#003C30","#01665E","#35978F","#80CDC1"),
      queries = list(
        list(
          query = elements, 
          params = list("updown", 'up'),
          color = "grey54", 
          active = T
        )
      )
)
#exported as 4x6 inch pdf
de_upset %>% filter(updown=='down')
## Venn diagram
W5v <- de_upset %>% filter(W5 %in% 1) %>% summarise(W5 = ID)
W6v <- de_upset %>% filter(W6 %in% 1) %>% summarise(W6 = ID)
W7v <- de_upset %>% filter(W7 %in% 1) %>% summarise(W7 = ID)
W8v <- de_upset %>% filter(W8 %in% 1) %>% summarise(W8 = ID)
Vn <- list(W5 = W5v$W5, W6 = W6v$W6, W7 = W7v$W7, W8 = W8v$W8)

ggVennDiagram(Vn, label_alpha = 0, col = "grey24", label = 'count', edge_size=1)+
  scale_fill_gradient(low="white",high = "grey24")+
  scale_color_manual(values = c("#68abb8","#4f90a6","#3b738f","#2a5674"))

#export
#write.table(de_upset, file= "DEGs_different_density.txt", sep="\t", quote=F, row.names=FALSE)

## Heatmap clustering
logCPM <- cpm(y, prior.count=2, log=TRUE) #prior.count reduces the variability of logCPM values of low count genes.
#logCPM <- cpm(y, log=TRUE)

logCPMgeneid <- rownames(logCPM)
logCPM <- as_tibble(logCPM)
logCPM$ID <- logCPMgeneid
logCPM <- right_join(logCPM, de, by='ID')
logCPM <- logCPM %>% dplyr::select('ID', 'MOCK_1','MOCK_2','MOCK_3',
                                   'W5_1', 'W5_2', 'W5_3',
                                   'W6_1', 'W6_2', 'W6_3',
                                   'W7_1', 'W7_2', 'W7_3',
                                   'W8_1', 'W8_2', 'W8_3')
logCPM <- as.data.frame(logCPM)
logCPMgeneid <- logCPM[,1]
logCPM <- logCPM[,2:16]
row.names(logCPM) <- logCPMgeneid
#scale each row/gene to have mean zero and standard deviation one
#logCPM <- t(scale(t(logCPM)))
#or
#scale the rows by subtracting the row median (see Vogel2016)

logCPM <- as_tibble(logCPM)
logCPM <- logCPM %>% rowwise() %>% mutate(m = mean(c(MOCK_1,MOCK_2,MOCK_3,
                                                       W5_1,W5_2,W5_3,
                                                       W6_1,W6_2,W6_3,
                                                       W7_1,W7_2,W7_3,
                                                       W8_1,W8_2,W8_3)))
logCPM <- logCPM %>% rowwise() %>% summarise(Mock1 = MOCK_1 - m, Mock2 = MOCK_2 - m, Mock3 = MOCK_3 -m,
                                             W5.1 = W5_1 - m, W5.2 = W5_2 - m, W5.3 = W5_3 - m,
                                             W6.1 = W6_1 - m, W6.2 = W6_2 - m, W6.3 = W6_3 - m,
                                             W7.1 = W7_1 - m, W7.2 = W7_2 - m, W7.3 = W7_3 - m,
                                             W8.1 = W8_1 - m, W8.2 = W8_2 - m, W8.3 = W8_3 - m)
####
logCPM <- as.matrix(logCPM)
row.names(logCPM) <- logCPMgeneid
coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)
an <- HeatmapAnnotation(treat = colnames(logCPM),
                        show_legend = FALSE,
                        border = TRUE,
                        col = list(treat = c("Mock1" = "grey24","Mock2" = "grey24","Mock3" = "grey24",
                                             "W5.1" =  "#68abb8","W5.2" =  "#68abb8","W5.3" =  "#68abb8",
                                             "W6.1" =  "#4f90a6","W6.2" =  "#4f90a6","W6.3" =  "#4f90a6", 
                                             "W7.1" =  "#3b738f","W7.2" =  "#3b738f","W7.3" =  "#3b738f",
                                             "W8.1" =  "#2a5674","W8.2" =  "#2a5674","W8.3" =  "#2a5674")))
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
        axis.title = element_text(colour = "black", size = 18),
        axis.ticks = element_blank(),
        # legend at the bottom 
        legend.text = element_text(colour = "black", size = 15),
        # plot title
        plot.title = element_text(colour = "black", size = 18),
       )+
  labs(x='# of clusters', y='total within sum of squares')
## exported as 4x6 inch portrait pdf
kclus.row <- kmeans(logCPM, 10, nstart=50, iter.max=100)
kclus.row$cluster
split.row <- factor(paste0("\n", kclus.row$cluster), 
                    levels=c("\n5","\n3","\n1","\n9", "\n7",
                             "\n8","\n6","\n4","\n2","\n10"))

kclus.col <- kmeans(t(logCPM), 5, nstart=50, iter.max=100)
kclus.col$cluster
split.col <- factor(paste0("\n", kclus.col$cluster), 
                    levels=c("\n2","\n3","\n1","\n5","\n4"))

ht2 <- Heatmap(logCPM, cluster_columns = FALSE,
               row_split=split.row, cluster_row_slices = FALSE, show_row_names=FALSE,
               name = "logcpm", top_annotation = an, show_column_names=FALSE)
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
#setwd("/Volumes/GoogleDrive/My Drive/Transcriptome/RNAseq_Willi_clean/willi.june")
#write.table(out, file= "gene_clusters10.txt", sep="\t", quote=F, row.names=FALSE)
#load 'out'
setwd("/Volumes/GoogleDrive/My Drive/Transcriptome/RNAseq_Willi_clean/willi.june")
out <- read_delim('gene_clusters10.txt')

out <- as_tibble(out)
lcpm.rn <- rownames(logCPM) 
lcpm <- as_tibble(logCPM)
lcpm$GeneID <- lcpm.rn
lcpm <- lcpm %>% pivot_longer(1:15,names_to = 'treat', values_to = 'logCPM') %>%
  mutate(treatment = if_else(treat %in% c("Mock1", "Mock2", "Mock3"), 'Mock',
                             if_else(treat %in% c("W5.1", "W5.2", "W5.3"), 'W5',
                                     if_else(treat %in% c("W6.1", "W6.2", "W6.3"), 'W6',
                                             if_else(treat %in% c("W7.1", "W7.2", "W7.3"), 'W7',
                                                     if_else(treat %in% c("W8.1", "W8.2", "W8.3"), 'W8',"x"))))),
         replicate = rep(1:3, times = n()/3)) %>%
  group_by(treatment, GeneID) %>% mutate(mean.logCPM = mean(logCPM),
                                         log.inoc = if_else(treatment == 'Mock', 0,
                                                            if_else(treatment == 'W5', 5,
                                                                    if_else(treatment == 'W6', 6,
                                                                            if_else(treatment == 'W7', 7,
                                                                                    if_else(treatment == 'W8', 8,9))))))

lcpm <- left_join(lcpm, out , by = "GeneID")

scaleFUN <- function(x) sprintf("%.1f", x)

lcpm %>% dplyr::select(-c(treat,logCPM,replicate)) %>% distinct() %>%
  filter(Cluster %in% "cluster2") %>% group_by(log.inoc) %>%
  mutate(m = mean(mean.logCPM), se = sd(mean.logCPM), n=n()) %>%
  ggplot()+
  geom_line(aes(log.inoc, mean.logCPM, group = GeneID), alpha = .3)+
  geom_line(aes(log.inoc, m), col = 'red2', size = 2)+
  geom_pointrange(aes(log.inoc, m, ymin=m-se, ymax=m+se), fill = 'white', shape=21, size = 1.5, col ='red3')+
  #geom_point(aes(log.inoc, mean.logCPM), shape=21, size=2, alpha=.3, fill = "grey24")+
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
        plot.title = element_text(colour = "black", size = 25)
  )+
  #scale_x_discrete(guide = guide_axis(n.dodge = 2))+
  scale_y_continuous(labels=scaleFUN)+
  #labs(y = expression(centered~log[2](CPM)), x = '') 
  labs(y='',x='', title='cluster2') 
  
#################################################################################################
## KEGG analysis
kegF <- function(x,y,i,j){
  z <- kegga(as.list(x), species="At")
  #z <- topKEGG(z, n=15)
  z$FC <- z$ID/(z$N/(27562+1883)*nrow(x)) #27562 protein genes and 1883 RNA genes in KEGG reference list
  z <- as_tibble(z[,c(1,3,i,j)]) %>% mutate(sample = rep(y, nrow(z)))
  return(z)
  kegW8$FC <- kegW8$ID/(kegW8$N/55865*1204)
}

### Functional enrichment (GO or KEGG) based on heatmap clusters
## KEGG based on heatmap clusters
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

kA %>% filter(ID > 1) %>%
  filter(P.ID < .01) %>% #filter(FC > 5) %>%
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
setwd("/Volumes/GoogleDrive/My Drive/Transcriptome/RNAseq_Willi_clean/willi.june") # new google drive instance

kegg.gene <- getGeneKEGGLinks(species.KEGG = "ath", convert = FALSE)
kegg.gene <- as_tibble(kegg.gene) %>% dplyr::rename(KEGGid = PathwayID, geneID = GeneID)

kegg.cl <- read_csv('KEGGwilli.csv')

# keep only KEGG pathways that were enriched in gene set

kegg.sub <- semi_join(kegg.gene, kegg.cl, by = 'KEGGid')

# which genes are DE in the enriched pathways?
de.cl <- lcpm %>% ungroup() %>%
  dplyr::select(GeneID, Cluster) %>% dplyr::rename(geneID = GeneID, clusterID = Cluster) %>%
  distinct()

de.kegg <- inner_join(de.cl, kegg.sub, by = 'geneID') 
k.cl <- left_join(kegg.cl, de.kegg, by = 'KEGGid')
# which genes are DE in the enriched pathway/cluster combination?
k.cl <- left_join(kegg.cl, de.kegg, by = c('KEGGid','clusterID'))
#write.table(k.cl, file= "kegg.geneid.june.txt", sep="\t", quote=F, row.names=FALSE)

gd <- read.delim('gene_descriptions.txt')
gdc <- left_join(out,gd,by='GeneID')
#write.table(gdc, file='gene.description.cluster.txt', sep="\t", quote=F, row.names=FALSE)
k.cl <- k.cl %>% dplyr::rename(GeneID = geneID)
k.cl.d <- left_join(k.cl,gd,by='GeneID')
#write.table(k.cl.d, file='"kegg.geneid.june.txt"', sep="\t", quote=F, row.names=FALSE)

ePlot(kA, kA$Pathway, kA$FC, kA$P.ID, kA$sample,kA$nlogPvalue,
      "KEGG pathway enrichment", .05)

## GO based on heatmap clusters
setwd("/Volumes/GoogleDrive/My Drive/Transcriptome/RNAseq_Willi_clean/GOcluster10")

go1 <- read_csv('GOcl1.csv') %>% dplyr::select(c('GObp', 'fold_Enrichment', 'Pvalue'))
go1$sample <- rep('cluster1', nrow(go1))
go2 <- read_csv('GOcl2.csv') %>% dplyr::select(c('GObp', 'fold_Enrichment', 'Pvalue'))
go2$sample <- rep('cluster2', nrow(go2))
go3 <- read_csv('GOcl3.csv') %>% dplyr::select(c('GObp', 'fold_Enrichment', 'Pvalue'))
go3$sample <- rep('cluster3', nrow(go3))
go4 <- read_csv('GOcl4.csv') %>% dplyr::select(c('GObp', 'fold_Enrichment', 'Pvalue'))
go4$sample <- rep('cluster4', nrow(go4))
go5 <- read_csv('GOcl5.csv') %>% dplyr::select(c('GObp', 'fold_Enrichment', 'Pvalue'))
go5$sample <- rep('cluster5', nrow(go5))
go6 <- read_csv('GOcl6.csv') %>% dplyr::select(c('GObp', 'fold_Enrichment', 'Pvalue'))
go6$sample <- rep('cluster6', nrow(go6))
go7 <- read_csv('GOcl7.csv') %>% dplyr::select(c('GObp', 'fold_Enrichment', 'Pvalue'))
go7$sample <- rep('cluster7', nrow(go7))
go8 <- read_csv('GOcl8.csv') %>% dplyr::select(c('GObp', 'fold_Enrichment', 'Pvalue'))
go8$sample <- rep('cluster8', nrow(go8))
go10 <- read_csv('GOcl10.csv') %>% dplyr::select(c('GObp', 'fold_Enrichment', 'Pvalue'))
go10$sample <- rep('cluster10', nrow(go10))

goA <- bind_rows(go1, go2, go3, go4, go5, go6, go7, go8, go10) %>% 
  mutate(nlogPvalue = ifelse(-1*log10(Pvalue) < Inf, -1*log10(Pvalue), 0)) %>%
  filter(fold_Enrichment > 4)

ePlot(goA, goA$GObp, goA$fold_Enrichment, goA$Pvalue, goA$sample, goA$nlogPvalue,
      "GO term enrichment", .0001)

### Functional enrichment based on groups from UpSet plot
w8sp <- de %>% filter(W8 %in% 1 & sum.de %in% 1) %>% dplyr::select(c('ID', 'W8'))
w8nsp <- de %>% filter(W8 %in% 1 & sum.de != 1) %>% dplyr::select(c('ID', 'W8'))

## KEGG
k8sp <- kegF(w8sp, "specific to Willi 10^8 cfu/ml", 5, 7)
k8nsp <- kegF(w8nsp, "non-specific to Willi 10^8 cfu/ml", 5, 7)

kA <- bind_rows(k8sp, k8nsp) %>% mutate(nlogPvalue = -1*log10(P.ID))

ePlot(kA, kA$Pathway, kA$FC, kA$P.ID, kA$sample,kA$nlogPvalue,
      "KEGG pathway enrichment", .05)

## GO
setwd("/Volumes/GoogleDrive/My Drive/Transcriptome/RNAseq_Willi_clean/GOupset")

write_csv(w8sp, file = 'w8sp.csv')
write_csv(w8nsp, file = 'w8nsp.csv')

go8sp <- read_csv('GO8sp.csv') %>% dplyr::select(c('GObp', 'fold_Enrichment', 'Pvalue'))
go8sp$sample <- rep('specific to Willi(10^8)', nrow(go8sp))
go8nsp <- read_csv('GO8nsp.csv') %>% dplyr::select(c('GObp', 'fold_Enrichment', 'Pvalue'))
go8nsp$sample <- rep('nonspecific to Willi(10^8)', nrow(go8nsp))

goA <- bind_rows(go8sp, go8nsp) %>% 
  mutate(nlogPvalue = ifelse(-1*log10(Pvalue) < Inf, -1*log10(Pvalue), 0)) %>%
  filter(fold_Enrichment > 4)

ePlot(goA, goA$GObp, goA$fold_Enrichment, goA$Pvalue, goA$sample, goA$nlogPvalue,
      "GO term enrichment", .0001)


#### add gene descriptions to the Willi 1e8 specific and non-specific GOterms 

setwd("/Volumes/GoogleDrive/My Drive/Transcriptome/RNAseq_Willi_clean/willi.june")

go.g <- read_csv('GOWilligenes.csv')

gd <- read.delim('gene_descriptions_willi.txt')
go.g.d <- left_join(go.g,gd,by='GeneID')

write_csv(go.g.d, 'GO_gene_decriptions.csv')

#### UpSet plot of enriched GO genes (> 6 fold increase)

go.up <- go.g.d %>% summarise(GOterm = GOterm,g.name=if_else(Primary.Gene.Symbol == '',GeneID,Primary.Gene.Symbol),presence = 1)
go.up <- go.up %>% pivot_wider(names_from = GOterm, values_from = presence)
go.up[is.na(go.up)] <- 0
go.up <- as.data.frame(go.up)
upset(go.up, order.by = "freq", nsets = 12, group.by = "sets",
      text.scale = 2, point.size = 4, line.size = 1.5, mb.ratio = c(0.65, 0.35),
      mainbar.y.label = "DEGenes", sets.x.label = "DEGenes",
  sets.bar.color=rev(c('#003C30','#003C30','#003C30','#003C30','#003C30','#E1AF00','#E1AF00','#E1AF00','#E1AF00','#003C30','#003C30','#E1AF00')))

upset(go.up, order.by = "freq", nsets = 12,
      text.scale = 2, point.size = 4, line.size = 1.5, mb.ratio = c(0.65, 0.35),
      mainbar.y.label = "DEGenes", sets.x.label = "DEGenes",
      sets = c('response to nitrogen compound', 'response to organonitrogen compound','response to light intensity',
               'response to oxygen levels', 'cellular response to hypoxia',
               'systemic acquired resistance', 'response to molecule of bacterial origin','regulation of salicylic acid mediated signaling pathway',
               'response to insect', 'cell wall thickening', 'indole glucosinolate metabolic process',
               'defense response by callose deposition in cell wall'),
      keep.order = TRUE,
      sets.bar.color=rev(c('#003C30','#003C30','#003C30','#003C30','#003C30','#003C30','#003C30','#E1AF00','#E1AF00','#E1AF00','#E1AF00','#E1AF00')))

#### Heatmap per GOterm
go.g.d$symbol <- sub("\\).*", "", sub(".*\\(", "", go.g.d$Primary.Gene.Symbol)) 

ht.go <- function(x){
xx <- go.g.d %>% filter(GOterm %in% x) %>%
  summarise(g.name=if_else(Primary.Gene.Symbol == '',GeneID,symbol), GeneID = GeneID)
yy <- as_tibble(logCPM)
yy$GeneID <- rownames(logCPM)
zz <- right_join(yy,xx, by='GeneID')
zz.rn <- zz$g.name
zz <- as.matrix(zz[,1:15])
rownames(zz) <- zz.rn
ht <- Heatmap(zz, cluster_columns=FALSE, #show_heatmap_legend = FALSE,
        show_row_names=TRUE,
        name = "logcpm", top_annotation = an, show_column_names=FALSE,
        column_title = x,
        col=colorRamp2(c(-4, 0, 4), c("blue", "white", "red")),
width = ncol(zz)*unit(5, "mm"), 
height = nrow(zz)*unit(3.5, "mm"))
return(ht)
}
draw(ht.go('response to nitrogen compound'))
draw(ht.go('response to organonitrogen compound'))
draw(ht.go('response to light intensity'))
draw(ht.go('response to oxygen levels'))
draw(ht.go('cellular response to hypoxia'))
draw(ht.go('systemic acquired resistance'))
draw(ht.go('response to molecule of bacterial origin'))
draw(ht.go('regulation of salicylic acid mediated signaling pathway'))
draw(ht.go('response to insect'))
draw(ht.go('cell wall thickening'))
draw(ht.go('indole glucosinolate metabolic process'))
draw(ht.go('defense response by callose deposition in cell wall'))


